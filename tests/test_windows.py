from pathlib import Path
import h5py
import numpy as np
import pytest
import tdt
from clean_dbs.flags import find_saturation, Saturation
from clean_dbs.save_h5 import save_h5, open_h5, StreamWriter, write_rest
from clean_dbs.read_bv import read_bv, bv_windows
from clean_dbs.read_tdt import read_tdt, tdt_windows
from clean_dbs.clean_block import clean_block


def whole_saturation(stream, min_ms=1.0):
    """
    Input: a stream with data, min_ms
    Returns: find_saturation as it was written before the windowed version (whole array at once), to check against
    """
    data, fs = stream["data"], float(stream["fs"])
    min_samples = max(2, int(round(min_ms / 1000 * fs)))
    channel, start, end = [], [], []
    for row, number in enumerate(stream["channels"]):
        x = data[row]
        if x.size == 0 or x.max() == x.min():
            continue
        for extreme in (x.max(), x.min()):
            hits = np.flatnonzero(x == extreme)
            if hits.size < min_samples:
                continue
            breaks = np.flatnonzero(np.diff(hits) > 1)
            run_starts = np.concatenate(([hits[0]], hits[breaks + 1]))
            run_ends = np.concatenate((hits[breaks], [hits[-1]]))
            keep = (run_ends - run_starts + 1) >= min_samples
            channel += [int(number)] * int(keep.sum())
            start += list(stream["start_time"] + run_starts[keep] / fs)
            end += list(stream["start_time"] + (run_ends[keep] + 1) / fs)
    return {"channel":np.asarray(channel, dtype=np.int64), "start":np.asarray(start, dtype=np.float64), "end":np.asarray(end, dtype=np.float64)}


def windowed_saturation(stream, sizes, min_ms=1.0):
    #feeds the data in windows of the given sizes (cycled), like clean_block does
    runs, a, k = Saturation(stream, min_ms), 0, 0
    n = stream["data"].shape[1]
    while a < n:
        b = min(n, a + sizes[k % len(sizes)])
        runs.add(stream["data"][:, a:b])
        a, k = b, k + 1
    return runs.result()


def assert_same(x, y):
    assert x.keys() == y.keys()
    for key in x:
        assert x[key].dtype == y[key].dtype
        np.testing.assert_array_equal(x[key], y[key])


def test_saturation_run_across_window_edge():
    fs = 10000.0 #min_ms=1 -> runs need 10 samples
    x = np.linspace(-0.5, 0.5, 1001, dtype=np.float32) #odd length
    x[95:105] = 2.0 #10 samples at the max, split 5 + 5 by a window edge at 100: only counts once joined
    x[300:304] = 2.0 #4 samples, too short even though it touches the edge at 300
    x[600:700] = -2.0 #at the min across two edges
    x[985:] = 2.0 #16 samples running into the very end
    stream = {"data":x[None, :], "fs":fs, "start_time":0.25, "channels":np.array([5])}
    expected = whole_saturation(stream)
    assert len(expected["start"]) == 3
    for sizes in ([100], [1], [7, 3, 300], [1001], [100, 200, 300, 296, 105]):
        assert_same(windowed_saturation(stream, sizes), expected)
    assert_same(find_saturation(stream), expected)


def test_saturation_new_max_in_later_window():
    #the first windows' max (and its runs) stop counting once a later window goes higher
    x = np.zeros(400, dtype=np.float32)
    x[10:40] = 1.0
    x[200:230] = 3.0
    x[390:400] = 3.0
    x[50:80] = -1.0
    stream = {"data":x[None, :], "fs":1000.0, "start_time":0.0, "channels":[1]}
    for sizes in ([100], [33], [1, 2, 3]):
        assert_same(windowed_saturation(stream, sizes, min_ms=5.0), whole_saturation(stream, min_ms=5.0))


def test_saturation_random_against_whole():
    rng = np.random.default_rng(0)
    for trial in range(30):
        n_ch, n = 3, int(rng.integers(50, 3000))
        data = rng.integers(-3, 4, size=(n_ch, n)).astype(np.float32) #few levels, so lots of runs at the max and min
        data[1] = 0.5 #flat channel: never flagged
        if trial % 5 == 0:
            data[2, int(rng.integers(0, n))] = np.nan #a NaN anywhere: never flagged, like the whole-array version
        stream = {"data":data, "fs":2000.0, "start_time":float(rng.random()), "channels":np.array([4, 9, 11])}
        sizes = list(rng.integers(1, 400, size=5))
        assert_same(windowed_saturation(stream, sizes), whole_saturation(stream))


def fake_streams(n):
    rng = np.random.default_rng(1)
    return {"ecos":{"data":rng.standard_normal((3, n)).astype(np.float32), "fs":48828.125, "start_time":0.0, "channels":np.array([1, 2, 3]),
                    "device":"ecogdbs", "units":"V", "type":"lfp", "subtype":"ecog+lead", "channel_names":["a", "b", "c"]},
            "mic_":{"data":rng.standard_normal((1, n // 2 + 1)).astype(np.float32), "fs":24414.0625, "start_time":0.0, "channels":np.array([1]),
                    "device":"microphone", "units":"V", "type":"audio", "subtype":""},
            "empty":{"data":np.zeros((2, 0), dtype=np.float32), "fs":1.0, "start_time":0.0, "channels":np.array([1, 2])}}


def h5_contents(path):
    #every dataset (values and storage layout) and attribute, except when the file was made
    out = {}
    def visit(name, obj):
        out[name + "@attrs"] = {k:repr(np.asarray(v).tolist()) for k, v in obj.attrs.items() if k != "created"} #repr so NaN matches NaN
        if isinstance(obj, h5py.Dataset):
            values = obj[()]
            values = values.tolist() if values.dtype.kind == "O" else values.tobytes() #bytes so NaN matches NaN
            out[name] = (values, obj.shape, obj.dtype.str, obj.chunks, obj.compression, obj.shuffle)
    with h5py.File(path, "r") as f:
        f.visititems(visit)
        out["/@attrs"] = {k:repr(np.asarray(v).tolist()) for k, v in f.attrs.items() if k != "created"}
    return out


def test_stream_writer_matches_save_h5(tmp_path):
    n = 2 * 65536 + 13 #odd length, last hdf5 chunk short
    streams = fake_streams(n)
    result = {"streams":streams, "epocs":{}, "scalars":{}, "info":{"block":"b"}, "flags":{"dnu":np.int8(0)}}
    source = tmp_path / "raw"
    source.mkdir()
    save_h5(result, tmp_path / "whole.h5", source)
    meta = {name:{k:v for k, v in s.items() if k != "data"} | {"n_samples":s["data"].shape[1]} for name, s in streams.items()}
    for sizes in ([65536], [1000, 70000, 5], [n]):
        out = tmp_path / f"windows_{len(sizes)}_{sizes[0]}.h5"
        with open_h5(out) as f:
            writer = StreamWriter(f, meta)
            a = dict.fromkeys(streams, 0)
            for k in range(1000):
                window = {}
                for name, s in streams.items():
                    b = min(s["data"].shape[1], a[name] + sizes[k % len(sizes)])
                    if b > a[name]:
                        window[name], a[name] = (a[name], s["data"][:, a[name]:b]), b
                if not window:
                    break
                writer.write(window)
            writer.close()
            write_rest(f, {**result, "streams":meta}, source)
        assert h5_contents(out) == h5_contents(tmp_path / "whole.h5")


def test_stream_writer_catches_gaps(tmp_path):
    meta = {"ecos":{"n_samples":100, "fs":1.0, "start_time":0.0, "channels":[1]}}
    with pytest.raises(ValueError, match="back to back"):
        with open_h5(tmp_path / "a.h5") as f:
            writer = StreamWriter(f, meta)
            writer.write({"ecos":(0, np.zeros((1, 40), dtype=np.float32))})
            writer.write({"ecos":(41, np.zeros((1, 40), dtype=np.float32))}) #skipped sample 40
    assert not (tmp_path / "a.h5").exists() and (tmp_path / "a.h5.part").exists()
    with pytest.raises(ValueError, match="samples arrived"):
        with open_h5(tmp_path / "b.h5") as f:
            writer = StreamWriter(f, meta)
            writer.write({"ecos":(0, np.zeros((1, 99), dtype=np.float32))}) #one short
            writer.close()


def write_bv(folder, data, fs=10000.0):
    """
    Input: folder (made like <patient>/<stage> so clean_block can name the output), data (channels x samples, microvolts), fs
    Returns: path of a small brainvision recording (.vhdr/.vmrk/.eeg, float32, multiplexed) MNE can read
    """
    folder.mkdir(parents=True, exist_ok=True)
    n_ch = data.shape[0]
    channels = "".join(f"Ch{k + 1}=DBS{k + 1},,0.1,µV\n" for k in range(n_ch))
    (folder / "rec.vhdr").write_text("Brain Vision Data Exchange Header File Version 1.0\n[Common Infos]\nCodepage=UTF-8\nDataFile=rec.eeg\n"
                                     f"MarkerFile=rec.vmrk\nDataFormat=BINARY\nDataOrientation=MULTIPLEXED\nNumberOfChannels={n_ch}\n"
                                     f"SamplingInterval={1e6 / fs:g}\n[Binary Infos]\nBinaryFormat=IEEE_FLOAT_32\n[Channel Infos]\n" + channels,
                                     encoding="utf-8")
    (folder / "rec.vmrk").write_text("Brain Vision Data Exchange Marker File, Version 1.0\n[Common Infos]\nCodepage=UTF-8\nDataFile=rec.eeg\n"
                                     "[Marker Infos]\nMk1=New Segment,,1,1,0,20240101120000000000\nMk2=Stimulus,S 13,5000,1,0\n", encoding="utf-8")
    (data.T.astype(np.float32) / np.float32(0.1)).tofile(folder / "rec.eeg") #stored values x resolution 0.1 = microvolts
    return folder / "rec.vhdr"


def test_bv_windows_match_read_bv(tmp_path):
    rng = np.random.default_rng(2)
    data = rng.standard_normal((3, 140001)) * 50 #odd sample count, not a whole number of windows
    data[1, 65530:65560] = 400 #saturated across the first window edge
    vhdr = write_bv(tmp_path / "ug9999" / "stage1 surgery", data)
    whole = read_bv(vhdr)
    meta, windows = bv_windows(vhdr, window_mb=0.5) #small enough for 3 windows
    joined = np.concatenate([w["lead"][1] for w in windows], axis=1)
    assert meta["streams"]["lead"]["n_samples"] == joined.shape[1] == 140001
    np.testing.assert_array_equal(joined, whole["streams"]["lead"]["data"])
    #the whole cleaned file is the same for any window size
    a = clean_block(vhdr, tmp_path / "small", window_mb=0.5)
    b = clean_block(vhdr, tmp_path / "big", window_mb=1000)
    assert h5_contents(a) == h5_contents(b)
    with h5py.File(a, "r") as f:
        assert f["flags/saturation_lead/channel"][:].tolist() == [2]


def write_tdt(block, streams, start=1.7e9):
    """
    Input: block folder to make, streams {name: (fs, channels x samples float32, samples per packet)}
    Returns: block, a minimal tdt block (.tsq index + .tev samples + Notes.txt) that tdt.read_block reads
    Packets go in time order with channels interleaved, like Synapse writes them, and get timestamps with the same
    ~1e-7 s of float rounding real blocks have (absolute time in s since 1970).
    """
    block.mkdir(parents=True)
    rows, tev = [], bytearray()
    def row(size, type_, code, chan, ts, offset, fs):
        r = np.zeros(1, dtype=[("size", "<i4"), ("type", "<i4"), ("code", "<u4"), ("chan", "<u2"), ("sort", "<u2"),
                               ("ts", "<f8"), ("offset", "<u8"), ("format", "<i4"), ("fs", "<f4")])
        r[0] = (size, type_, code, chan, 0, ts, offset, 0, fs)
        return r.tobytes()
    packets = []
    for name, (fs, data, per) in streams.items():
        code = int.from_bytes(name.encode(), "little")
        for p in range(data.shape[1] // per):
            for c in range(data.shape[0], 0, -1):
                packets.append((start + p * per / fs, name, code, c, fs, data[c - 1, p * per:(p + 1) * per], per))
    packets.sort(key=lambda q: q[0])
    end = max(q[0] for q in packets) + 1
    rows.append(row(0, 0, 0, 0, 0.0, 0, 0)) #file header, not read
    rows.append(row(10, tdt.EVTYPE_MARK, tdt.EVMARK_STARTBLOCK, 0, start, 0, 0))
    for ts, name, code, c, fs, samples, per in packets:
        rows.append(row(10 + per, tdt.EVTYPE_STREAM, code, c, ts, len(tev), fs))
        tev += samples.astype(np.float32).tobytes()
    rows.append(row(10, tdt.EVTYPE_MARK, tdt.EVMARK_STOPBLOCK, 0, end, 0, 0))
    (block / "b.tsq").write_bytes(b"".join(rows))
    (block / "b.tev").write_bytes(bytes(tev))
    (block / "b.tnt").write_text("1\n") #empty notes file, so tdt doesn't warn
    (block / "Notes.txt").write_text("Experiment: test\nSubject: uh9999\nUser: lab\nStart: 10:00:00am 01/01/2024\nStop: 10:01:00am 01/01/2024\n")
    return block


def test_tdt_windows_match_read_tdt(tmp_path):
    rng = np.random.default_rng(3)
    fs = 48828.125 #a tdt rate: window edges never land on whole samples
    ecos = rng.standard_normal((4, 1024 * 301)).astype(np.float32) #301 packets of 1024
    ecos[2, 65536 - 60:65536 + 60] = 9.0 #saturated across a window edge (1 ms = 49 samples at this rate)
    mic = rng.standard_normal((1, 512 * 151)).astype(np.float32)
    block = write_tdt(tmp_path / "uh9999" / "stage1 surgery" / "block1", {"ecos":(fs, ecos, 1024), "mic_":(fs / 2, mic, 512)})
    whole = read_tdt(block)
    np.testing.assert_array_equal(whole["streams"]["ecos"]["data"], ecos) #the made-up block reads back as written
    for mb in (0.3, 3.0, 1000):
        meta, windows = tdt_windows(block, window_mb=mb)
        joined = {}
        for window in windows:
            for name, (a, x) in window.items():
                assert a == sum(y.shape[1] for y in joined.get(name, [])) #back to back
                joined.setdefault(name, []).append(x)
        for name, s in whole["streams"].items():
            np.testing.assert_array_equal(np.concatenate(joined[name], axis=1), s["data"])
            assert meta["streams"][name]["n_samples"] == s["data"].shape[1]
            assert meta["streams"][name]["start_time"] == s["start_time"]
            np.testing.assert_array_equal(meta["streams"][name]["channels"], s["channels"])
    a = clean_block(block, tmp_path / "small", window_mb=0.3)
    b = clean_block(block, tmp_path / "big", window_mb=1000)
    assert h5_contents(a) == h5_contents(b)
    with h5py.File(a, "r") as f:
        assert f["flags/saturation_ecos/channel"][:].tolist() == [3]

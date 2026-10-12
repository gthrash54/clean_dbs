from pathlib import Path
import numpy as np
from clean_dbs.flags import find_saturation, new_segment_gaps, make_flags


def test_saturation_long_run_found_short_run_ignored():
    fs = 10000.0 #min_ms=1 -> runs need 10 samples
    x = np.linspace(-0.5, 0.5, 1000, dtype=np.float32)
    x[100:150] = 2.0 #50 samples pinned at the max: a real saturation
    x[400:405] = 2.0 #5 samples at the same max: too short to count
    y = np.linspace(-1, 1, 1000, dtype=np.float32) #max and min only once, nothing flat
    stream = {"data":np.vstack([x, y]), "fs":fs, "start_time":0.5, "channels":np.array([3, 7])}
    flags = find_saturation(stream, min_ms=1.0)
    np.testing.assert_array_equal(flags["channel"], [3]) #channel number, not row index
    np.testing.assert_allclose(flags["start"], [0.5 + 100 / fs])
    np.testing.assert_allclose(flags["end"], [0.5 + 150 / fs]) #end is just after the last flat sample


def test_saturation_at_min_and_max():
    x = np.zeros(200, dtype=np.float32)
    x[:] = np.sin(np.arange(200)) * 0.1
    x[10:30] = 1.0
    x[100:130] = -1.0
    flags = find_saturation({"data":x[None, :], "fs":1000.0, "start_time":0.0, "channels":[1]}, min_ms=5.0)
    assert len(flags["channel"]) == 2
    np.testing.assert_allclose(sorted(flags["start"]), [0.010, 0.100])


def test_saturation_empty_stream():
    flags = find_saturation({"data":np.zeros((2, 0), dtype=np.float32), "fs":1000.0, "start_time":0.0, "channels":[1, 2]})
    assert flags["channel"].size == 0 and flags["start"].dtype == np.float64


def write_bv(folder, marker_lines):
    #a tiny brainvision header + marker pair; only the lines new_segment_gaps reads are filled in
    vhdr = folder / "rec.vhdr"
    vhdr.write_text("Brain Vision Data Exchange Header File Version 1.0\n[Common Infos]\nDataFile=rec.eeg\nMarkerFile=rec.vmrk\n"
                    "SamplingInterval=100\n") #100 us = 10 kHz
    (folder / "rec.vmrk").write_text("Brain Vision Data Exchange Marker File, Version 1.0\n[Marker Infos]\n" + "\n".join(marker_lines) + "\n")
    return vhdr


def test_new_segment_gaps(tmp_path):
    vhdr = write_bv(tmp_path, [
        "Mk1=New Segment,,1,1,0,20240101120000000000",
        "Mk2=Stimulus,S 13,5000,1,0", #ordinary marker, ignored
        "Mk3=New Segment,,10001,1,0,20240101120001000000", #1 s of samples and 1 s of clock: not a gap
        "Mk4=New Segment,,20001,1,0,20240101120002500000", #1 s of samples but 1.5 s of clock: 0.5 s missing
    ])
    gaps = new_segment_gaps(vhdr)
    np.testing.assert_allclose(gaps["time"], [1.0, 2.0])
    np.testing.assert_allclose(gaps["missing"], [0.0, 0.5], atol=1e-9)
    np.testing.assert_array_equal(gaps["is_gap"], [0, 1])


def test_new_segment_gaps_no_markers(tmp_path):
    vhdr = write_bv(tmp_path, ["Mk1=New Segment,,1,1,0,20240101120000000000"]) #only one segment: nothing to compare
    assert new_segment_gaps(vhdr)["time"].size == 0
    (tmp_path / "rec.vmrk").unlink()
    assert new_segment_gaps(vhdr)["is_gap"].size == 0 #missing .vmrk gives empty, not a crash


def test_make_flags_dnu(tmp_path):
    stream = {"data":np.zeros((1, 10), dtype=np.float32), "fs":1000.0, "start_time":0.0, "channels":[1]}
    flags = make_flags({"streams":{"ecos":stream}}, Path(tmp_path) / "block3_DNU")
    assert flags["dnu"] == 1 and "saturation_ecos" in flags and "new_segments" not in flags
    assert make_flags({"streams":{}}, Path(tmp_path) / "block3")["dnu"] == 0

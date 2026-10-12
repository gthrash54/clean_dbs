import tdt
import numpy as np
from pathlib import Path
import io
import os
import contextlib
from datetime import datetime
from zoneinfo import ZoneInfo

#about how many MB of samples (float32, all streams together) clean_block holds at once; keeps a job's memory flat whatever the block size
WINDOW_MB = 128

#the rigs are in Birmingham; block times are stored as seconds since 1970 (UTC), so they're turned into local clock time here,
#not in the time zone of whatever computer runs this (tdt's own start_date does that, and Cheaha runs on UTC)
RIG_TZ = ZoneInfo("America/Chicago")


def read_tdt(block_dir, t1=0, t2=0):
    """
    This function is used to read tdt files into workable data and is compatable with both python and Matlab
    block_dir = path to the folder holding .tev/.tsq files
    t1 = start time in s (0=beginning)
    t2 = end time in seconds (0=read to the end), defaults to 0
    Returns: streams, epocs, info
        streams:data(channels x samples, in volts), fs, start time, channels
        epocs: onset, offset, data. camera, stim info etc
        info: subject, block, experiment, duration, start
    Contact Garrett Thrash for any issues garrettthrash54@gmail.com
    """
    block_dir = Path(block_dir)
    _check_block(block_dir, t1, t2)
    with contextlib.redirect_stdout(io.StringIO()): #tdt prints the full path, which includes the subject folder
        headers = tdt.read_block(block_dir, headers=1)
        _check_tev(block_dir, headers)
        raw = tdt.read_block(block_dir, t1=t1, t2=t2, headers=headers)
    out = _events(raw, headers)
    settings = _store_settings(block_dir)
    #ascontiguousarray: tdt hands back multi-channel arrays in a memory order MATLAB refuses to convert
    for name, s in raw.streams.items():
        data, units = _to_volts(np.atleast_2d(s.data), settings.get(name, {}), name, block_dir)
        out["streams"][name] = {"data":np.ascontiguousarray(data), "fs":s.fs, "start_time":s.start_time, "channels":np.asarray(s.channel, dtype = np.int64),
                                "device":settings.get(name, {}).get("device", "unknown"), "units":units, "rail":_rail(s.data.dtype, settings.get(name, {}))}
    if np.isnan(out["info"]["duration"]) and t1 == 0 and t2 == 0: #.sev blocks: the whole streams were just read, so use them
        out["info"]["duration"] = max([s["data"].shape[1] / s["fs"] for s in out["streams"].values()], default=np.nan)
    return out


def tdt_windows(block_dir, window_mb=WINDOW_MB):
    """
    Reads a whole tdt block one time window at a time, so a block bigger than memory can still be cleaned (clean_block uses it).
    Input:
    block_dir: the block folder
    window_mb: about how many MB of samples (float32, all streams together) to hold at once
    Returns: (result, windows)
      result: read_tdt's dictionary for the whole block, except each stream has n_samples instead of data
      windows: yields {stream name: (first sample, channels x samples in volts)} for back-to-back windows, in time order
    Why not tdt.read_block(t1=, t2=) per window: it turns t1 into a sample from the packet timestamps, which carry ~1e-7 s of
    rounding, so some windows came back one sample early (checked on uh3021 clinical_dbs1). Instead the samples are read straight
    from the .tev, packet by packet, using the .tsq index tdt already parsed; a stream is each channel's packets end to end,
    which is exactly what a whole-block read_tdt gives.
    """
    block_dir = Path(block_dir)
    _check_block(block_dir, 0, 0)
    with contextlib.redirect_stdout(io.StringIO()):
        headers = tdt.read_block(block_dir, headers=1) #the .tsq index: per store, every packet's time, channel and place in the .tev
        _check_tev(block_dir, headers)
        out = _events(tdt.read_block(block_dir, headers=headers, evtype=["epocs", "scalars"]), headers)
        first = tdt.read_block(block_dir, headers=headers, evtype=["streams"], t2=1) #one second, for each stream's start_time as tdt gives it
    stores = {name:s for name, s in headers.stores.items() if s.type_str == "streams"}
    if any(s.ucf for s in stores.values()): #.sev stores keep their samples outside the .tev, so read those blocks whole
        out = read_tdt(block_dir)
        whole = {name:(0, s.pop("data")) for name, s in out["streams"].items()}
        for name, (_, x) in whole.items():
            out["streams"][name]["n_samples"] = x.shape[1]
        return out, iter([whole])
    settings = _store_settings(block_dir)
    packets, bytes_per_s = {}, 0.0
    for name, s in stores.items():
        packets[name] = _packets(s)
        n_ch, n_packets, per_packet = packets[name][0].shape[0], packets[name][0].shape[1], packets[name][2]
        bytes_per_s += n_ch * s.fs * 4
        start_time = first.streams[name].start_time if name in first.streams.keys() else 0.0
        out["streams"][name] = {"n_samples":n_packets * per_packet, "fs":s.fs, "start_time":start_time, "channels":np.asarray(packets[name][1], dtype=np.int64),
                                "device":settings.get(name, {}).get("device", "unknown"), "units":"V", "rail":_rail(packets[name][3], settings.get(name, {}))}
    del first
    #same length of time per window for every stream, rounded to whole 65536-sample hdf5 chunks (see save_h5)
    seconds = window_mb * 1e6 / max(bytes_per_s, 1.0)
    step = {name:max(1, int(seconds * s.fs) // 65536) * 65536 for name, s in stores.items()}
    n_windows = max([-(-out["streams"][name]["n_samples"] // step[name]) for name in stores] + [1])
    def windows():
        with open(headers.tev_path, "rb") as tev:
            for k in range(n_windows):
                window = {}
                for name in stores:
                    a, b = k * step[name], min((k + 1) * step[name], out["streams"][name]["n_samples"])
                    if b > a:
                        x = _read_packets(tev, *packets[name], a, b)
                        window[name] = (a, np.ascontiguousarray(_to_volts(x, settings.get(name, {}), name, block_dir)[0]))
                yield window
    return out, windows()


def _packets(store):
    """
    Input: one stream store from tdt.read_block(headers=1)
    Returns: (offsets, channels, per_packet, dtype): offsets is channels x packets, where each packet sits in the .tev (bytes),
             in the order a whole-block read joins them
    Keeps the packets a whole-block read keeps: time >= 0, plus one packet per channel just before that if there is one.
    """
    dtype = np.dtype(tdt.ALLOWED_FORMATS[store.dform])
    per_packet = (int(store.size) - 10) * 4 // dtype.itemsize #a packet is 10 words of header then samples, size counts 4-byte words
    ts, offsets = np.asarray(store.ts), np.asarray(store.data, dtype=np.uint64)
    chan = np.asarray(store.chan if len(store.chan) == len(ts) else np.repeat(store.chan, len(ts)), dtype=np.int64)
    keep = np.flatnonzero(ts >= 0)
    if keep.size and keep[0] - int(chan.max()) > -1:
        keep = np.concatenate((np.arange(keep[0] - int(chan.max()), keep[0]), keep))
    chan, offsets = chan[keep], offsets[keep]
    channels = sorted(set(chan.tolist()))
    counts = [int((chan == c).sum()) for c in channels]
    if len(set(counts)) > 1: #a whole-block read can't line these up either
        raise ValueError(f"{store.name}: channels have different numbers of packets {dict(zip(channels, counts))}")
    return np.vstack([offsets[chan == c] for c in channels]), channels, per_packet, dtype


def _read_packets(tev, offsets, channels, per_packet, dtype, a, b):
    #samples a to b (whole-block numbering) of one stream, channels x samples, raw values as stored
    first, last = a // per_packet, -(-b // per_packet)
    out = np.zeros((len(channels), (last - first) * per_packet), dtype=dtype)
    size = per_packet * dtype.itemsize
    for row in range(len(channels)):
        for j, offset in enumerate(offsets[row, first:last]):
            tev.seek(int(offset))
            if tev.readinto(out[row, j * per_packet:(j + 1) * per_packet]) < size: #packet missing from the .tev: zeros, like tdt
                out[row, j * per_packet:(j + 1) * per_packet] = 0
    return out[:, a - first * per_packet:b - first * per_packet]


def tev_problem(block_dir, headers=None):
    """
    Checks the .tev holds every stream packet the .tsq index points to. A copy cut short (Box download stopped partway) still
    opens, and tdt fills the missing packets with zeros without a word, so this has to be checked before reading.
    Input:
    block_dir: a tdt block folder (with its .tsq and .tev)
    headers: tdt.read_block(block_dir, headers=1), read here if left out
    Returns: "" if the .tev is complete (or the streams are in .sev files), otherwise what's wrong
    The last byte needed is the largest packet offset + (size - 10) * 4 over the stream stores (size counts 4-byte words,
    10 of them header); checked against the size of complete and truncated copies of uh3021 clinical_dbs1.
    """
    block_dir = Path(block_dir)
    if headers is None:
        with contextlib.redirect_stdout(io.StringIO()):
            headers = tdt.read_block(block_dir, headers=1)
    needed = 0
    for _, s in headers.stores.items(): #tdt's StructType.values() comes back empty, items() works
        if s.type_str == "streams" and not s.ucf and len(s.data):
            needed = max(needed, int(np.max(np.asarray(s.data, dtype=np.uint64))) + (int(s.size) - 10) * 4)
    if needed == 0:
        return ""
    tev = Path(headers.tev_path)
    size = os.path.getsize(tev) if tev.is_file() else 0
    if size < needed:
        return f"{tev.name} is cut short: {size} bytes but the .tsq index needs {needed} ({needed - size} bytes missing), copy the block again from Box"
    return ""


def _check_tev(block_dir, headers):
    problem = tev_problem(block_dir, headers)
    if problem:
        raise ValueError(problem)


def _check_block(block_dir, t1, t2):
    if not block_dir.is_dir():
        raise FileNotFoundError(f"Folder Not Found: {block_dir}, check path and spelling")
    #a block folder holds its own .tsq file; a surgery folder (e.g. stage1 surgery) only holds block folders
    if not any(block_dir.glob("*.tsq")):
        raise FileNotFoundError(f"No .tsq file in {block_dir}, point to a single block folder (e.g. .../stage1 surgery/clinical_dbs1)")
    #the .tev holds the signal; partial copies of the data sometimes only include the small index files
    if not any(block_dir.glob("*.tev")) and not any(block_dir.glob("*.sev")):
        raise FileNotFoundError(f"No .tev/.sev signal file in {block_dir}, only the index files were copied; get the full block from Box")
    if t1 < 0 or (t2 != 0 and t2 <= t1):
        raise ValueError(f"Bad time window t1={t1}, t2={t2}: need t1 >= 0 and t2 = 0 (to the end) or t2 > t1")


def _events(raw, headers):
    #everything except the streams, from a tdt.read_block result and its .tsq headers
    out = {"streams":{},"epocs":{},"scalars":{},"info":{} }
    #raw epocs are the event markers ie Cam 1
    warnings = []
    for name, e in raw.epocs.items():
        onset, offset, data = _match_epoc(name, e, headers, warnings)
        out["epocs"][name] = {"onset":onset, "offset":offset, "data":data}# data is value tied to each event ie frame number for cam 1
    #scalars: one vector of values per event, e.g. eS1p = the stimulator's 24 settings for every stim pulse (amplitude, contact, ...)
    for name, sc in raw.scalars.items():
        out["scalars"][name] = {"ts":np.asarray(sc.ts, dtype=np.float64), "data":np.atleast_2d(np.asarray(sc.data, dtype=np.float64))}
    start, stop = float(np.ravel(headers.start_time)[0]), float(np.ravel(headers.stop_time)[0])
    #the header's stop time is missing on blocks that didn't end cleanly and nonsense on some (uh3003 impedence_test gave about -1.7e9 s)
    header_s = stop - start if np.isfinite(stop) and 0 < stop - start < 7 * 86400 else np.nan
    #the streams say how long the recording really is, so duration comes from them; the header's is kept as duration_header
    streams_s = _stream_seconds(headers)
    out["info"] = {"subject":raw.info.subject, "block": raw.info.blockname, "experiment":raw.info.experiment,
                   "start":datetime.fromtimestamp(start, RIG_TZ).replace(tzinfo=None).isoformat(), #rig clock time, like StoresListing.txt
                   "duration":streams_s, "duration_header":header_s}
    if warnings:
        out["info"]["epoc_warning"] = "; ".join(warnings)
    return out


def _match_epoc(name, e, headers, warnings):
    """
    Input: epoc name, tdt's epoc, the .tsq headers, list to add a warning to
    Returns: onset, offset, data with one entry per event
    tdt sometimes gives a settings store one onset more than it has values: on uh3001 and uh3006 CntA (14 local blocks, 2026-10-12)
    it adds the 0 s onset of the PeA_ store it rides on (2602 onsets, 2601 values). Onsets the .tsq doesn't hold for this store are
    dropped; if the counts still don't match, the extra events at the end are cut. Either way it's written in info epoc_warning.
    """
    onset, offset, data = np.asarray(e.onset), np.asarray(e.offset), np.asarray(e.data)
    if len(onset) == len(data) == len(offset):
        return onset, offset, data
    n_on = len(onset)
    if len(onset) > len(data) and name in headers.stores.keys():
        mine = np.isin(onset, np.ravel(headers.stores[name].onset))
        if mine.sum() == len(data):
            onset = onset[mine]
    n = min(len(onset), len(offset), len(data))
    onset, offset, data = onset[:n], offset[:n], data[:n]
    warnings.append(f"{name}: tdt gave {n_on} onsets for {len(e.data)} values, kept {n} events")
    return onset, offset, data


def _stream_seconds(headers):
    #length of the longest stream (n_samples / fs) from the .tsq index, the same samples a whole-block read gives; NaN if the samples are in .sev files
    seconds = []
    for _, s in headers.stores.items(): #tdt's StructType.values() comes back empty, items() works
        if s.type_str == "streams" and not s.ucf and len(s.data):
            offsets, _, per_packet, _ = _packets(s)
            seconds.append(offsets.shape[1] * per_packet / s.fs)
    return max(seconds, default=np.nan)

#Synapse's "Scale" setting for Integer-16 stores: the whole numbers are in these units (checked on uh3002, where ecog sits around 33 = microvolts)
SCALES = {"Unity":1.0, "Milli":1e-3, "Micro":1e-6, "Nano":1e-9}


def _store_settings(block_dir):
    """
    Reads StoresListing.txt (written by Synapse next to every block) for each store's device, format and scale.
    tdt.read_block doesn't read these, but older blocks store streams as whole numbers in uV or mV, so we need them.
    Returns: {store name: {"device":..., "format":..., "scale":...}}, empty if there's no StoresListing.txt
    """
    listing = block_dir / "StoresListing.txt"
    settings = {}
    if not listing.is_file():
        return settings
    device, store = "", None
    for line in listing.read_text(errors="ignore").splitlines():
        key, _, value = line.partition(":")
        key, value = key.strip(), value.strip()
        if key == "Object ID":
            device = value.split(" - ")[0] #"ecog - Neural Stream Processor" -> "ecog"
        elif key == "Store ID":
            store = tdt.fix_var_name(value) if hasattr(tdt, "fix_var_name") else value.replace("/", "_") #match tdt's names ("PeA/" -> "PeA_")
            settings[store] = {"device":device}
        elif store is not None and key in ("Format", "Scale"):
            settings[store][key.lower()] = value
    return settings


#the amplifier's ceiling: Float-32 stores pin at 0.53687084 V (seen on uh3021 clinical_dbs1 ecos), Integer-16 stores at +-32000 counts
FLOAT_RAIL, INT_RAIL = 0.53687084, 32000


def _rail(dtype, setting):
    """
    Input: dtype the store is saved as, its StoresListing.txt setting
    Returns: the rail in volts (flags.Saturation flags a channel sitting there the whole block), NaN if unknown
    """
    if np.dtype(dtype).kind == "f":
        return FLOAT_RAIL
    return INT_RAIL * SCALES[setting["scale"]] if setting.get("scale") in SCALES else np.nan


def _to_volts(data, setting, name, block_dir):
    """
    Float-32 stores are already volts. Integer-16 stores are whole numbers in the unit Synapse's Scale says, so multiply into volts.
    Returns: (float32 data in volts, units label)
    """
    scale = setting.get("scale")
    if data.dtype.kind == "f" and scale in (None, "Unity"): #every Float-32 store in the inventory so far is Unity
        return data.astype(np.float32, copy=False), "V"
    if scale not in SCALES:
        raise ValueError(f"{name} in {block_dir.name} is stored as {data.dtype} with scale {scale!r}, can't convert to volts without StoresListing.txt")
    return (data.astype(np.float32) * np.float32(SCALES[scale])), "V"

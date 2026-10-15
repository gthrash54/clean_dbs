import re
import tempfile
from datetime import datetime
import numpy as np
import mne
from pathlib import Path
from clean_dbs.read_tdt import WINDOW_MB
from clean_dbs.bv_channels import channel_types, cohort_of, read_vhdr_channels, VOLT_UNITS

#units MNE scales to SI on its own (it divides by 1e6 for uS etc); any unit not here (BIP, %, ...) comes through as the raw number
MNE_SI_UNITS = {"µS":"S", "uS":"S", "S":"S", "C":"C", "°C":"C", "N":"N", "ARU":"ARU", "n/a":"n/a"}
#per-channel labels read_bv puts on each stream, all lists with one entry per channel
CHANNEL_KEYS = ("group", "contact", "side", "muscle", "subtype", "units", "source")
#actiCHamp input range: 24 bit ADC at 0.048828125 uV per step = +-409.6 mV per electrode. Seen on ug3003 stage2 erp001 (2026-10-12): DBS channels
#bottom out at -814 to -819 mV, about twice that, because the recorder subtracts the reference electrode (Ch17) in software, so a referenced
#channel can reach +-819 mV when both electrodes rail. A channel flat for the whole block at or past 409.6 mV is railed.
ACTICHAMP_RAIL = 0.4096


def read_bv(vhdr_path, t1=0, t2=0, split=True):
    """
    This function read_bv, can be used in python or matlab to read brainvision files
      (.vhdr, .vmrk) in simple arrays and metadata. Built for Harrison Walker's lab
      Inputs:
      vhdr_Path: path to vhdr file(.eeg,.vmrk must sit next to it)
      t1 = start time in s (0=beginning)
      t2 = end time in seconds (0=read to the end), defaults to 0
      split = True (default) gives one stream per channel type, False gives everything in one stream "dbs" like before
      Returns:
      streams: one per channel type from the channel names (see bv_channels.py): lead, ecog, eeg, emg, audio, mer, accel, joystick,
               force, aux, unused, unknown. Each has data (channels x samples, volts), fs, start_time, channels (position in the .vhdr, 1-based),
               channel_names, type, subtype, units, rail, and per channel: channel_group, channel_contact, channel_side, channel_muscle,
               channel_subtype, channel_units, channel_source. Channels whose units aren't volts get their own stream (e.g. "aux_BIP") with
               the units as written, never relabeled V.
      epocs: one per marker code ("S15", "R11", "BtnPressed", "NewSegment"); what each code means is worked out in markers.py
      info: block, start (local clock time, like read_tdt), duration, cohort (lowercase study ID prefix like cohorts.py, "" if none)
      Contact Garrett Thrash for any issues garrettthrash54@gmail.com
    """
    vhdr_path = Path(vhdr_path)
    raw, header = _open(vhdr_path, t1, t2)
    t2 = None if t2 == 0 else t2 #bc tmax=0 means stop at 0 seconds in mne
    out, picks = _result(raw, vhdr_path, header, t1, t2, split)
    for name, p in picks.items():
        out["streams"][name]["data"] = raw.get_data(picks=p, tmin=t1, tmax=t2).astype(np.float32)
    return out


def bv_windows(vhdr_path, window_mb=WINDOW_MB, split=True):
    """
    Reads a whole brainvision recording one time window at a time, so a recording bigger than memory can still be cleaned (clean_block uses it).
    Input:
    vhdr_path: the .vhdr (.eeg, .vmrk next to it)
    window_mb: about how many MB of samples (float32) to hold at once
    split: like read_bv, True = one stream per channel type
    Returns: (result, windows)
      result: read_bv's dictionary for the whole recording, except each stream has n_samples instead of data
      windows: yields {stream name: (first sample, channels x samples in volts)} for back-to-back windows, in time order
    MNE reads by sample number here, so the windows join up exactly; each window is read for all channels at once, then split into the
    typed streams, and is float32 straight away: the same numbers read_bv gives.
    """
    vhdr_path = Path(vhdr_path)
    raw, header = _open(vhdr_path, 0, 0)
    out, picks = _result(raw, vhdr_path, header, 0, None, split)
    for name in picks:
        out["streams"][name]["n_samples"] = raw.n_times
    #a whole number of 65536-sample hdf5 chunks per window (see save_h5), so each chunk is written in one go
    step = max(1, int(window_mb * 1e6 / (4 * len(raw.ch_names))) // 65536) * 65536
    def windows():
        for start in range(0, raw.n_times, step):
            x = raw.get_data(start=start, stop=min(start + step, raw.n_times)).astype(np.float32)
            yield {name:(start, x[p]) for name, p in picks.items()}
            x = None #so the window is freed before the next one is read
    return out, windows()


def _open(vhdr_path, t1, t2):
    #checks the files and the window, then opens the recording without reading the samples. Returns: (mne raw, .vhdr text)
    if not vhdr_path.is_file():
        raise FileNotFoundError(f"File Not Found: {vhdr_path}, check path and spelling")
    header = vhdr_path.read_text(errors="ignore")
    #the .vhdr names its signal file on the "DataFile=" line; partial copies of the data sometimes leave it out
    data_file = next((line.split("=", 1)[1].strip() for line in header.splitlines() if line.startswith("DataFile=")), vhdr_path.stem + ".eeg")
    if not (vhdr_path.parent / data_file).is_file():
        raise FileNotFoundError(f"Signal file {data_file} is missing next to {vhdr_path}; only the header was copied, get the full recording from Box")
    if (vhdr_path.parent / data_file).stat().st_size == 0: #DLPFC085 burstpp_trial2: the .eeg on Box is 0 bytes, the recording never saved
        raise ValueError(f"Signal file {data_file} is empty (0 bytes, also on Box): the recording has no data, nothing to clean")
    if t1 < 0 or (t2 != 0 and t2 <= t1):
        raise ValueError(f"Bad time window t1={t1}, t2={t2}: need t1 >= 0 and t2 = 0 (to the end) or t2 > t1")
    raw = _read_raw(vhdr_path, header)
    if raw.n_times == 0: #a data file too short to hold one sample
        raise ValueError(f"Signal file {data_file} has no samples: the recording has no data, nothing to clean")
    if t1 >= raw.times[-1]:
        raise ValueError(f"t1={t1} is past the end of the recording ({raw.times[-1]:.1f} s)")
    return raw, header


def _read_raw(vhdr_path, header):
    """
    Opens the recording with MNE (samples stay on disk). A header that gives two channels the same name (BI002 12 Month Baseline has
    LFP1 on Ch20 and Ch25) makes MNE fail ("Channel LFP1-0 has no associated original unit"), so MNE gets a copy of the header with
    the repeats renamed, in a temporary folder, pointing at the real .eeg and .vmrk. Only MNE sees those names: the labels in the
    cleaned file come from the original header by channel position (bv_channels.channel_types).
    """
    names = [n.strip() for n, _ in read_vhdr_channels(vhdr_path)]
    if len(set(names)) == len(names):
        return mne.io.read_raw_brainvision(vhdr_path, preload=False, verbose="error") #verbose="error" hides MNE's progress messages
    text = vhdr_path.read_bytes().decode("utf-8", "surrogateescape") #keeps any byte as it was, whatever the codepage
    seen = {}
    def rename(m):
        name = m.group(2).split(",")[0]
        seen[name.strip()] = seen.get(name.strip(), 0) + 1
        return m.group(0) if seen[name.strip()] == 1 else f"{m.group(1)}{name.strip()}_dup{seen[name.strip()]}{m.group(2)[len(name):]}"
    text = re.sub(r"^(Ch\d+=)([^\r\n]*)", rename, text, flags=re.M)
    for key in ("DataFile", "MarkerFile"): #full paths, so the copy finds the files next to the real header
        text = re.sub(rf"^({key}=)([^\r\n]*)", lambda m: m.group(1) + str(vhdr_path.parent.resolve() / m.group(2).strip()), text, flags=re.M)
    with tempfile.TemporaryDirectory() as folder:
        copy = Path(folder) / vhdr_path.name
        copy.write_bytes(text.encode("utf-8", "surrogateescape"))
        return mne.io.read_raw_brainvision(copy, preload=False, verbose="error") #the header and markers are read here, the .eeg later


def _result(raw, vhdr_path, header, t1, t2, split):
    """
    Everything read_bv returns except the samples.
    Returns: (result, picks): picks = {stream name: 0-based channel positions in the file}, so read_bv and bv_windows cut the same streams
    """
    out = {"streams":{}, "epocs":{}, "scalars":{}, "info":{}}
    fs = raw.info["sfreq"]
    rails = _rails(header, len(raw.ch_names))
    if not split:
        picks = {"dbs":list(range(len(raw.ch_names)))}
        out["streams"]["dbs"] = {"fs":fs, "start_time":t1, "channels":np.arange(1, len(raw.ch_names) + 1), "device":"brainvision", "units":"V",
                                 "rail":_one(rails)}
    else:
        chans = channel_types(vhdr_path)
        if len(chans) != len(raw.ch_names): #header and MNE disagree: don't label anything, keep one unknown stream
            chans = [{"number":i + 1, "name":n, "units":"", "stream":"unknown", "type":"unknown", **{k:"" for k in CHANNEL_KEYS if k != "units"}}
                     for i, n in enumerate(raw.ch_names)]
        picks = {}
        for stream in dict.fromkeys(c["stream"] for c in chans): #keeps header order
            rows = [c for c in chans if c["stream"] == stream]
            picks[stream] = [c["number"] - 1 for c in rows]
            units = [_units_out(c["units"]) for c in rows]
            subtypes = {c["subtype"] for c in rows}
            out["streams"][stream] = {"fs":fs, "start_time":t1, "channels":np.array([c["number"] for c in rows]), "device":"brainvision",
                                      "units":units[0] if len(set(units)) == 1 else "mixed", "type":rows[0]["type"],
                                      "subtype":subtypes.pop() if len(subtypes) == 1 else "", "channel_names":[c["name"] for c in rows],
                                      "rail":_one([rails[p] if u == "V" else np.nan for p, u in zip(picks[stream], units)])}
            for key in CHANNEL_KEYS:
                out["streams"][stream][f"channel_{key}"] = units if key == "units" else [c[key] for c in rows]
    full = {}
    for name in np.unique(raw.annotations.description): #epocs loop
        pick = (raw.annotations.description == name) & (raw.annotations.onset >= t1)
        if t2 is not None: #t2 is None when reading to the end, and comparing with None crashes
            pick = pick & (raw.annotations.onset < t2)
        key = name.split("/")[-1].replace(" ", "") #"Response/R 15" -> "R15", slashes/spaces aren't allowed in MATLAB field names
        if key == "": #markers with no description (e.g. "New Segment/") fall back to their type
            key = name.split("/")[0].replace(" ", "")
        if key in full and full[key] != name: #e.g. "Stimulus/S 15" and "S15" would both be S15: keep them apart instead of overwriting
            key = re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")
        full[key] = name
        out["epocs"][key] = {"onset":raw.annotations.onset[pick], "offset":raw.annotations.onset[pick] + raw.annotations.duration[pick], "marker":name}
    out["info"] = {"block":vhdr_path.stem, "start":_start(vhdr_path, header), "duration":raw.n_times / fs, "cohort":cohort_of(vhdr_path).lower()}
    return out, picks


def _units_out(units):
    #what the numbers are in after MNE: volts for any voltage unit, SI for the units MNE knows, otherwise the unit as written (not scaled)
    return "V" if units in VOLT_UNITS else MNE_SI_UNITS.get(units, units)


def _rails(header, n_channels):
    """
    Input: .vhdr text, channel count
    Returns: the amplifier ceiling in volts per channel (flags.Saturation flags a channel sitting there the whole block), NaN where unknown
    INT_16 files: the biggest number the file can hold (32767) x the channel's resolution. Float files from an actiCHamp: ACTICHAMP_RAIL.
    Any other amplifier writing floats: NaN, we don't know its range (every .vhdr in the local inventory is an actiCHamp, 2026-10-12).
    """
    binary = re.search(r"^BinaryFormat=(\S+)", header, re.M)
    binary = binary.group(1).upper() if binary else "INT_16" #brainvision's default
    rails = np.full(n_channels, np.nan)
    for m in re.finditer(r"^Ch(\d+)=([^\r\n]*)$", header, re.M):
        k, fields = int(m.group(1)) - 1, m.group(2).split(",")
        if k >= n_channels:
            continue
        units = fields[3].strip() if len(fields) > 3 and fields[3].strip() else "µV"
        if binary == "INT_16" and units in VOLT_UNITS:
            resolution = float(fields[2]) if len(fields) > 2 and fields[2].strip() else 1.0
            rails[k] = 32767 * resolution * VOLT_UNITS[units]
        elif binary == "IEEE_FLOAT_32" and "actichamp" in header.lower():
            rails[k] = ACTICHAMP_RAIL
    return rails


def _one(rails):
    #a stream's rail: the same for every channel, else NaN (unknown)
    rails = np.asarray(rails, dtype=np.float64)
    return float(rails[0]) if len(rails) and np.all(rails == rails[0]) else np.nan


def _start(vhdr_path, header):
    """
    Returns: when the recording started, local clock time, in read_tdt's format ("2025-01-31T09:15:02.123456")
    The .vmrk's first dated New Segment is the recorder computer's clock, which is local time already, so it's taken as written, no time
    zone shift (MNE labels it UTC; using it straight would be off by 5-6 h if anything converted it). Without one (no .vmrk, ug3012
    thum_propofol) the .vhdr impedance line's "recording started at hh:mm:ss" gives the time of day only ("12:17:05"), or blank.
    """
    marker_file = next((l.split("=", 1)[1].strip() for l in header.splitlines() if l.startswith("MarkerFile=")), vhdr_path.stem + ".vmrk")
    vmrk = vhdr_path.parent / marker_file
    if vmrk.is_file():
        m = re.search(r"^Mk\d+=New Segment,[^,]*,\d+,\d+,[^,]*,(\d{20})", vmrk.read_text(errors="ignore"), re.M)
        if m:
            return datetime.strptime(m.group(1), "%Y%m%d%H%M%S%f").isoformat()
    m = re.search(r"recording started at (\d\d:\d\d:\d\d)", header)
    return datetime.strptime(m.group(1), "%H:%M:%S").time().isoformat() if m else ""

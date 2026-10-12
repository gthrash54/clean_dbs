import numpy as np
import mne
from pathlib import Path

def read_bv(vhdr_path, t1=0, t2=0):
    """
    This function read_bv, can be used in python or matlab to read brainvision files
      (.vhdr, .vmrk) in simple arrays and metadata. Built for Harrison Walker's lab
      Inputs:
      vhdr_Path: path to vhdr file(.eeg,.vmrk must sit next to it)
      t1 = start time in s (0=beginning)
      t2 = end time in seconds (0=read to the end), defaults to 0
      Returns:
      streams: dbs, channels x samples in volts
      epocs: marker names depend on the rig (R15/R11 on uh3021 = stim pulse start / +60 µs; other rigs use S13/S15 etc.)
      info: block: path file, start time, duration
      Contact Garrett Thrash for any issues garrettthrash54@gmail.com
    """
    vhdr_path = Path(vhdr_path)
    if not vhdr_path.is_file():
        raise FileNotFoundError(f"File Not Found: {vhdr_path}, check path and spelling")
    #the .vhdr names its signal file on the "DataFile=" line; partial copies of the data sometimes leave it out
    data_file = next((line.split("=", 1)[1].strip() for line in vhdr_path.read_text(errors="ignore").splitlines() if line.startswith("DataFile=")), vhdr_path.stem + ".eeg")
    if not (vhdr_path.parent / data_file).is_file():
        raise FileNotFoundError(f"Signal file {data_file} is missing next to {vhdr_path}; only the header was copied, get the full recording from Box")
    if t1 < 0 or (t2 != 0 and t2 <= t1):
        raise ValueError(f"Bad time window t1={t1}, t2={t2}: need t1 >= 0 and t2 = 0 (to the end) or t2 > t1")
    raw = mne.io.read_raw_brainvision(vhdr_path, preload=False, verbose="error") #verbose="error" hides MNE's progress messages
    if t1 >= raw.times[-1]:
        raise ValueError(f"t1={t1} is past the end of the recording ({raw.times[-1]:.1f} s)")
    if t2 == 0:
        t2 = None                #bc tmax=0 means stop at 0 seconds in mne
    data = raw.get_data(tmin =t1, tmax=t2)
    out = {"streams":{},"epocs":{},"info":{} }
    out["streams"]["dbs"] = {"data":data.astype(np.float32), "fs":raw.info["sfreq"], "start_time":t1, "channels":np.arange(1, len(raw.ch_names)+1)}
    for name in np.unique(raw.annotations.description): #epocs loop
        pick = (raw.annotations.description == name) & (raw.annotations.onset >= t1)
        if t2 is not None: #t2 is None when reading to the end, and comparing with None crashes
            pick = pick & (raw.annotations.onset < t2)
        key = name.split("/")[-1].replace(" ", "") #"Response/R 15" -> "R15", slashes/spaces aren't allowed in MATLAB field names
        if key == "": #markers with no description (e.g. "New Segment/") fall back to their type
            key = name.split("/")[0].replace(" ", "")
        out["epocs"][key] ={"onset":raw.annotations.onset[pick], "offset":raw.annotations.onset[pick]+raw.annotations.duration[pick],}
    #MNE labels the start time as UTC, but the rig writes local clock time (matches the .vhdr "recording started at"), so drop the label to match read_tdt
    out["info"] = {"block":vhdr_path.stem, "start":raw.info["meas_date"].replace(tzinfo=None).isoformat(), "duration": raw.n_times/raw.info["sfreq"]}
    return out
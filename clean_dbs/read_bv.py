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
      epocs: named after markers (ie R15) R15 is start of a pulse and R11 is 60us later
      info: block: path file, start time, duration
      Contact Garrett Thrash for any issues garrettthrash54@gmail.com
    """
    vhdr_path = Path(vhdr_path)
    if not vhdr_path.is_file():
        raise FileNotFoundError(f"File Not Found: {vhdr_path}, check path and spelling")
    raw = mne.io.read_raw_brainvision(vhdr_path, preload=False, verbose="error") #verbose="error" hides MNE's progress messages
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
        out["epocs"][key] ={"onset":raw.annotations.onset[pick], "offset":raw.annotations.onset[pick]+raw.annotations.duration[pick],}
    #MNE labels the start time as UTC, but the rig writes local clock time (matches the .vhdr "recording started at"), so drop the label to match read_tdt
    out["info"] = {"block":vhdr_path.stem, "start":raw.info["meas_date"].replace(tzinfo=None).isoformat(), "duration": raw.n_times/raw.info["sfreq"]}
    return out
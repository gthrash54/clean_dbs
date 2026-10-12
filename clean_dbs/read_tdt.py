import tdt
import numpy as np
from pathlib import Path
import io
import contextlib
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
    out = {"streams":{},"epocs":{},"info":{} }
    with contextlib.redirect_stdout(io.StringIO()): #tdt prints the full path, which includes the subject folder
        raw = tdt.read_block(block_dir,t1=t1,t2=t2)
    #ascontiguousarray: tdt hands back multi-channel arrays in a memory order MATLAB refuses to convert
    for name, s in raw.streams.items():
        out["streams"][name] = {"data":np.ascontiguousarray(np.atleast_2d(s.data)), "fs":s.fs, "start_time":s.start_time, "channels":np.asarray(s.channel, dtype = np.int64)}
    #raw epocs are the event markers ie Cam 1
    for name, e in raw.epocs.items():
        out["epocs"][name] = {"onset":e.onset, "offset":e.offset, "data":e.data}# data is value tied to each event ie frame number for cam 1
    out["info"] = {"subject":raw.info.subject, "block": raw.info.blockname, "experiment":raw.info.experiment, "start": raw.info.start_date.isoformat(),
                          "duration":raw.info.duration.total_seconds()}
    return out
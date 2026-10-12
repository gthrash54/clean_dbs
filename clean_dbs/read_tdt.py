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
    settings = _store_settings(block_dir)
    #ascontiguousarray: tdt hands back multi-channel arrays in a memory order MATLAB refuses to convert
    for name, s in raw.streams.items():
        data, units = _to_volts(np.atleast_2d(s.data), settings.get(name, {}), name, block_dir)
        out["streams"][name] = {"data":np.ascontiguousarray(data), "fs":s.fs, "start_time":s.start_time, "channels":np.asarray(s.channel, dtype = np.int64),
                                "device":settings.get(name, {}).get("device", "unknown"), "units":units}
    #raw epocs are the event markers ie Cam 1
    for name, e in raw.epocs.items():
        out["epocs"][name] = {"onset":e.onset, "offset":e.offset, "data":e.data}# data is value tied to each event ie frame number for cam 1
    out["info"] = {"subject":raw.info.subject, "block": raw.info.blockname, "experiment":raw.info.experiment, "start": raw.info.start_date.isoformat(),
                          "duration":raw.info.duration.total_seconds()}
    return out

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

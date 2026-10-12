import os
import subprocess
from datetime import datetime
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path

import h5py
import numpy as np


def save_h5(result, out_path, source_dir):
    """
    Saves a reader result (from read_tdt or read_bv) to one HDF5 file that opens in Python and MATLAB.
    result = the dictionary from read_tdt / read_bv (streams, epocs, info)
    out_path = where the .h5 goes (folders are made if needed)
    source_dir = the raw block folder or .vhdr it came from, stored inside the file so you can trace it back
    Returns: out_path

    Layout inside the file:
      /streams/<name>/data       channels x samples, float32, volts
                                 (MATLAB's h5read flips it: you get samples x channels)
      /streams/<name>/channels   channel numbers
      /streams/<name>/channel_names  contact names when known (brainvision now, tdt once we have the wiring)
      /streams/<name>            attributes: fs, start_time, device, units, type, subtype
      /epocs/<name>/onset, offset, data   attribute: type (video, stim, marker, other)
      /stim/time, source, voice, channel, amplitude_ua, period_ms, count, confirmed (+ time_clock for brainvision)
      /scalars/<name>/ts, data          data is parameters x events (e.g. eS1p stim settings per pulse)
      /info                      attributes: whatever the reader put in info
      /flags                     saturation_<stream>/channel,start,end; new_segments/time,missing,is_gap; impedance_<stream>/channel,time,frequency,impedance,connected; attribute dnu (notes only, data unchanged)
      /                          attributes: provenance (source, clean_dbs version, package versions, when)

    It writes to <name>.h5.part first and only renames to .h5 once everything is written,
    so a crashed or killed job never leaves a half-written .h5 that looks finished.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    part_path = out_path.with_name(out_path.name + ".part")

    with h5py.File(part_path, "w") as f:
        for name, s in result["streams"].items():
            g = f.create_group(f"streams/{name}")
            _write_array(g, "data", np.asarray(s["data"], dtype=np.float32), signal=True)
            _write_array(g, "channels", np.asarray(s["channels"], dtype=np.int64))
            if "channel_names" in s: #h5py needs its own string type, plain numpy text arrays dont save
                g.create_dataset("channel_names", data=np.asarray(s["channel_names"], dtype=h5py.string_dtype()))
            g.attrs["fs"] = float(s["fs"])
            g.attrs["start_time"] = float(s["start_time"])
            g.attrs["device"] = str(s.get("device", "unknown")) #which Synapse gizmo / rig recorded it, e.g. ecogdbs, emg, microphone
            g.attrs["units"] = str(s.get("units", "V"))
            g.attrs["type"] = str(s.get("type", "unknown")) #lfp, mer, emg, audio, stim_monitor (see datatypes.py)
            g.attrs["subtype"] = str(s.get("subtype", ""))
        f.create_group("epocs")
        for name, e in result["epocs"].items():
            g = f.create_group(f"epocs/{name}")
            for key, values in e.items():
                if isinstance(values, str): #labels like type go on as attributes
                    g.attrs[key] = values
                else:
                    _write_array(g, key, np.asarray(values, dtype=np.float64))
        f.create_group("scalars")
        for name, sc in result.get("scalars", {}).items():
            g = f.create_group(f"scalars/{name}")
            _write_array(g, "ts", np.asarray(sc["ts"], dtype=np.float64))
            _write_array(g, "data", np.asarray(sc["data"], dtype=np.float64)) #rows = parameters, columns = events
            g.attrs["type"] = str(sc.get("type", "unknown"))
        stim = f.create_group("stim") #one row per stim pulse, see stim_events.py
        for key, values in result.get("stim", {}).items():
            values = np.asarray(values)
            if values.dtype.kind == "U": #text columns (source, voice) need h5py's string type
                stim.create_dataset(key, data=values.astype(object), dtype=h5py.string_dtype())
            else:
                _write_array(stim, key, values)
        info = f.create_group("info")
        for key, value in result["info"].items():
            info.attrs[key] = _attr_value(value)
        flags = f.create_group("flags") #notes about the data (saturation, new segments, dnu), the data itself isn't changed
        for name, value in result.get("flags", {}).items():
            if isinstance(value, dict):
                g = flags.create_group(name)
                for key, values in value.items():
                    _write_array(g, key, np.asarray(values))
            else:
                flags.attrs[name] = _attr_value(value)
        for key, value in _provenance(source_dir).items():
            f.attrs[key] = value

    os.replace(part_path, out_path)  #the rename is all-or-nothing on Windows, Mac and Linux
    return out_path


def _write_array(group, name, array, signal=False):
    #compress the big signal arrays; chunks span all channels so reading a time window stays fast
    if signal and array.size > 0:
        chunk = (array.shape[0], min(array.shape[1], 65536))
        group.create_dataset(name, data=array, chunks=chunk, compression="gzip", compression_opts=1, shuffle=True)
    else:
        group.create_dataset(name, data=array)


def _attr_value(value):
    #MATLAB reads Python True/False back as the text 'TRUE', so flags are stored as 0/1 numbers instead
    if isinstance(value, (bool, np.bool_)):
        return np.int8(value)
    if isinstance(value, (int, float, np.integer, np.floating)):
        return value
    return str(value)


def _provenance(source_dir):
    source_dir = Path(source_dir)
    repo = Path(__file__).resolve().parent.parent
    try:
        commit = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(repo), "status", "--porcelain"], capture_output=True, text=True, check=True).stdout.strip() != ""
    except (OSError, subprocess.CalledProcessError):
        commit, dirty = "unknown", True
    versions = {}
    for package in ("numpy", "h5py", "tdt", "mne"):
        try:
            versions[f"version_{package}"] = version(package)
        except PackageNotFoundError:
            versions[f"version_{package}"] = "not installed"
    #keep the source path from the study ID down (e.g. uh3021/stage1 surgery/clinical_dbs1), never the full path on someone's computer
    source_files = [source_dir] if source_dir.is_file() else list(source_dir.iterdir()) if source_dir.is_dir() else []
    return {
        "source": "/".join(source_dir.parts[-3:]),
        "source_bytes": int(sum(p.stat().st_size for p in source_files if p.is_file())),
        "clean_dbs_commit": commit,
        "clean_dbs_uncommitted_changes": np.int8(dirty),  #1 means the code that made this file wasn't committed yet
        "created": datetime.now().isoformat(timespec="seconds"),
        **versions,
    }

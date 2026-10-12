import os
import subprocess
from contextlib import contextmanager
from datetime import datetime
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path

import h5py
import numpy as np
from clean_dbs.block_info import source_path
from clean_dbs.version import CLEAN_VERSION


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
      /streams/<name>/channel_names, channel_group, channel_target, channel_contact, channel_source (+ brainvision: channel_subtype,
                                 channel_side, channel_muscle, channel_units): what each channel is (e.g. lead1_c3 / lead1 / STN / 3), and where that came from
      /streams/<name>            attributes: fs, start_time, device, units, type, subtype, rail (amplifier ceiling in volts, NaN = unknown)
                                 brainvision streams are one per channel type: lead, ecog, eeg, emg, audio, mer, accel, joystick, force, aux, unused, unknown
      /epocs/<name>/onset, offset, data   attributes: type (video, stim, task, marker, other); brainvision also meaning, evidence, kind (see markers.py)
      /stim/time, source, voice, channel, amplitude_ua, period_ms, count, confirmed (+ brainvision: time_clock, marker_lag_ms, and sidecar
                                 columns like pulse_width_us, amplitude_mv, ttl_ms, condition, protocol_* when a stimulator file lines up)
      /sidecars/<name>/<column>  stimulator / task files found next to a brainvision recording, as read; attributes kind, file, matched, match_*
      /scalars/<name>/ts, data          data is parameters x events (e.g. eS1p stim settings per pulse)
      /info                      attributes: whatever the reader put in info
      /flags                     saturation_<stream>/channel,start,end; new_segments/time,missing,is_gap; impedance_<stream>/channel,time,frequency,impedance,connected; attribute dnu (notes only, data unchanged)
      /                          attributes: clean_version (see version.py), provenance (source, source_bytes, clean_dbs commit, package versions, when)

    It writes to <name>.h5.part first and only renames to .h5 once everything is written,
    so a crashed or killed job never leaves a half-written .h5 that looks finished.
    Recordings too big for memory go through open_h5 + StreamWriter + write_rest instead (see clean_block); the file comes out the same.
    """
    with open_h5(out_path) as f:
        writer = StreamWriter(f, result["streams"])
        writer.write({name:(0, s["data"]) for name, s in result["streams"].items()})
        writer.close()
        write_rest(f, result, source_dir)
    return out_path


@contextmanager
def open_h5(out_path):
    """
    Input: out_path, the .h5 to make (folders are made if needed)
    Returns: (in a with-block) the open h5py file, which is really <name>.h5.part; it's renamed to .h5 only if the
             with-block finishes without an error, so a crashed or killed job never leaves a half-written .h5 that looks finished
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    part_path = out_path.with_name(out_path.name + ".part")
    with h5py.File(part_path, "w") as f:
        yield f
    os.replace(part_path, out_path)  #the rename is all-or-nothing on Windows, Mac and Linux


#samples per hdf5 chunk; chunks span all channels so reading a time window stays fast
CHUNK = 65536


class StreamWriter:
    """
    Writes the streams into an open .h5 one window at a time, so a recording never has to be in memory all at once.
    Input:
    f: the open file (from open_h5)
    streams: result["streams"]; each needs data, or n_samples when the samples come in windows (tdt_windows / bv_windows)
    Use: write({name: (first sample, channels x samples)}) for back-to-back windows, then close()
    The datasets are made at full size up front. Samples are held back until a whole 65536-sample chunk is ready, so every
    chunk is compressed and written once (writing half a chunk and coming back to it costs time and leaves dead space in the file).
    """
    def __init__(self, f, streams):
        self.datasets, self.held, self.done = {}, {}, {}
        for name, s in streams.items():
            n_ch, n = np.shape(s["data"]) if "data" in s else (len(s["channels"]), s["n_samples"])
            g = f.create_group(f"streams/{name}")
            if n_ch * n > 0: #compress the big signal arrays
                self.datasets[name] = g.create_dataset("data", shape=(n_ch, n), dtype=np.float32, chunks=(n_ch, min(n, CHUNK)),
                                                       compression="gzip", compression_opts=1, shuffle=True)
                self.held[name], self.done[name] = np.zeros((n_ch, 0), dtype=np.float32), 0
            else:
                g.create_dataset("data", shape=(n_ch, n), dtype=np.float32)
            _write_array(g, "channels", np.asarray(s["channels"], dtype=np.int64))
            for key in s:
                if key.startswith("channel_"): #h5py needs its own string type, plain numpy text arrays dont save
                    g.create_dataset(key, data=np.asarray([str(v) for v in s[key]], dtype=object), dtype=h5py.string_dtype())
            g.attrs["fs"] = float(s["fs"])
            g.attrs["start_time"] = float(s["start_time"])
            g.attrs["device"] = str(s.get("device", "unknown")) #which Synapse gizmo / rig recorded it, e.g. ecogdbs, emg, microphone
            g.attrs["units"] = str(s.get("units", "V"))
            g.attrs["type"] = str(s.get("type", "unknown")) #lfp, mer, emg, audio, stim_monitor (see datatypes.py)
            g.attrs["subtype"] = str(s.get("subtype", ""))
            g.attrs["rail"] = float(s.get("rail", np.nan)) #amplifier ceiling in volts the saturation check used, NaN = unknown

    def write(self, window):
        for name, (start, x) in window.items():
            if name not in self.datasets: #nothing to write into an empty stream
                continue
            ds, held, done = self.datasets[name], self.held[name], self.done[name]
            if start != done + held.shape[1]:
                raise ValueError(f"{name}: window starts at sample {start}, expected {done + held.shape[1]} (windows have to be back to back)")
            x = np.asarray(x, dtype=np.float32)
            if held.shape[1]:
                x = np.concatenate((held, x), axis=1)
            end = done + x.shape[1]
            if end > ds.shape[1]:
                raise ValueError(f"{name}: samples run past the end of the stream ({end} > {ds.shape[1]})")
            whole = x.shape[1] if end == ds.shape[1] else (x.shape[1] // CHUNK) * CHUNK #the stream's last chunk can be short
            if whole:
                ds[:, done:done + whole] = x[:, :whole]
            self.held[name], self.done[name] = x[:, whole:].copy(), done + whole

    def close(self):
        #every sample has to have arrived; a short or long stream means a reader lost or repeated samples
        for name, ds in self.datasets.items():
            if self.done[name] != ds.shape[1] or self.held[name].shape[1]:
                raise ValueError(f"{name}: {self.done[name] + self.held[name].shape[1]} samples arrived but the stream has {ds.shape[1]}")


def write_rest(f, result, source_dir):
    """
    Writes everything but the stream samples (epocs, scalars, stim, info, flags, provenance) into an open .h5, see save_h5 for the layout.
    Input: f (the open file), result (reader dictionary plus whatever clean_block added), source_dir (raw block folder or .vhdr)
    """
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
        _write_column(stim, key, values)
    f.create_group("sidecars") #stimulator / task files found next to a brainvision recording, kept as read (see sidecars.py)
    for name, entry in result.get("sidecars", {}).items():
        g = f.create_group(f"sidecars/{name}")
        for key, values in entry["table"].items():
            _write_column(g, key, values)
        g.attrs["kind"], g.attrs["file"], g.attrs["matched"] = entry["kind"], entry["file"], np.int8(entry["matched"])
        for key, value in entry.get("match", {}).items():
            if np.isscalar(value):
                g.attrs[f"match_{key}"] = _attr_value(value)
        if "matched_code" in entry:
            g.attrs["matched_code"] = entry["matched_code"]
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


def stale_reason(out_path):
    """
    Input: out_path, a cleaned .h5
    Returns: "" if it opens and was made by the current CLEAN_VERSION (or newer), otherwise why it should be cleaned again
    """
    try:
        with h5py.File(out_path, "r") as f:
            made_by = f.attrs.get("clean_version")
    except Exception as error: #half-copied or corrupt file; the error text has the full path, so only its type is kept
        return f"file won't open ({type(error).__name__})"
    if made_by is None:
        return "made before clean_version existed"
    if int(made_by) < CLEAN_VERSION:
        return f"clean_version {int(made_by)} is older than {CLEAN_VERSION}"
    return ""


def _write_array(group, name, array):
    group.create_dataset(name, data=array)


def _write_column(group, name, values):
    #one table column: text (numpy text or python objects) as h5py strings, numbers as they are
    values = np.asarray(values)
    if values.dtype.kind in "UOS":
        group.create_dataset(name, data=np.asarray([str(v) for v in values.ravel()], dtype=object), dtype=h5py.string_dtype())
    else:
        _write_array(group, name, values)


def _attr_value(value):
    #MATLAB reads Python True/False back as the text 'TRUE', so flags are stored as 0/1 numbers instead
    if isinstance(value, (bool, np.bool_)):
        return np.int8(value)
    if isinstance(value, (int, float, np.integer, np.floating)):
        return value
    return str(value)


def _source_files(source_dir):
    #a tdt block is its whole folder; a brainvision recording is the .vhdr plus the .eeg and .vmrk it names
    if source_dir.is_dir():
        return list(source_dir.iterdir())
    if not source_dir.is_file():
        return []
    files = [source_dir]
    if source_dir.suffix.lower() == ".vhdr":
        header = source_dir.read_text(errors="ignore").splitlines()
        for key, default in (("DataFile=", ".eeg"), ("MarkerFile=", ".vmrk")):
            name = next((line.split("=", 1)[1].strip() for line in header if line.startswith(key)), source_dir.stem + default)
            files.append(source_dir.parent / name)
    return files


def code_version():
    """
    Returns: (commit, dirty): the git commit of this clean_dbs copy ("unknown" if it isn't a git checkout or git isn't there)
             and whether it has uncommitted changes (True when unknown), so a cleaned file can be traced to the exact code
    """
    repo = Path(__file__).resolve().parent.parent
    try:
        commit = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(repo), "status", "--porcelain"], capture_output=True, text=True, check=True).stdout.strip() != ""
    except (OSError, subprocess.CalledProcessError):
        commit, dirty = "unknown", True
    return commit, dirty


def _provenance(source_dir):
    source_dir = Path(source_dir)
    commit, dirty = code_version()
    versions = {}
    for package in ("numpy", "h5py", "tdt", "mne"):
        try:
            versions[f"version_{package}"] = version(package)
        except PackageNotFoundError:
            versions[f"version_{package}"] = "not installed"
    #keep the source path from the (outermost) study ID down (e.g. uh3021/stage1 surgery/clinical_dbs1), never the full path on someone's computer
    return {
        "clean_version": np.int64(CLEAN_VERSION), #clean.py re-cleans files with an older one
        "source": source_path(source_dir),
        "source_bytes": int(sum(p.stat().st_size for p in _source_files(source_dir) if p.is_file())),
        "clean_dbs_commit": commit,
        "clean_dbs_uncommitted_changes": np.int8(dirty),  #1 means the code that made this file wasn't committed yet
        "created": datetime.now().isoformat(timespec="seconds"),
        **versions,
    }

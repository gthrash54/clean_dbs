import csv
import re
from pathlib import Path
import numpy as np

#the lab's stimulator for brainvision recordings is a Multi Channel Systems STG driven by MC_Stimulus II; every sidecar below is
#either the ASCII program it played or a table the lab's MATLAB/R scripts wrote next to it. Surveyed on Box 2026-10-12:
#  MC_Stimulus ASCII (.txt / .dat): results.txt, stim_results.dat, pp_results.dat, results_itbs.dat, "stimulation protocol.dat",
#     *_stimulationProtocol.dat, erp*_events.dat, buzz_stim.dat, PNS_DBS.txt, CHR/PP "2A_4A_6A_Ramp_100fixedISI.txt", MC001.txt (ECOG) ...
#  TTL tables (MATLAB text): TTL.txt, stim_TTL.txt, pp_TTL.txt, *_TTL.txt: condition names + the TTL width (ms) that codes each one
#  stim_struct.mat / pp_struct.mat: one struct per condition (stim_amp mA, pw ms, num_pulses, TTL ms, paired, interval ms, ...)
#  per-pulse protocol tables (R): stimuli_ordered.csv, *_stimulationProtocol.csv, erp*_events.csv (protocol_time_cum s, amplitude mA, pulsewidth us, ...)
#  ttl code tables: ttl_code.csv / *_ttlCode.csv (condition -> ttl_pw1 us)
#  E-Prime task logs (.txt, utf-16): Go_NoGo, AXCPT, Simon, InstDelReach, ... one LogFrame per trial
MCS_HEADER = "Multi Channel Systems MC_Stimulus"
MATCH_TOL_S = 0.002 #intervals between markers and sidecar events must agree within 2 ms (the STG and amplifier clocks drift ~0.01%)


def read_mcs(path):
    """
    Reads an MC_Stimulus II ASCII program (the file the lab loads into the STG stimulator).
    Input: path to the .txt / .dat
    Returns: dictionary: mode (current / voltage), units (uA / mV: MC_Stimulus writes current in uA, voltage in mV),
             n_channels (stimulation channels; channel numbers above it are the sync / TTL outputs),
             format (3 = value,time x3 + repeat; 5 = pulse,value,value,time), channels {channel number: rows (n x columns float array)}
             or None if the file isn't an MC_Stimulus program
    """
    lines = Path(path).read_text(errors="ignore").splitlines()
    if not lines or not lines[0].startswith(MCS_HEADER):
        return None
    out = {"mode":"", "units":"", "n_channels":0, "format":3, "channels":{}}
    channel = None
    for line in lines[1:]:
        cells = line.replace(":", ": ").split()
        if not cells:
            continue
        if cells[0].startswith("output"):
            out["mode"] = cells[-1].lower()
            out["units"] = "uA" if out["mode"] == "current" else "mV"
        elif cells[0].startswith("format"):
            out["format"] = int(cells[-1])
        elif cells[0] == "channels:":
            out["n_channels"] = int(cells[-1])
        elif cells[0] == "channel:":
            channel = int(cells[-1])
            out["channels"][channel] = []
        elif channel is not None and re.fullmatch(r"-?\d+(\.\d+)?", cells[0]):
            out["channels"][channel].append([float(c) for c in cells])
    width = 7 if out["format"] == 3 else 4
    out["channels"] = {k:np.array([r[:width] + [0] * (width - len(r)) for r in v], dtype=np.float64).reshape(-1, width) for k, v in out["channels"].items()}
    return out


def mcs_events(mcs):
    """
    Turns an MC_Stimulus program into one row per played row (repeats expanded), per channel.
    Input: the dictionary from read_mcs
    Returns: {channel: dictionary of arrays}: onset (s from the start of the program), duration (s), amplitude (first phase with output, in mcs["units"];
             1 = high on a sync channel), width_us (that phase), amplitude2 / width2_us (the phase after it, biphasic pulses), active (1 = output, 0 = gap row)
    Format 5 rows (pulse, value, value, time) aren't documented in the files: amplitude = first nonzero value, width_us = the row's time.
    """
    out = {}
    for ch, rows in mcs["channels"].items():
        if len(rows) == 0:
            continue
        if mcs["format"] == 3: #value,time,value,time,value,time,repeat; repeat = plays of the row (0 = once)
            rows = np.repeat(rows, np.maximum(rows[:, 6], 1).astype(int), axis=0)
            dur = (rows[:, 1] + rows[:, 3] + rows[:, 5]) / 1e6
            first = np.where(rows[:, 0] != 0, 0, np.where(rows[:, 2] != 0, 2, 4)) #first phase with output
            idx = np.arange(len(rows))
            table = {"amplitude":rows[idx, first], "width_us":rows[idx, first + 1],
                     "amplitude2":np.where(first < 4, rows[idx, np.minimum(first + 2, 4)], 0.0), "width2_us":np.where(first < 4, rows[idx, np.minimum(first + 3, 5)], 0.0),
                     "active":(rows[:, [0, 2, 4]] != 0).any(axis=1).astype(np.int8)}
        else:
            dur = rows[:, 3] / 1e6
            value = np.where(rows[:, 1] != 0, rows[:, 1], rows[:, 2])
            table = {"amplitude":value, "width_us":rows[:, 3], "amplitude2":np.zeros(len(rows)), "width2_us":np.zeros(len(rows)),
                     "active":(rows[:, 1:3] != 0).any(axis=1).astype(np.int8)}
        table["onset"] = np.concatenate(([0.0], np.cumsum(dur)[:-1]))
        table["duration"] = dur
        out[ch] = table
    return out


def read_ttl_txt(path):
    """
    Reads the MATLAB-style TTL tables the lab's stim scripts write (TTL.txt, stim_TTL.txt, pp_TTL.txt, buzzy_bee_stim_TTL.txt).
    Each condition list (name={'a' 'b' ...};) is followed by the TTL widths that code it (name_ms=[0.1 0.2 ...];).
    Input: path
    Returns: dictionary of arrays: condition (text), ttl_ms (the TTL width coding it), list (which list it came from)
             plus "settings": any other name=[...] / name = word lines (e.g. amp, pulse_width_ms, stim_mode)
    """
    text = Path(path).read_text(errors="ignore")
    cells = [(m.group(1), re.findall(r"'([^']*)'", m.group(2))) for m in re.finditer(r"^\s*(\w+)\s*=\s*\{(.*?)\}\s*;?", text, re.M | re.S)]
    #a list written as rows (BI016 STIM/TTL.txt: "pulesp_ms=[0.2 0.4 ; 1.2 1.4 ; ...]") has ";" between rows, read it as one list
    nums = [(m.group(1), [float(x) for x in m.group(2).replace(",", " ").replace(";", " ").split()]) for m in re.finditer(r"^\s*(\w+)\s*=\s*\[([^\]]*)\]\s*;?", text, re.M)]
    words = {m.group(1):m.group(2).strip() for m in re.finditer(r"^\s*(\w+)\s*=\s*([A-Za-z]\w*)\s*;?\s*$", text, re.M)}
    cond, ttl, src, used = [], [], [], set()
    order = sorted([(text.find(n + "="), "cell", n, v) for n, v in cells] + [(text.find(n + "="), "num", n, v) for n, v in nums])
    for i, (_, kind, name, values) in enumerate(order):
        if kind != "cell":
            continue
        #the TTL list is the next numeric list with the same length
        nxt = next(((n, v) for _, k, n, v in order[i + 1:] if k == "num" and len(v) == len(values) and n not in used), None)
        if nxt is None:
            continue
        used.add(nxt[0])
        cond += values
        ttl += nxt[1]
        src += [name] * len(values)
    settings = {n:v for n, v in nums if n not in used}
    settings.update(words)
    return {"condition":np.array(cond, dtype=object), "ttl_ms":np.array(ttl, dtype=np.float64), "list":np.array(src, dtype=object), "settings":settings}


def read_stim_struct(path):
    """
    Reads stim_struct.mat / pp_struct.mat / *_struct.mat (one struct per stim condition, written by the lab's MATLAB stim scripts).
    Input: path
    Returns: dictionary of arrays, one entry per condition, one key per struct field (e.g. stim_amp mA, pw ms, num_pulses, TTL ms,
             paired, anodal_cathodal, interval ms); None if it can't be read
    """
    try:
        import scipy.io
        m = scipy.io.loadmat(path, squeeze_me=True, struct_as_record=False)
    except Exception: #matlab v7.3 files are hdf5, scipy can't read them
        return _struct_v73(path)
    for key, value in m.items():
        if key.startswith("__"):
            continue
        items = np.atleast_1d(value)
        if items.size and hasattr(items.flat[0], "_fieldnames"):
            fields = items.flat[0]._fieldnames
            return {f:_column([getattr(s, f) for s in items.flat]) for f in fields}
    return None


def _struct_v73(path):
    try:
        import h5py
        with h5py.File(path, "r") as f:
            key = next(k for k in f if not k.startswith("#"))
            g = f[key]
            out = {}
            for name, ds in g.items():
                refs = ds[()]
                if refs.dtype == h5py.ref_dtype:
                    out[name] = _column([np.squeeze(f[r][()]) for r in refs.flat])
                else:
                    out[name] = np.squeeze(refs)
            return out
    except Exception:
        return None


def _column(values):
    #numbers become a float array (empty -> nan), anything else text
    try:
        return np.array([np.nan if np.size(v) == 0 else float(np.ravel(v)[0]) for v in values], dtype=np.float64)
    except (TypeError, ValueError):
        return np.array([str(v) for v in values], dtype=object)


def read_csv_table(path):
    """
    Reads a sidecar csv (stimuli_ordered.csv, *_stimulationProtocol.csv, erp*_events.csv, ttl_code.csv, *_channels.csv).
    Input: path
    Returns: dictionary of arrays, one per column; numeric columns as float (NA -> nan), the rest as text
    """
    with open(path, newline="", errors="ignore") as f:
        rows = list(csv.reader(f))
    if not rows:
        return {}
    header = [h.strip() for h in rows[0]]
    out = {}
    for i, name in enumerate(header):
        values = [r[i] if i < len(r) else "" for r in rows[1:]]
        try:
            out[name or f"col{i}"] = np.array([np.nan if v in ("", "NA", "NaN") else float(v) for v in values], dtype=np.float64)
        except ValueError:
            out[name or f"col{i}"] = np.array(values, dtype=object)
    return out


def read_eprime(path):
    """
    Reads an E-Prime text log (one LogFrame per trial / block; usually utf-16).
    Input: path
    Returns: dictionary of text arrays, one entry per LogFrame, one key per attribute seen (blank where a frame doesn't have it);
             Subject / Session / date / name fields are left out so no identifiers get copied. None if it isn't an E-Prime log
    """
    raw = Path(path).read_bytes()
    text = next((raw.decode(e, errors="ignore") for e in ("utf-16", "utf-8") if "LogFrame Start" in raw.decode(e, errors="ignore")), None)
    if text is None:
        return None
    frames = [dict((k.strip(), v.strip()) for k, v in re.findall(r"^\s*([^:\r\n]+):\s?(.*?)\r?$", fr, re.M))
              for fr in re.findall(r"\*\*\* LogFrame Start \*\*\*(.*?)\*\*\* LogFrame End \*\*\*", text, re.S)]
    skip = re.compile(r"^(Subject|Session.*|.*Name|.*Date.*|SessionTime|Clock\.Info.*|DataFile\.Basename|Handedness|Age|Sex|Group|RandomSeed|Display\.RefreshRate|Runtime.*|StudioVersion)$", re.I) #identifiers and setup, not task data
    keys = [k for k in dict.fromkeys(k for fr in frames for k in fr) if not skip.search(k)]
    return {k:np.array([fr.get(k, "") for fr in frames], dtype=object) for k in keys}


def read_sidecar(path):
    """
    Reads any sidecar this module knows.
    Input: path
    Returns: (kind, table) where kind is mcs, mcs_binary (found, can't be read), ttl_txt, stim_struct, protocol_csv, ttl_code_csv, channels_csv, eprime or None;
             table is a dictionary of equal-length arrays (for mcs: the stimulation channel's rows, see mcs_table)
    """
    path = Path(path)
    name, ext = path.name.lower(), path.suffix.lower()
    if ext == ".mat":
        t = read_stim_struct(path)
        return ("stim_struct", t) if t else (None, None)
    if ext == ".csv":
        t = read_csv_table(path)
        if "protocol_time_cum" in t:
            return "protocol_csv", t
        if "ttl_pw1" in t:
            return "ttl_code_csv", t
        if "channel" in t and "name" in t:
            return "channels_csv", t
        return None, None
    if ext in (".txt", ".dat"):
        head = path.read_bytes()[:200]
        if head.startswith(MCS_HEADER.encode()):
            mcs = read_mcs(path)
            return "mcs", mcs_table(mcs)
        if b"MC_Stimulus" in head[:40]: #MC_Stimulus's own binary save (MC012 results.txt), not the ASCII export: listed, not read
            return "mcs_binary", {"note":np.array(["binary MC_Stimulus II file, not the ASCII export; open it in MC_Stimulus and export it as ASCII to read it"], dtype=object)}
        if b"\xff\xfe" == head[:2] or b"LogFrame" in path.read_bytes()[:4000]:
            t = read_eprime(path)
            return ("eprime", t) if t else (None, None)
        if re.search(r"ttl", name):
            t = read_ttl_txt(path)
            return ("ttl_txt", {k:v for k, v in t.items() if k != "settings"}) if len(t["condition"]) else (None, None)
    return None, None


def mcs_table(mcs):
    """
    One table for an MC_Stimulus program: one row per sync (TTL) output, with what the stimulation channels started at that moment.
    Without a sync channel: one row per active row of the stimulation channels.
    Input: the dictionary from read_mcs
    Returns: dictionary of arrays: onset (s from program start), channel (stimulation channel, -1 = none started then), amplitude, width_us,
             amplitude2, width2_us (0 V / 0 uA rows are kept: sham conditions still send their TTL), ttl_ms (sync width, nan without sync),
             units (uA or mV) as a one-element array
    """
    events = mcs_events(mcs)
    n_stim = mcs["n_channels"] or max(events, default=0)
    sync = next((e for ch, e in sorted(events.items()) if ch > n_stim and e["active"].any()), None) #STG sync outputs come after the stim channels
    stim = [(ch, e) for ch, e in sorted(events.items()) if ch <= n_stim]
    if sync is not None:
        on = sync["active"] == 1
        t, ttl = sync["onset"][on], sync["width_us"][on] / 1000
    else:
        t = np.unique(np.concatenate([e["onset"][e["active"] == 1] for _, e in stim])) if stim else np.array([])
        ttl = np.full(len(t), np.nan)
    out = {"onset":t, "channel":np.full(len(t), -1, dtype=np.int64), "ttl_ms":ttl}
    for key in ("amplitude", "width_us", "amplitude2", "width2_us"):
        out[key] = np.full(len(t), np.nan)
    for ch, e in reversed(stim): #lowest channel wins when several start together
        i = np.searchsorted(e["onset"], t - 1e-7)
        ok = (i < len(e["onset"])) & (np.abs(e["onset"][np.clip(i, 0, len(e["onset"]) - 1)] - t) < 1e-6)
        ok &= (e["active"][np.clip(i, 0, len(e["onset"]) - 1)] == 1) | (out["channel"] < 0)
        for key in ("amplitude", "width_us", "amplitude2", "width2_us"):
            out[key][ok] = e[key][i[ok]]
        out["channel"][ok] = ch
    out["units"] = np.array([mcs["units"]], dtype=object)
    return out


def event_times(kind, table):
    #the times (s) a sidecar says something happened, used to line it up with the markers
    if kind == "mcs":
        return table["onset"]
    if kind == "protocol_csv":
        t = table["protocol_time_cum"]
        return t[np.isfinite(t)]
    return None


def match_times(markers, events, tol=MATCH_TOL_S):
    """
    Lines up sidecar event times with marker times by their spacing (the sidecar's clock starts when the stimulator starts, not the recording).
    Walks both lists together: each event is predicted from the last matched pair, so slow clock drift and a few missing markers or
    events (a recording stopped early, or started after the program) don't break it.
    Input: markers (s, recording clock, use the gap-corrected clock for brainvision), events (s, sidecar clock), tol (s) allowed error
    Returns: dictionary: matched (1/0), pairs (n x 2 int array: marker index, event index), n_markers, n_events, n_pairs, max_err_ms
             matched = 1 when at least 90% of the markers inside the matched stretch pair up, there are 3 or more pairs, and the pairs
             cover at least half of the shorter of the two lists
    """
    m, e = np.asarray(markers, float), np.asarray(events, float)
    empty = {"matched":0, "pairs":np.zeros((0, 2), dtype=np.int64), "n_markers":len(m), "n_events":len(e), "n_pairs":0, "max_err_ms":np.nan}
    if len(m) < 3 or len(e) < 3:
        return empty
    best = empty
    #start candidates: marker 0..4 against any event whose next two intervals agree
    dm, de = np.diff(m), np.diff(e)
    for i0 in range(min(5, len(m) - 2)):
        cands = np.flatnonzero((np.abs(de[:-1] - dm[i0]) < tol) & (np.abs(de[1:] - dm[i0 + 1]) < tol))[:50]
        for j0 in cands:
            pairs, err = [(i0, j0)], []
            for j in range(j0 + 1, len(e)):
                lm, le = pairs[-1]
                guess = m[lm] + (e[j] - e[le])
                k = np.searchsorted(m, guess)
                k = k if k < len(m) and (k == 0 or abs(m[k] - guess) < abs(m[k - 1] - guess)) else k - 1
                if k > lm and abs(m[k] - guess) <= tol:
                    pairs.append((k, j))
                    err.append(abs(m[k] - guess))
                if guess > m[-1] + tol or (j - j0 == 200 and len(pairs) < 50): #past the recording, or clearly the wrong start
                    break
            if len(pairs) > best["n_pairs"]:
                span = pairs[-1][0] - pairs[0][0] + 1
                #and the pairs have to cover at least half of the shorter list: a few chance pairs out of thousands (TM001 pp: 4 of 5165
                #markers) aren't a match
                enough = len(pairs) >= 0.5 * min(len(m), len(e))
                best = {"matched":int(len(pairs) >= 3 and len(pairs) >= 0.9 * span and enough), "pairs":np.array(pairs, dtype=np.int64),
                        "n_markers":len(m), "n_events":len(e), "n_pairs":len(pairs), "max_err_ms":float(max(err, default=0) * 1e3)}
        if best["matched"] and best["n_pairs"] >= 0.9 * len(m):
            break
    return best


#names that are sidecars, checked on the Box survey; anything else next to a recording is left alone
SIDECAR_PATTERN = re.compile(r"(ttl|results|stim|protocol|stimuli_ordered|_events|_channels|struct|buzz|ramp|pns|reversal|^mc\d+|go_?nogo|axcpt|simon|reach|choice|exp\d|prac\d)", re.I)


def _is_mcs(path):
    #an MC_Stimulus program says so on its first line, whatever it's called (CHR "0-1_80us_M.txt", PT "pulse_directionality_ch1.txt")
    if path.suffix.lower() not in (".txt", ".dat"):
        return False
    with open(path, "rb") as f:
        return f.read(len(MCS_HEADER)) == MCS_HEADER.encode()


def find_sidecars(vhdr):
    """
    Input: a brainvision .vhdr
    Returns: sidecar files in its folder and in stim or task folders next to it or one level up ("Stim", "stim files", "Task Files": ug3001 keeps
             buzz_stim.dat in <study>/stim files), small ones only; files named after the recording (erp1_events.dat for erp1.vhdr) come first
    """
    vhdr = Path(vhdr)
    up = [d for d in vhdr.parent.parent.iterdir() if d.is_dir() and re.search(r"stim|task", d.name, re.I)] if vhdr.parent.parent != vhdr.parent else []
    folders = [vhdr.parent] + [d for d in vhdr.parent.iterdir() if d.is_dir() and re.search(r"stim|task", d.name, re.I)] + up
    out = []
    for d in dict.fromkeys(folders):
        for f in sorted(d.iterdir()):
            if (f.is_file() and f.suffix.lower() in (".txt", ".dat", ".csv", ".mat") and f.stat().st_size < 50_000_000
                    and (SIDECAR_PATTERN.search(f.stem) or _is_mcs(f))):
                out.append(f)
    stem = vhdr.stem.lower()
    return sorted(out, key=lambda f: not f.stem.lower().startswith(stem)) #stable: named-after-the-recording first, then folder order



def attach_sidecars(stim, vhdr):
    """
    Reads every sidecar next to a brainvision recording, and where one lines up with the recording's stim markers, copies its
    per-pulse settings into the stim table. Sidecars that don't line up are kept as they are (nothing is dropped).
    Input:
    stim: the stim table from stim_events (needs time); changed in place when a sidecar matches
    vhdr: the .vhdr
    Returns: {sidecar name: {"kind", "file", "matched" (1/0), "match" (see match_times), "table"}} to save under /sidecars/<name>
    Matched columns added to stim: amplitude (uA for current mode, mV for voltage), width_us, ttl_ms, condition / polarity / frequency etc. from
      protocol csvs, and sidecar (file name). TTL tables (TTL.txt, stim_struct.mat) turn a ttl width into a condition name when the stim table has
      marker_lag_ms (the TTL width seen by the amplifier).
    """
    out, copied = {}, set()
    times = np.asarray(stim.get("time", []), dtype=np.float64)
    n = len(times)
    for f in find_sidecars(vhdr):
        try:
            kind, table = read_sidecar(f)
        except Exception as err: #a broken sidecar shouldn't stop cleaning; say why in the file
            kind, table = "unreadable", {"error":np.array([str(err)[:200]], dtype=object)}
        if kind is None:
            continue
        entry = {"kind":kind, "file":f.name, "matched":0, "match":{}, "table":table}
        t = event_times(kind, table) if kind in ("mcs", "protocol_csv") else None
        if t is not None and n >= 3:
            entry["match"] = match_times(times, t)
            if not entry["match"]["matched"] and "source" in stim: #a program can follow one code only (TM002 Burst1: 200 S15 burst onsets
                for code in dict.fromkeys(np.asarray(stim["source"])): #among 3,300 S14 pulses), so try each stim code on its own too
                    rows = np.flatnonzero(np.asarray(stim["source"]) == code)
                    m = match_times(times[rows], t) if len(rows) >= 3 else {"matched":0}
                    if m["matched"]:
                        m["pairs"] = np.column_stack((rows[m["pairs"][:, 0]], m["pairs"][:, 1]))
                        entry["match"] = dict(m, source=str(code))
                        break
            if entry["match"]["matched"]:
                entry["matched"] = 1
                if kind not in copied: #the first match of each kind fills the stim table (the file named after the recording wins); later ones are only marked matched
                    _copy_matched(stim, kind, table, entry["match"]["pairs"], f.name)
                    copied.add(kind)
        out[_key(f, out)] = entry
    #ttl tables: turn the ttl width measured from the markers into the condition it codes
    lag = np.asarray(stim.get("marker_lag_ms", np.full(n, np.nan)), dtype=np.float64)
    decodes = {}
    for key, entry in out.items():
        if entry["kind"] in ("ttl_txt", "stim_struct") and n and np.isfinite(lag).any():
            ttl = entry["table"]["ttl_ms"] if entry["kind"] == "ttl_txt" else entry["table"].get("TTL")
            if ttl is None:
                continue
            if entry["kind"] == "ttl_txt":
                names = entry["table"]["condition"]
            else: #one text per struct row, e.g. "stim_amp=1.4 pw=0.06 num_pulses=2 ..."
                names = np.array([" ".join(f"{k}={entry['table'][k][i]:g}" if isinstance(entry["table"][k][i], float) else f"{k}={entry['table'][k][i]}"
                                           for k in entry["table"] if k != "TTL") for i in range(len(ttl))], dtype=object)
            hit = np.abs(lag[:, None] - np.asarray(ttl, float)[None, :]) <= 0.025 #ttl widths step 0.1 ms or more; amplifier resolves 0.01-0.04 ms
            one = hit.sum(axis=1) == 1
            if one.sum() >= 0.9 * np.isfinite(lag).sum() and one.any():
                decoded = np.full(n, "", dtype=object)
                decoded[one] = names[hit[one].argmax(axis=1)]
                decodes[key] = decoded
                entry["match"] = {"matched":1, "n_markers":int(np.isfinite(lag).sum()), "n_events":len(ttl), "coded":int(one.sum())}
    #several TTL tables can decode the same widths (CHR001 has one per program, all 0.1 ms steps); if they name the same pulse
    #differently we can't tell which program ran, so no condition is filled in from any of them (never guess)
    differ = len({tuple(d) for d in decodes.values()}) > 1
    if differ: #...unless one of them was written with the program that lined up with the markers (0-1_80us_M.txt + 0-1_80us_M_TTL.txt)
        programs = {Path(e["file"]).stem.lower() for e in out.values() if e["kind"] == "mcs" and e["matched"]}
        paired = [k for k in decodes if re.sub(r"[_ -]?ttl$", "", Path(out[k]["file"]).stem, flags=re.I).lower() in programs]
        if len(paired) == 1:
            for k in decodes:
                if k != paired[0]:
                    out[k]["match"].update(matched=0, ambiguous=1)
            decodes, differ = {paired[0]:decodes[paired[0]]}, False
            out[paired[0]]["match"]["chosen_by"] = "named after the stimulator program that matched the markers"
    for key, decoded in decodes.items():
        if differ:
            out[key]["match"].update(matched=0, ambiguous=1)
            continue
        cond = np.asarray(stim.get("condition", np.full(n, "", dtype=object)), dtype=object)
        cond = np.where(cond == "", decoded, cond)
        stim["condition"] = cond
        out[key]["matched"] = 1
    return out


def _key(f, existing):
    key = re.sub(r"[^A-Za-z0-9_]+", "_", f.stem).strip("_") or "sidecar"
    while key in existing:
        key += "_"
    return key


def _copy_matched(stim, kind, table, pairs, name):
    #copies the sidecar's values onto the stim rows that paired up; unpaired rows get nan / ""
    n = len(stim["time"])
    rows, idx = pairs[:, 0], pairs[:, 1]
    def put(col, values):
        values = np.asarray(values)
        full = np.full(n, "", dtype=object) if values.dtype.kind in "OUS" else np.full(n, np.nan)
        full[rows] = values[idx]
        stim[col] = full
    if kind == "mcs":
        units = str(table["units"][0])
        if units == "uA":
            amp = np.asarray(stim.get("amplitude_ua", np.full(n, np.nan)), dtype=np.float64)
            amp[rows] = table["amplitude"][idx]
            stim["amplitude_ua"] = amp
        else:
            put("amplitude_mv", table["amplitude"])
        put("pulse_width_us", table["width_us"])
        put("amplitude2", table["amplitude2"])
        put("width2_us", table["width2_us"])
        put("ttl_ms", table["ttl_ms"])
        put("sidecar_channel", table["channel"].astype(np.float64))
    else:
        finite = np.isfinite(table["protocol_time_cum"])
        for col, values in table.items():
            if len(values) == len(finite):
                put(f"protocol_{col}", np.asarray(values)[finite])
    side = np.asarray(stim.get("sidecar", np.full(n, "", dtype=object)), dtype=object)
    side[rows] = np.where(side[rows] == "", name, side[rows] + ";" + name) #two sidecars can describe the same pulses (.dat and .csv)
    stim["sidecar"] = side

from pathlib import Path
import h5py
import numpy as np
import pandas as pd

#columns copied from the case tables, only for rows Garrett has checked (verified not blank)
CASE_COLUMNS = ["implant", "hemisphere", "target1", "target2"]
BLOCK_COLUMNS = ["target", "leads_in_place", "stim_contacts", "microdrive_mm", "notes"]


def cleaned_files(cleaned_root):
    """
    Input:
    cleaned_root: top of the cleaned set (out_root from clean.py)
    Returns: sorted list of finished .h5 files, skipping _cases (and any folder starting with _) and half-written .part files
    """
    cleaned_root = Path(cleaned_root)
    return sorted(p for p in cleaned_root.rglob("*.h5") if not any(part.startswith("_") for part in p.relative_to(cleaned_root).parts))


def load_cases(cleaned_root):
    """
    Input:
    cleaned_root: top of the cleaned set, the tables are in <cleaned_root>/_cases/cases.csv and blocks.csv
    Returns: (cases, blocks) tables with only the verified rows; empty tables if a file isn't there
    Everything is read as text so "2-3-4" or "11.910" come back exactly as typed
    """
    tables = []
    for name in ("cases.csv", "blocks.csv"):
        path = Path(cleaned_root) / "_cases" / name
        t = pd.read_csv(path, dtype=str, keep_default_na=False) if path.is_file() else pd.DataFrame(columns=["study_id", "block", "verified"])
        tables.append(t[t["verified"].str.strip() != ""] if "verified" in t else t.iloc[0:0])
    return tables[0], tables[1]


def build_index(cleaned_root):
    """
    One row per cleaned .h5 file, built from the file's attributes and array shapes only (no signal is loaded), so it's fast on hundreds of files.
    Input:
    cleaned_root: top of the cleaned set (out_root from clean.py)
    Returns: pandas DataFrame with
      path, file (path from cleaned_root down), patient, stage, block, format (tdt or bv), duration (s), start, size_mb
      streams: stream names; subtypes: type:subtype pairs, e.g. "lfp:ecog+lead;mer:spike"
      n_<type>: channel count per type (lfp, mer, emg, audio, stim_monitor, unknown), 0 if absent
      fs_<type>: sampling rate per type; text like "12207.03125;48828.125" if streams of one type differ
      events: event types present (video, stim, stim_params, marker, other); n_video_frames, n_stim_events
      video: where the .avi is (from the study ID down), blank if none
      saturated_s: seconds where any channel sits at the amplifier ceiling (stim_monitor not counted); saturated_s_<type> per type
      n_gaps, gap_s: real BrainVision recording gaps and the seconds missing (NaN for tdt, no gap check there)
      dnu, clean_dbs_commit, uncommitted_changes, created
      implant, hemisphere, target1, target2 (cases.csv) and target, leads_in_place, stim_contacts, microdrive_mm, block_notes (blocks.csv),
      only filled where that row is verified; case_verified / block_verified are 1 if it was
    Save it with index.to_csv(...) to open it in MATLAB (readtable) or Excel.
    """
    cleaned_root = Path(cleaned_root)
    cases, blocks = load_cases(cleaned_root)
    rows = []
    for path in cleaned_files(cleaned_root):
        row = _read_row(path, cleaned_root)
        case = cases[cases["study_id"] == row["patient"]]
        match = blocks[(blocks["study_id"] == row["patient"]) & (blocks["block"] == row["block"])]
        if "stage" in blocks: #if the table ever gets a stage column, use it so stage1/stage2 blocks with the same name don't mix
            match = match[match["stage"].str.replace(" ", "") == row["stage"]]
        row["case_verified"] = int(len(case) == 1)
        row["block_verified"] = int(len(match) == 1) #more than one matching row = ambiguous, so it's left out
        for col in CASE_COLUMNS:
            row[col] = case.iloc[0][col] if len(case) == 1 and col in case else ""
        for col in BLOCK_COLUMNS:
            row["block_notes" if col == "notes" else col] = match.iloc[0][col] if len(match) == 1 and col in match else ""
        rows.append(row)
    index = pd.DataFrame(rows)
    if len(index) == 0:
        return index
    counts = [c for c in index if c.startswith("n_") and c not in ("n_gaps",)]
    index[counts] = index[counts].fillna(0).astype(int) #types missing from a file are 0 channels, not NaN
    sats = [c for c in index if c.startswith("saturated_s")]
    index[sats] = index[sats].fillna(0.0)
    index.loc[index["format"] == "tdt", ["n_gaps", "gap_s"]] = np.nan
    #keep the per-type columns together after the stream list so the table reads left to right
    per_type = sorted(c for c in index if c.startswith(("n_", "fs_")) and c not in ("n_gaps", "n_video_frames", "n_stim_events"))
    per_sat = sorted(c for c in sats if c != "saturated_s")
    first = [c for c in index if c not in per_type + per_sat]
    at, at_sat = first.index("subtypes") + 1, first.index("saturated_s") + 1
    return index[first[:at] + per_type + first[at:at_sat] + per_sat + first[at_sat:]]


def _read_row(path, cleaned_root):
    #everything here is attributes or .shape, plus the small flag/epoc arrays; the signal arrays are never read
    rel = path.relative_to(cleaned_root)
    row = {"path":str(path), "file":rel.as_posix(), "patient":rel.parts[0], "stage":rel.parts[1], "block":path.stem}
    with h5py.File(path, "r") as f:
        info, flags = f["info"].attrs, f["flags"]
        row["format"] = "bv" if str(f.attrs.get("source", "")).lower().endswith(".vhdr") else "tdt"
        row["duration"] = float(info.get("duration", np.nan))
        row["start"] = str(info.get("start", ""))
        row["size_mb"] = round(path.stat().st_size / 1e6, 1)
        types, fs, subtypes, sat = {}, {}, [], {}
        for name, g in f["streams"].items():
            t = str(g.attrs.get("type", "unknown"))
            types[t] = types.get(t, 0) + g["data"].shape[0]
            fs.setdefault(t, set()).add(float(g.attrs["fs"]))
            if g.attrs.get("subtype", ""):
                subtypes.append(f"{t}:{g.attrs['subtype']}")
            if f"saturation_{name}" in flags:
                s = flags[f"saturation_{name}"]
                sat.setdefault(t, []).extend(zip(s["start"][:], s["end"][:]))
        row["streams"] = ";".join(f["streams"].keys())
        row["subtypes"] = ";".join(sorted(set(subtypes)))
        for t in types:
            row[f"n_{t}"] = types[t]
            row[f"fs_{t}"] = fs[t].pop() if len(fs[t]) == 1 else ";".join(str(x) for x in sorted(fs[t]))
        events = {}
        for group in ("epocs", "scalars"):
            for name, g in f.get(group, {}).items():
                t = str(g.attrs.get("type", "other"))
                events[t] = events.get(t, 0) + (g["onset"].shape[0] if "onset" in g else g["ts"].shape[0])
        row["events"] = ";".join(sorted(events))
        row["n_video_frames"] = events.get("video", 0)
        row["n_stim_events"] = events.get("stim", 0)
        row["video"] = str(info.get("video", ""))
        #the stim monitor is left out of the total: it's the stimulator's output, and its unused channels sit constant (which find_saturation counts as saturated)
        row["saturated_s"] = _union_seconds([iv for t, ivs in sat.items() if t != "stim_monitor" for iv in ivs])
        for t in sat:
            row[f"saturated_s_{t}"] = _union_seconds(sat[t])
        if "new_segments" in flags:
            gap = flags["new_segments"]["is_gap"][:].astype(bool)
            row["n_gaps"] = int(gap.sum())
            row["gap_s"] = float(flags["new_segments"]["missing"][:][gap].sum())
        else:
            row["n_gaps"], row["gap_s"] = 0, 0.0
        row["dnu"] = int(flags.attrs.get("dnu", 0))
        row["clean_dbs_commit"] = str(f.attrs.get("clean_dbs_commit", ""))
        row["uncommitted_changes"] = int(f.attrs.get("clean_dbs_uncommitted_changes", 1))
        row["created"] = str(f.attrs.get("created", ""))
    return row


def _union_seconds(intervals):
    #total time covered by (start, end) stretches, counting overlaps once (two channels saturating together is still one second)
    total, last_end = 0.0, -np.inf
    for start, end in sorted(intervals):
        if end > last_end:
            total += end - max(start, last_end)
            last_end = end
    return round(float(total), 4)


def find_clean(index, type=None, subtype=None, max_saturated_s=None, max_gap_s=None, min_duration=None, include_dnu=False, **columns):
    """
    Picks cleaned files out of an index from build_index.
    Input:
    index: the table from build_index
    type: only files with channels of this type, e.g. "lfp", "mer", "emg"
    subtype: only files with this subtype, e.g. "spike" or "lead"
    max_saturated_s: at most this many saturated seconds (of that type if type is given, otherwise any stream); 0 = none
    max_gap_s: at most this many seconds missing in BrainVision gaps (tdt files always pass)
    min_duration: at least this many seconds long
    include_dnu: False drops blocks marked do not use
    **columns: any other index column must equal the value (text is not case sensitive), or be one of a list,
        e.g. target="STN", implant="dual", patient=["uh3021", "uh3003"], block_verified=1
        case/block columns are blank unless verified, so these filters only ever match verified rows
    Returns: list of .h5 paths
    Example: find_clean(idx, type="lfp", target="STN", implant="dual", max_saturated_s=0)
    """
    keep = pd.Series(True, index=index.index)
    if type is not None:
        keep &= index.get(f"n_{type}", 0) > 0
    if subtype is not None:
        keep &= index["subtypes"].apply(lambda s: subtype in [p.split(":", 1)[1] for p in s.split(";") if ":" in p])
    if max_saturated_s is not None:
        column = f"saturated_s_{type}" if type is not None else "saturated_s"
        keep &= index.get(column, 0) <= max_saturated_s
    if max_gap_s is not None:
        keep &= index["gap_s"].fillna(0) <= max_gap_s
    if min_duration is not None:
        keep &= index["duration"] >= min_duration
    if not include_dnu:
        keep &= index["dnu"] == 0
    for column, wanted in columns.items():
        if column not in index:
            raise ValueError(f"{column!r} isn't a column of the index; columns are {list(index.columns)}")
        wanted = wanted if isinstance(wanted, (list, tuple, set)) else [wanted]
        wanted = [str(w).lower() for w in wanted]
        keep &= index[column].astype(str).str.lower().isin(wanted)
    return [Path(p) for p in index.loc[keep, "path"]]


#TDT rates are the 195312.5 Hz base clock divided by a power of 2 (48828.125, 24414.0625, 12207.03125, ...)
TDT_RATES = [195312.5 / 2 ** k for k in range(12)]


def qc_block(h5_path, sample_s=1.0):
    """
    Quick checks on one cleaned file. Reads at most 3 x sample_s seconds per stream (start, middle, end), so it stays fast.
    Input:
    h5_path: a cleaned .h5 at <cleaned_root>/<patient>/<stage>/<block>.h5 (the case tables are found from there)
    sample_s: seconds read at each of the three spots
    Returns: list of issues (text), empty list = nothing found
    Checks: empty streams, NaN/inf in the sample, channels that are perfectly flat in the sample, tdt sampling rates
    that aren't a TDT rate, streams much shorter or longer than the recording, events outside the recording,
    Cam1 frames with no .avi (or the other way round), dnu, and lead streams vs leads_in_place in a verified blocks.csv row
    """
    h5_path = Path(h5_path)
    issues = []
    with h5py.File(h5_path, "r") as f:
        is_bv = str(f.attrs.get("source", "")).lower().endswith(".vhdr")
        duration = float(f["info"].attrs.get("duration", np.nan))
        lead_streams = 0
        for name, g in f["streams"].items():
            data, fs = g["data"], float(g.attrs["fs"])
            n_ch, n = data.shape
            if n_ch == 0 or n == 0:
                issues.append(f"{name}: empty stream ({n_ch} channels x {n} samples)")
                continue
            if g.attrs.get("type") == "lfp" and g.attrs.get("subtype") in ("lead", "ecog+lead", ""):
                lead_streams += 1 #brainvision ("") is the lead in stage 2
            if not is_bv and not any(np.isclose(fs, r) for r in TDT_RATES):
                issues.append(f"{name}: fs {fs} Hz is not a TDT rate")
            if abs(n / fs - duration) > 1: #one stream can end a few ms early, more than a second means something got cut
                issues.append(f"{name}: {n / fs:.2f} s long but the recording is {duration:.2f} s")
            w = max(1, int(sample_s * fs))
            spots = sorted({0, max(0, n // 2 - w // 2), max(0, n - w)})
            sample = np.concatenate([data[:, s:s + w] for s in spots], axis=1)
            bad = np.flatnonzero(~np.isfinite(sample).all(axis=1))
            if bad.size:
                issues.append(f"{name}: NaN or inf on channels {g['channels'][:][bad].tolist()} in the sample")
            flat = np.flatnonzero(np.ptp(np.nan_to_num(sample), axis=1) == 0)
            if flat.size:
                issues.append(f"{name}: channels {g['channels'][:][flat].tolist()} are perfectly flat in the sample")
        for group in ("epocs", "scalars"):
            for name, g in f.get(group, {}).items():
                t = g["onset"][:] if "onset" in g else g["ts"][:]
                if t.size and (not np.isfinite(t).all() or t.min() < -1 or t.max() > duration + 1):
                    issues.append(f"{name}: event times run outside the recording (0 to {duration:.2f} s)")
        has_cam = "epocs/Cam1" in f and f["epocs/Cam1/onset"].shape[0] > 0
        video = str(f["info"].attrs.get("video", ""))
        if has_cam and not video:
            issues.append("Cam1 frames but no .avi was found next to the block")
        if video and not has_cam:
            issues.append("an .avi is listed but there are no Cam1 frames to line it up")
        if int(f["flags"].attrs.get("dnu", 0)):
            issues.append("block is marked dnu (do not use)")
    #compare with the verified OR data sheet row, when there is one
    root, patient, block = h5_path.parents[2], h5_path.parents[1].name, h5_path.stem
    _, blocks = load_cases(root)
    match = blocks[(blocks["study_id"] == patient) & (blocks["block"] == block)]
    if len(match) == 1 and match.iloc[0].get("leads_in_place", "").strip().isdigit():
        leads = int(match.iloc[0]["leads_in_place"])
        if leads > 0 and lead_streams == 0:
            issues.append(f"blocks.csv says {leads} lead(s) in place but the file has no lead stream")
    return issues


def qc_report(cleaned_root, sample_s=1.0):
    """
    Runs qc_block on every cleaned file, and lists leftover .part files (a cleaning job that crashed mid-write).
    Input:
    cleaned_root: top of the cleaned set
    sample_s: seconds read at each of the three spots per stream, see qc_block
    Returns: pandas DataFrame, one row per file: file, patient, stage, block, n_issues, issues (joined with " | ")
    """
    cleaned_root = Path(cleaned_root)
    rows = []
    for path in cleaned_files(cleaned_root):
        rel = path.relative_to(cleaned_root)
        try:
            issues = qc_block(path, sample_s)
        except Exception as error: #a file that won't even open is the biggest issue of all
            issues = [f"could not read: {error}"]
        rows.append({"file":rel.as_posix(), "patient":rel.parts[0], "stage":rel.parts[1], "block":path.stem, "n_issues":len(issues), "issues":" | ".join(issues)})
    for path in sorted(cleaned_root.rglob("*.h5.part")):
        rel = path.relative_to(cleaned_root)
        rows.append({"file":rel.as_posix(), "patient":rel.parts[0], "stage":rel.parts[1] if len(rel.parts) > 2 else "", "block":path.name[:-len(".h5.part")],
                     "n_issues":1, "issues":"leftover .part file, cleaning crashed or is still running; delete it and clean again"})
    return pd.DataFrame(rows, columns=["file", "patient", "stage", "block", "n_issues", "issues"])

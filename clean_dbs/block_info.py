import re
from pathlib import Path
from clean_dbs.cohorts import COHORTS, ID_IN_NAME_RE, NUMBER_RE, cohort_of, condition_of
from clean_dbs import privacy

#a Box number folder on its own (digits and dots, sometimes a short suffix); never written into a source path or report
NUMBER_FOLDER_RE = re.compile(r"^\d+(\.\d+)*(_[A-Za-z]+)?$")


def describe_block(block):
    """
    Input:
    block: a tdt block folder or a brainvision .vhdr file (one item from find_all_blocks)
    Returns: dictionary with patient, stage, block (name), format ("tdt" or "bv")
        uh/UG (and any recording in a stage folder with no study ID above it) give just those four, exactly as before
        every other cohort also gives cohort, session (subfolder path below the ID, lowercase, spaces as _, "" at the ID root),
        parent_id (the outer ID when this ID's folder sits inside another one, else ""), condition and state (see cohorts.condition_of);
        stage is "stageN" there only if a stage folder is in the session, else ""
    Folder and file names go through <out>/_cases/renames.csv first (see privacy.py), so a name written there never reaches the outputs
    Raises ValueError when no patient can be worked out (no stage folder for uh/UG, or no study ID above it or in its name)
    """
    block = Path(block)
    is_bv = block.suffix.lower() == ".vhdr"
    name, fmt = block.stem if is_bv else block.name, "bv" if is_bv else "tdt"
    folders = block.parent.parts
    ids = [i for i, part in enumerate(folders) if cohort_of(part)]
    #in a number folder with no ID folder (straight in it, or in a cohort folder with no number like MOTIF/): only the file name can say whose it is
    if not ids and any(NUMBER_RE.match(f) for f in folders) and not re.search(r"stage\s*\d", block.parent.name, re.IGNORECASE):
        return _from_name(block, name, fmt)
    #uh/UG keep the original stage layout and output path, so cleaned sets made before cohorts existed still line up
    if not ids or COHORTS[cohort_of(folders[ids[-1]])]["layout"] == "stage":
        return _stage_layout(block, name, fmt)
    #innermost ID wins, so LFP03/CHR003/x.vhdr belongs to CHR003, with LFP03 kept as parent_id. The same ID again further down
    #(DLPFC110/Contents of DLPFC110 folder.../DLPFC110/x.vhdr) is the same person, so the outermost copy of the innermost ID starts the session
    people = [canonical_id(folders[i]) for i in ids]
    first = ids[people.index(people[-1])]
    below = tuple(privacy.rename(f, people[-1]) for f in folders[first + 1:])
    name = privacy.rename(name, people[-1])
    stage = next((f"stage{m.group(1)}" for f in below if (m := re.search(r"stage\s*(\d)", f, re.IGNORECASE))), "")
    cohort = cohort_of(folders[ids[-1]])
    parents = [p for p in people if p != people[-1]]
    condition, state = condition_of(_without_ids("/".join(below + (name,)), people), cohort)
    return {"patient":people[-1], "stage":stage, "block":name, "format":fmt, "cohort":cohort,
            "session":"/".join(f.strip().lower().replace(" ", "_") for f in below),
            "parent_id":parents[-1] if parents else "", "condition":condition, "state":state}


def _without_ids(text, people):
    #a study ID in a folder or file name (BPP032.vhdr, LFP/BI020.vhdr) names the person, not what was recorded, so "BPP" isn't read as paired pulses
    return ID_IN_NAME_RE.sub(lambda m: "" if canonical_id(m.group(0)) in people else m.group(0), text)


def canonical_id(name):
    """
    Input:
    name: a study ID as written on a folder or file ("dlpfc045", "SEP_03", "Burst012")
    Returns: one spelling per person, prefix in capitals and no separator ("DLPFC045", "SEP03", "BURST012"); uh/UG IDs never go through this
    """
    match = ID_IN_NAME_RE.search(str(name))
    return match.group(1).upper() + match.group(2) if match else str(name)


def _stage_layout(block, name, fmt):
    #folder names vary ("stage1 surgery", "stage 2 surgery", "stage2 surgery (sedated, not GA)"), so pull out just "stage" + the number
    match = re.search(r"stage\s*(\d)", block.parent.name, re.IGNORECASE)
    if match is None:
        raise ValueError(f"Can't tell the stage from {block.parent.name!r}, expected a folder like 'stage1 surgery' above {block.name}")
    patient = block.parent.parent.name
    return {"patient":patient, "stage":"stage" + match.group(1), "block":privacy.rename(name, canonical_id(patient)), "format":fmt}


def _from_name(block, name, fmt):
    #a few recordings sit straight in a number folder; the ID is only usable if the file name carries one
    match = ID_IN_NAME_RE.search(name)
    if match is None:
        raise ValueError(f"No study ID folder above {block.name} and none in its name, can't tell whose recording it is")
    cohort = match.group(1).lower()
    name = privacy.rename(name, canonical_id(match.group(0)))
    condition, state = condition_of(_without_ids(name, [canonical_id(match.group(0))]), cohort)
    return {"patient":canonical_id(match.group(0)), "stage":"", "block":name, "format":fmt, "cohort":cohort, "session":"",
            "parent_id":"", "condition":condition, "state":state}


def recording_condition(block):
    """
    Input:
    block: a tdt block folder or a brainvision .vhdr file, from any cohort including uh/UG
    Returns: (condition, state) read from the folder and file names, see cohorts.condition_of; ("", "") if they don't say
    """
    try:
        info = describe_block(block)
    except ValueError:
        return condition_of(Path(block).stem)
    if "condition" in info:
        return info["condition"], info["state"]
    #uh/UG: the stage folder isn't a condition, so only the block name counts
    return condition_of(info["block"], cohort_of(info["patient"]) or "")


def source_path(block):
    """
    Input:
    block: a tdt block folder or a brainvision .vhdr file
    Returns: the raw path from the outermost study ID down (e.g. "uh3021/stage1 surgery/clinical_dbs1", "PP012/LFP034/lfp1.vhdr"),
        never the full path on someone's computer; uh/UG keep the last 3 parts as before. Without a study ID folder above it, the path
        below the Box number folder (at most 3 parts), so a number folder is never in it. Names go through renames.csv like describe_block
    """
    parts = Path(block).parts
    out, scope = [], ""
    for i in range(_source_start(parts), len(parts)):
        if cohort_of(parts[i]):
            scope = canonical_id(parts[i])
            out.append(parts[i])
        elif i == len(parts) - 1 and Path(block).suffix.lower() == ".vhdr":
            out.append(privacy.rename_file(parts[i], scope or _id_in_name(parts[i])))
        else:
            out.append(privacy.rename(parts[i], scope or _id_in_name(parts[i])))
    return "/".join(out)


def _source_start(parts):
    #where source_path starts: the outermost study ID; uh/UG and recordings with no ID folder: the last 3 parts, but below any number folder
    ids = [i for i, part in enumerate(parts[:-1]) if cohort_of(part)]
    if ids and COHORTS[cohort_of(parts[ids[-1]])]["layout"] != "stage":
        return ids[0]
    start = max(len(parts) - 3, 0)
    numbers = [i for i in range(start, len(parts) - 1) if NUMBER_FOLDER_RE.match(parts[i])]
    return numbers[-1] + 1 if numbers else start


def has_study_id(block):
    """
    Input: a tdt block folder or a brainvision .vhdr file
    Returns: True if a study ID is in a folder above it or in its own name (False: MOTIF/x.vhdr, PR/polarity02.vhdr, can't be cleaned)
    """
    block = Path(block)
    return any(cohort_of(part) for part in block.parts[:-1]) or ID_IN_NAME_RE.search(block.name) is not None


def _id_in_name(name):
    match = ID_IN_NAME_RE.search(str(name))
    return canonical_id(match.group(0)) if match else ""


def names_to_check(block):
    """
    Every name of this recording that would go into the cleaned set (folders from the study ID down, its own name, the signal and
    marker files its header names, the stimulator/task files cleaned with it, tdt videos), after renames.csv
    Input:
    block: a tdt block folder or a brainvision .vhdr file
    Returns: list of {"study_id", "kind", "original"} for the names that don't pass privacy.is_clear (empty list = nothing to check)
    """
    from clean_dbs.sidecars import find_sidecars #here, not at the top: sidecars isn't needed to describe a block
    block = Path(block)
    parts, is_bv = block.parts, block.suffix.lower() == ".vhdr"
    names, scope = [], ""
    for i in range(_source_start(parts), len(parts)): #the same parts as source_path
        if cohort_of(parts[i]):
            scope = canonical_id(parts[i])
            continue
        last = i == len(parts) - 1
        names.append((scope or _id_in_name(parts[i]), "recording" if last else "folder", block.stem if last and is_bv else parts[i]))
    scope = scope or _id_in_name(block.name)
    if is_bv and block.is_file():
        header = block.read_text(errors="ignore").splitlines()
        for key in ("DataFile=", "MarkerFile="):
            name = next((line.split("=", 1)[1].strip() for line in header if line.startswith(key)), "")
            if name:
                names.append((scope, "data file", Path(name).stem))
        names += [(scope, "sidecar", f.stem) for f in find_sidecars(block)]
    elif block.is_dir():
        names += [(scope, "video", v.stem) for v in sorted(block.glob("*.avi"))]
    out, seen = [], set()
    for study_id, kind, original in names:
        if (study_id, original.lower()) not in seen and not privacy.is_clear(privacy.rename(original, study_id)):
            out.append({"study_id":study_id, "kind":kind, "original":original})
        seen.add((study_id, original.lower()))
    return out


def missing_files(block):
    """
    Input:
    block: a tdt block folder or a brainvision .vhdr file
    Returns: list of the signal files that should be there but aren't (empty list = everything is there)
    """
    block = Path(block)
    if block.suffix.lower() == ".vhdr":
        #the .vhdr names its signal file on the "DataFile=" line
        data_file = next((line.split("=", 1)[1].strip() for line in block.read_text(errors="ignore").splitlines() if line.startswith("DataFile=")), block.stem + ".eeg")
        path = block.parent / data_file
        return [] if path.is_file() and path.stat().st_size > 0 else [data_file] #a 0-byte .eeg (DLPFC085 burstpp_trial2) has no data either
    #tdt keeps the signal in .tev (or .sev files on some rigs)
    return [] if any(block.glob("*.tev")) or any(block.glob("*.sev")) else [".tev"]


def output_path(block, out_root):
    """
    Input:
    block: a tdt block folder or a brainvision .vhdr file
    out_root: top of the cleaned set
    Returns: where that block's cleaned file goes
        uh/UG: out_root/<patient>/<stage>/<block>.h5 (unchanged)
        other cohorts: out_root/<patient>/<session>/<block>.h5, the session keeps LFP1.vhdr in two subfolders apart
    """
    info = describe_block(block)
    if "session" not in info:
        return Path(out_root) / info["patient"] / info["stage"] / (info["block"] + ".h5")
    return Path(out_root) / info["patient"] / Path(*info["session"].split("/")) / (info["block"] + ".h5")


def recording_size(block):
    """
    Input:
    block: a tdt block folder or a brainvision .vhdr file
    Returns: bytes of the signal (the .eeg named in the .vhdr, or the .tev/.sev/.tsq files), or of the .vhdr when its signal isn't there
    """
    block = Path(block)
    if block.is_dir():
        return sum(f.stat().st_size for f in block.iterdir() if f.suffix.lower() in (".tev", ".sev", ".tsq"))
    names = [line.split("=", 1)[1].strip() for line in block.read_text(errors="ignore").splitlines() if line.startswith("DataFile=")]
    data = block.parent / (names[0] if names else block.stem + ".eeg")
    return data.stat().st_size if data.is_file() else block.stat().st_size


def find_duplicates(blocks, size=recording_size):
    """
    The same recording is often copied under several number folders on Box, so two recordings can land on the same output path.
    Input:
    blocks: list of recordings (from find_all_blocks)
    size: function giving a recording's size, to tell a true copy from a different recording with the same name and place
    Returns: dictionary {skipped recording: (the one that's kept, reason)}. reason is "duplicate" (same size) or "conflict" (sizes differ;
        the biggest is kept since a smaller one is usually an incomplete copy, check these by hand). Equal sizes keep the first by path,
        so the choice is the same every run. Recordings that can't be described (no ID) are left out, they fail later with the reason.
    """
    groups = {}
    for block in sorted(blocks):
        try:
            where = str(output_path(block, "")).lower() #lowercase so it's a clash on Mac/Windows too
        except ValueError:
            continue
        groups.setdefault(where, []).append((block, size(block)))
    skipped = {}
    for copies in groups.values():
        kept, kept_size = max(copies, key=lambda c: c[1]) #max returns the first of equal sizes, and copies are sorted by path
        for block, n in copies:
            if block != kept:
                skipped[block] = (kept, "duplicate" if n == kept_size else "conflict")
    return skipped

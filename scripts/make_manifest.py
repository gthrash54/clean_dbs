import sys
import argparse
from pathlib import Path
import pandas as pd
from clean_dbs import find_all_blocks, describe_block, missing_files, privacy
from clean_dbs.block_info import find_duplicates, recording_condition, source_path, output_path, recording_size, names_to_check, has_study_id

#a cluster job array needs a fixed numbered list of blocks, so task 7 always means the same block even if new data shows up mid-run
parser = argparse.ArgumentParser(description="Write one row per recording (path, patient, stage, block, format, cohort, session, condition, missing, duplicate, excluded) for clean.py --manifest")
parser.add_argument("raw_folders", nargs="+", help="one or more raw folders (whole inventory, one cohort's copy, one patient, or one surgery folder)")
parser.add_argument("manifest", help="csv to write, keep it in the out folder (never in the repo) since it has the raw paths")
parser.add_argument("--out", help="the out folder clean.py will write to (for _cases/renames.csv and other manifests there); default: the manifest's folder")
parser.add_argument("--skip-missing", action="store_true", help="leave out blocks missing their signal files (.tev / .eeg)")
parser.add_argument("--exclude-unchecked", action="store_true", help="instead of stopping, list recordings whose names aren't cleared as excluded (not cleaned) and go on")
parser.add_argument("--chunk", type=int, default=1, help="blocks per array task, only used to print the --array range")
args = parser.parse_args()

out_root = Path(args.out) if args.out else Path(args.manifest).resolve().parent
try:
    n_renames, n_words = privacy.use_out_folder(out_root) #renames.csv and allowed_words.txt from <out>/_cases
except ValueError as error:
    sys.exit(str(error))
if n_renames or n_words:
    print(f"using {n_renames} renames and {n_words} allowed words from {out_root / privacy.CASES}")

#already sorted, so the row numbers are the same every time on the same data; a recording under two of the folders given is listed once
blocks = sorted({Path(b).resolve(): b for root in args.raw_folders for b in find_all_blocks(root)}.values())

#names first: nothing below may print or save a name that isn't cleared (see clean_dbs/privacy.py)
unchecked = {block:names_to_check(block) for block in blocks}
unchecked = {block:names for block, names in unchecked.items() if names}
if unchecked:
    todo = out_root / privacy.CASES / privacy.TO_CHECK
    rows = pd.DataFrame([{**n, "renamed":"", "example":source_path(block)} for block, names in unchecked.items() for n in names])
    rows = rows.drop_duplicates(["study_id", "original"]).sort_values(["study_id", "kind", "original"])
    todo.parent.mkdir(parents=True, exist_ok=True)
    rows[["study_id", "kind", "original", "renamed", "example"]].to_csv(todo, index=False)
    per_id = rows.groupby("study_id").size()
    print(f"{len(unchecked)} recordings have {len(rows)} folder or file names that aren't made of known words (could be a person's name): "
          + ", ".join(f"{study_id or 'no ID'} {n}" for study_id, n in per_id.items()))
    print(f"they're listed in {todo}. Fill in 'renamed' and copy those rows into {privacy.RENAMES} next to it, or add words that "
          f"describe recordings to {privacy.ALLOWED}, then run this again (docs/cheaha.md, 'Cleaning the other cohorts')")
    if not args.exclude_unchecked:
        sys.exit("stopped, nothing written (add --exclude-unchecked to leave those recordings out and go on)")

#recordings that can't be cleaned at all: no study ID above them or in their name (MOTIF/, PR/ ...), or names not cleared
excluded = {}
for block in blocks:
    if block in unchecked:
        excluded[block] = f"names not cleared ({len(unchecked[block])}), see _cases/{privacy.TO_CHECK}"
        continue
    if not has_study_id(block):
        excluded[block] = f"no study ID above it or in its name ({privacy.mask(source_path(block))}), file it under a study ID first"
blocks_ok = [b for b in blocks if b not in excluded]
#copies of one recording under several number folders all land on one output path; the extra copies stay in the manifest
#(so they show up in the report) but clean.py skips them
skipped = find_duplicates(blocks_ok)

rows = []
for block in blocks:
    if block in excluded:
        #the descriptive columns stay blank: they'd carry the names that weren't cleared, and the report copies them
        row = {"patient":"", "stage":"", "block":privacy.WITHHELD if block in unchecked else Path(block).stem, "format":"bv" if Path(block).suffix.lower() == ".vhdr" else "tdt"}
    else:
        try:
            row = describe_block(block)
        except ValueError: #bad stage folder name, kept so clean.py reports it as failed instead of it silently vanishing
            row = {"patient":"", "stage":"", "block":Path(block).stem, "format":"bv" if Path(block).suffix.lower() == ".vhdr" else "tdt"}
        row["condition"], row["state"] = recording_condition(block)
    row["excluded"] = excluded.get(block, "")
    kept, why = skipped.get(block, ("", ""))
    row["duplicate"] = why #"duplicate" = same recording as duplicate_of, "conflict" = same name and place but a different size
    row["duplicate_of"] = source_path(kept) if kept else "" #from the study ID down, the full path is in that row's path
    row["missing"] = ";".join(missing_files(block)) if block not in unchecked else ("signal file" if missing_files(block) else "") #";" not "," so it stays one csv cell
    row["path"] = str(Path(block).resolve()) #full path so the job finds it no matter what folder it starts in
    row["signal_bytes"] = recording_size(block) #the signal files (.tev/.sev/.tsq or the .eeg); exact, to tell copies apart across manifests
    row["signal_gb"] = round(row["signal_bytes"] / 1e9, 3)
    #cleaning reads one window at a time, ~0.5 GB peak whatever the size (measured); .sev blocks are read whole (raw + volts + a copy)
    row["est_peak_gb"] = round(0.5 + (3 * row["signal_gb"] if Path(block).is_dir() and any(Path(block).glob("*.sev")) else 0), 1)
    #where clean.py will write it, from the out folder down; clean.py checks it again, so a renames.csv changed afterwards can't go unnoticed
    try:
        row["out"] = output_path(block, "").as_posix() if block not in excluded else ""
    except ValueError:
        row["out"] = ""
    rows.append(row)

columns = ["path", "patient", "stage", "block", "format", "cohort", "session", "parent_id", "condition", "state", "excluded", "missing", "duplicate",
           "duplicate_of", "out", "signal_gb", "signal_bytes", "est_peak_gb"]
manifest = pd.DataFrame(rows, columns=columns).fillna("") #uh/UG rows have no cohort/session/parent_id, blank not NaN
#two recordings that aren't marked as copies of each other must never write the same .h5 (lowercase: one file on Mac/Windows)
live = manifest[(manifest["out"] != "") & (manifest["duplicate"] == "")]
outs = live["out"].str.lower()
clashes = sorted(set(outs[outs.duplicated()]))
if clashes:
    sys.exit(f"{len(clashes)} output path(s) would be written by more than one recording, fix the folder names first: {clashes[:10]}")


def other_manifests(out_root, this):
    #the other manifests in the same out folder (one per cohort, e.g. manifest.csv for uh/UG and manifest_dlpfc.csv),
    #as {out path: (raw path, (bytes, gb), manifest name)}
    seen = {}
    for m in sorted(Path(out_root).glob("manifest*.csv")):
        if m.resolve() == Path(this).resolve():
            continue
        other = pd.read_csv(m, keep_default_na=False, dtype=str)
        if "out" not in other: #made before manifests had the out column: uh/UG layout, or patient/session/block
            session = other["session"] if "session" in other else pd.Series("", index=other.index)
            cohort = other["cohort"] if "cohort" in other else pd.Series("", index=other.index)
            other["out"] = [f"{p}/{st}/{b}.h5" if not c else "/".join(x for x in (p, s, b + ".h5") if x)
                            for p, st, b, s, c in zip(other["patient"], other["stage"], other["block"], session, cohort)]
        for row in other.to_dict("records"):
            if row["out"] and not row.get("duplicate", "") and not row.get("excluded", ""):
                seen[row["out"].lower()] = (row["path"], (row.get("signal_bytes", ""), row.get("signal_gb", "")), m.name)
    return seen


#cohorts share one out folder, so a recording in this manifest must not land on a file another manifest writes: a nested ID
#copied with its parent's cohort (LFP034 inside PP012) and again on its own, or the same recording under two number folders
others = other_manifests(out_root, args.manifest)
cross = []
for i, row in manifest.iterrows():
    hit = others.get(row["out"].lower()) if row["out"] and not row["duplicate"] else None
    if hit is None or hit[0] == row["path"]: #the same raw file in a remade manifest isn't a clash
        continue
    #the same size (to the byte; to the MB for a manifest made before signal_bytes): another copy, cleaned once by the other manifest
    if (hit[1][0] and hit[1][0] == str(row["signal_bytes"])) or (not hit[1][0] and hit[1][1] and float(hit[1][1]) == row["signal_gb"]):
        manifest.loc[i, ["duplicate", "duplicate_of"]] = ["duplicate", f"{hit[2]}: {row['out']}"]
    else:
        cross.append(f"{row['out']} (also in {hit[2]})")
if cross:
    sys.exit(f"{len(cross)} output path(s) are written by another manifest in {out_root} from a different recording, sort these out first: {cross[:10]}")

n_missing = int((manifest["missing"] != "").sum())
n_dup = int((manifest["duplicate"] != "").sum())
print(f"{len(manifest)} recordings found, {n_missing} missing signal files, {n_dup} extra copies that clean.py will skip")
if excluded:
    reasons = pd.Series([r.split(" (")[0] for r in excluded.values()]).value_counts()
    print(f"{len(excluded)} excluded (listed in the manifest and the report, never cleaned): " + ", ".join(f"{n} {r}" for r, n in reasons.items()))
print(f"{manifest['signal_gb'].sum():.0f} GB of signal, biggest peak memory estimate {manifest['est_peak_gb'].max() if len(manifest) else 0} GB")
if args.skip_missing:
    manifest = manifest[manifest["missing"] == ""].reset_index(drop=True) #reset so row numbers stay 0, 1, 2... with no gaps
    print(f"skipped {n_missing}, {len(manifest)} left in the manifest")

Path(args.manifest).parent.mkdir(parents=True, exist_ok=True)
manifest.to_csv(args.manifest, index=False)
print(f"saved {args.manifest}")
#slurm array numbers are 0 to (number of tasks - 1), each task does --chunk rows
if len(manifest) > 0:
    print(f"with --chunk {args.chunk} submit with --array=0-{(len(manifest) - 1) // args.chunk}")
    #the case tables only hold uh/UG; a cohort without rows there runs with --no-cases (targets blank, case_status "missing")
    cases = out_root / "_cases" / "cases.csv"
    ids = set(manifest["patient"].str.lower()) - {""}
    has_cases = cases.is_file() and bool(ids & set(pd.read_csv(cases, dtype=str, keep_default_na=False)["study_id"].str.strip().str.lower()))
    extra = "" if has_cases else " --no-cases"
    print(f"  sbatch --array=0-{(len(manifest) - 1) // args.chunk} scripts/clean_slurm.sh {args.manifest} {out_root} {args.chunk}{extra}")

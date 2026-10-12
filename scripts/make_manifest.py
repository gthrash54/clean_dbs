import sys
import argparse
from pathlib import Path
import pandas as pd
from clean_dbs import find_all_blocks, describe_block, missing_files
from clean_dbs.block_info import find_duplicates, recording_condition, source_path, output_path, recording_size

#a cluster job array needs a fixed numbered list of blocks, so task 7 always means the same block even if new data shows up mid-run
parser = argparse.ArgumentParser(description="Write one row per recording (path, patient, stage, block, format, cohort, session, condition, missing, duplicate) for clean.py --manifest")
parser.add_argument("raw_folder", help="whole inventory, one patient, or one surgery folder")
parser.add_argument("manifest", help="csv to write, keep it outside the repo since it has patient folder names")
parser.add_argument("--skip-missing", action="store_true", help="leave out blocks missing their signal files (.tev / .eeg)")
parser.add_argument("--chunk", type=int, default=1, help="blocks per array task, only used to print the --array range")
args = parser.parse_args()

blocks = find_all_blocks(args.raw_folder) #already sorted, so the row numbers are the same every time on the same data
#copies of one recording under several number folders all land on one output path; the extra copies stay in the manifest
#(so they show up in the report) but clean.py skips them
skipped = find_duplicates(blocks)
rows = []
for block in blocks:
    try:
        row = describe_block(block)
    except ValueError: #bad stage folder name, keep it so clean.py reports it as failed instead of it silently vanishing
        row = {"patient":"", "stage":"", "block":Path(block).name, "format":"bv" if Path(block).suffix.lower() == ".vhdr" else "tdt"}
    row["condition"], row["state"] = recording_condition(block)
    kept, why = skipped.get(block, ("", ""))
    row["duplicate"] = why #"duplicate" = same recording as duplicate_of, "conflict" = same name and place but a different size
    row["duplicate_of"] = source_path(kept) if kept else "" #from the study ID down, the full path is in that row's path
    row["missing"] = ";".join(missing_files(block)) #";" not "," so it stays one csv cell
    row["path"] = str(Path(block).resolve()) #full path so the job finds it no matter what folder it starts in
    row["signal_gb"] = round(recording_size(block) / 1e9, 3) #the signal files (.tev/.sev/.tsq or the .eeg)
    #cleaning reads one window at a time, ~0.5 GB peak whatever the size (measured); .sev blocks are read whole (raw + volts + a copy)
    row["est_peak_gb"] = round(0.5 + (3 * row["signal_gb"] if Path(block).is_dir() and any(Path(block).glob("*.sev")) else 0), 1)
    try:
        row["out"] = str(output_path(block, "")).lower() #only for the check below, not saved
    except ValueError:
        row["out"] = ""
    rows.append(row)

columns = ["path", "patient", "stage", "block", "format", "cohort", "session", "parent_id", "condition", "state", "missing", "duplicate", "duplicate_of",
           "signal_gb", "est_peak_gb"]
#two recordings that aren't marked as copies of each other must never write the same .h5 (lowercase: one file on Mac/Windows)
outs = pd.Series([r["out"] for r in rows if r["out"] and not r["duplicate"]], dtype=str)
clashes = sorted(set(outs[outs.duplicated()]))
if clashes:
    sys.exit(f"{len(clashes)} output path(s) would be written by more than one recording, fix the folder names first: {clashes[:10]}")
manifest = pd.DataFrame(rows, columns=columns).fillna("") #uh/UG rows have no cohort/session/parent_id, blank not NaN
n_missing = int((manifest["missing"] != "").sum())
print(f"{len(manifest)} recordings found, {n_missing} missing signal files, {len(skipped)} extra copies that clean.py will skip")
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

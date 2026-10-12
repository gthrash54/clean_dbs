import argparse
from pathlib import Path
import pandas as pd
from clean_dbs import find_all_blocks, describe_block, missing_files

#a cluster job array needs a fixed numbered list of blocks, so task 7 always means the same block even if new data shows up mid-run
parser = argparse.ArgumentParser(description="Write one row per recording (path, patient, stage, block, format, missing) for clean.py --manifest")
parser.add_argument("raw_folder", help="whole inventory, one patient, or one surgery folder")
parser.add_argument("manifest", help="csv to write, keep it outside the repo since it has patient folder names")
parser.add_argument("--skip-missing", action="store_true", help="leave out blocks missing their signal files (.tev / .eeg)")
parser.add_argument("--chunk", type=int, default=1, help="blocks per array task, only used to print the --array range")
args = parser.parse_args()

rows = []
for block in find_all_blocks(args.raw_folder): #already sorted, so the row numbers are the same every time on the same data
    try:
        row = describe_block(block)
    except ValueError: #bad stage folder name, keep it so clean.py reports it as failed instead of it silently vanishing
        row = {"patient":"", "stage":"", "block":Path(block).name, "format":"bv" if Path(block).suffix.lower() == ".vhdr" else "tdt"}
    row["missing"] = ";".join(missing_files(block)) #";" not "," so it stays one csv cell
    row["path"] = str(Path(block).resolve()) #full path so the job finds it no matter what folder it starts in
    rows.append(row)

manifest = pd.DataFrame(rows, columns=["path", "patient", "stage", "block", "format", "missing"])
n_missing = int((manifest["missing"] != "").sum())
print(f"{len(manifest)} recordings found, {n_missing} missing signal files")
if args.skip_missing:
    manifest = manifest[manifest["missing"] == ""].reset_index(drop=True) #reset so row numbers stay 0, 1, 2... with no gaps
    print(f"skipped {n_missing}, {len(manifest)} left in the manifest")

Path(args.manifest).parent.mkdir(parents=True, exist_ok=True)
manifest.to_csv(args.manifest, index=False)
print(f"saved {args.manifest}")
#slurm array numbers are 0 to (number of tasks - 1), each task does --chunk rows
if len(manifest) > 0:
    print(f"with --chunk {args.chunk} submit with --array=0-{(len(manifest) - 1) // args.chunk}")

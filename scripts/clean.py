import os
import sys
import argparse
from pathlib import Path
import pandas as pd
from clean_dbs import find_all_blocks, clean_block, describe_block, output_path

#three ways to run it:
#  python scripts/clean.py <raw folder> <out folder>                          cleans everything, one after another
#  python scripts/clean.py --manifest m.csv --task 3 --chunk 10 <out folder>  cleans rows 30-39 of the manifest (one cluster array task)
#  python scripts/clean.py --merge-reports <out folder>                        combines the per-task reports into clean_report.csv
parser = argparse.ArgumentParser(description="Clean DBS recordings into <out>/<patient>/<stage>/<block>.h5")
parser.add_argument("folders", nargs="*", help="<raw folder> <out folder>, or just <out folder> with --manifest")
parser.add_argument("--manifest", help="csv from scripts/make_manifest.py")
parser.add_argument("--task", type=int, help="which array task this is, defaults to $SLURM_ARRAY_TASK_ID")
parser.add_argument("--chunk", type=int, default=1, help="manifest rows per task")
parser.add_argument("--merge-reports", metavar="OUT_FOLDER", help="combine <out>/_reports/*.csv into <out>/clean_report.csv")
args = parser.parse_args()

if args.merge_reports:
    out_folder = Path(args.merge_reports)
    parts = sorted((out_folder / "_reports").glob("clean_report_*.csv"))
    if not parts:
        sys.exit(f"no reports in {out_folder / '_reports'}")
    report = pd.concat([pd.read_csv(p, keep_default_na=False) for p in parts], ignore_index=True)
    if args.manifest: #a task killed before its first block (time limit, out of memory) leaves no rows, so list those as "not run"
        manifest = pd.read_csv(args.manifest, keep_default_na=False)
        not_run = manifest.drop(index=report["row"], errors="ignore").drop(columns=["path", "missing"]) #path left out so the report has no full paths
        not_run["row"], not_run["status"], not_run["reason"] = not_run.index, "not run", "no report, check the slurm log for this task"
        report = pd.concat([report, not_run], ignore_index=True)
    report = report.sort_values("row").reset_index(drop=True)
    print(report["status"].value_counts().to_string())
    report.to_csv(out_folder / "clean_report.csv", index=False)
    print(f"merged {len(parts)} task reports into {out_folder / 'clean_report.csv'}")
    sys.exit()

if args.manifest:
    if len(args.folders) != 1:
        parser.error("with --manifest give just the out folder")
    out_folder = args.folders[0]
    task = args.task if args.task is not None else os.environ.get("SLURM_ARRAY_TASK_ID") #slurm sets this to the array number in each task
    if task is None:
        parser.error("give --task or run it inside a slurm array job")
    task = int(task)
    manifest = pd.read_csv(args.manifest, keep_default_na=False) #keep_default_na so an empty "missing" cell stays "" not NaN
    rows_to_do = list(range(task * args.chunk, min((task + 1) * args.chunk, len(manifest))))
    blocks = [manifest["path"][i] for i in rows_to_do]
    #each task writes its own file so parallel tasks never write the same csv at the same time
    report_path = Path(out_folder) / "_reports" / f"clean_report_{task}.csv"
    if not rows_to_do: #an --array range bigger than the manifest, nothing to do and nothing to report
        sys.exit(f"task {task} is past the end of the manifest ({len(manifest)} rows, chunk {args.chunk})")
    print(f"task {task}: manifest rows {rows_to_do}")
else:
    if len(args.folders) != 2:
        parser.error("give <raw folder> <out folder>")
    raw_folder, out_folder = args.folders
    rows_to_do = []
    blocks = find_all_blocks(raw_folder)
    report_path = Path(out_folder) / "clean_report.csv"

rows = []
for n, block in enumerate(blocks):
    try: #one bad block (missing file, corrupt data) shouldnt stop the rest
        already_done = output_path(block, out_folder).exists() #checked before cleaning so we can tell "skipped" from "cleaned"
        print(clean_block(block,out_folder))
        row = describe_block(block)
        if already_done: #skipped ones were cleaned on an earlier run, so adding new patients only cleans whats new
            row["status"] = "skipped"
        else:
            row["status"] = "cleaned"
        row["reason"] = ""
        rows.append(row)
    except Exception as error:
        print(f"FAILED {block}: {error}")
        try:
            row = describe_block(block)
        except ValueError: #if the stage folder name was the problem, describe_block fails too, so keep just the name
            row = {"block":Path(block).name}
        row["status"] = "failed"
        row["reason"] = str(error)
        rows.append(row)
    if rows_to_do: #manifest row number lets --merge-reports sort the table and spot rows that never ran
        row["row"] = rows_to_do[n]
        #saved after every block so a task killed partway still reports the blocks it finished
        report_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(report_path, index=False)

report = pd.DataFrame(rows)
print(report)
#saved next to the cleaned data so it outlives the script; pd.read_csv brings it back as a table
report_path.parent.mkdir(parents=True, exist_ok=True)
report.to_csv(report_path, index=False)
#nonzero exit makes sacct show FAILED for tasks with a failed block, so theyre easy to find
if args.manifest and not report.empty and (report["status"] == "failed").any():
    sys.exit(1)

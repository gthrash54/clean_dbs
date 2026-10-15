import os
import sys
import argparse
import hashlib
from datetime import datetime
from pathlib import Path
import shutil
import pandas as pd
from clean_dbs import find_all_blocks, clean_block, describe_block, output_path, privacy
from clean_dbs.block_info import find_duplicates, source_path, names_to_check, has_study_id
from clean_dbs.save_h5 import stale_reason, code_version

#three ways to run it:
#  python scripts/clean.py <raw folder> [<raw folder> ...] <out folder>      cleans everything, one after another
#  python scripts/clean.py --manifest m.csv --task 3 --chunk 10 <out folder>  cleans rows 30-39 of the manifest (one cluster array task)
#  python scripts/clean.py --merge-reports <out folder> --manifest m.csv       combines the per-task reports into clean_report.csv
parser = argparse.ArgumentParser(description="Clean DBS recordings into <out>/<patient>/<stage>/<block>.h5 (uh/UG) or <out>/<patient>/<session>/<block>.h5 (other cohorts)")
parser.add_argument("folders", nargs="*", help="<raw folder> [<raw folder> ...] <out folder>, or just <out folder> with --manifest")
parser.add_argument("--manifest", help="csv from scripts/make_manifest.py")
parser.add_argument("--task", type=int, help="which array task this is, defaults to $SLURM_ARRAY_TASK_ID")
parser.add_argument("--chunk", type=int, default=1, help="manifest rows per task")
parser.add_argument("--merge-reports", metavar="OUT_FOLDER", help="combine <out>/_reports/*.csv into <out>/clean_report.csv (needs --manifest)")
parser.add_argument("--force", action="store_true", help="clean every block again, even ones already cleaned by this version")
parser.add_argument("--no-cases", action="store_true", help="run without <out>/_cases/cases.csv (cohorts with no case tables); targets stay blank")
parser.add_argument("--no-dnu", action="store_true", help="run without <out>/_cases/dnu_channels.csv (the lab's do-not-use list); every channel_dnu is 0 and dnu_list says none")
parser.add_argument("--allow-dirty", action="store_true", help="let --manifest run on code with uncommitted changes (the files can't be traced to a commit)")
args = parser.parse_args()


def copy_guide(out_folder):
    #the guide travels with the data, so whoever opens the cleaned folder finds how it was made without asking
    guide = Path(__file__).resolve().parent.parent / "docs" / "user_guide.pdf"
    if guide.is_file():
        (Path(out_folder) / "_about").mkdir(parents=True, exist_ok=True)
        shutil.copy2(guide, Path(out_folder) / "_about" / "clean_dbs_user_guide.pdf")


def manifest_tag(manifest_path):
    #first 8 characters of the manifest's hash: a task's report is named after the manifest it ran, so reports from another one never mix in
    return hashlib.sha1(Path(manifest_path).read_bytes()).hexdigest()[:8]


def short_reason(text, block):
    #error text can hold the full local path (with the Box number folder above the study ID); keep paths from the study ID down, like source,
    #and with the names from renames.csv swapped in
    text, block = privacy.scrub(text), Path(block)
    depth = len(source_path(block).split("/"))
    for b in {block, block.resolve()}:
        if len(b.parts) > depth:
            above = Path(*b.parts[:-depth])
            for prefix in {str(above) + os.sep, above.as_posix() + "/"}:
                text = text.replace(prefix, "")
    return text


def report_name(manifest_path):
    #one merged report per manifest, so cohorts can share an out folder: manifest.csv -> clean_report.csv, manifest_dlpfc.csv -> clean_report_dlpfc.csv
    stem = Path(manifest_path).stem
    return "clean_report.csv" if stem == "manifest" else f"clean_report_{stem.removeprefix('manifest_')}.csv"


def held_back(block, excluded=""):
    #why a block isn't cleaned at all ("" = clean it): excluded in the manifest, no study ID, or a name that isn't cleared (privacy.py)
    if excluded:
        return excluded
    if not has_study_id(block):
        return f"no study ID above it or in its name ({privacy.mask(source_path(block))}), file it under a study ID first"
    names = names_to_check(block)
    return f"names not cleared ({len(names)}), see _cases/{privacy.TO_CHECK}" if names else ""


if args.merge_reports:
    if not args.manifest:
        parser.error("--merge-reports needs --manifest, the manifest the tasks ran")
    out_folder = Path(args.merge_reports)
    manifest = pd.read_csv(args.manifest, keep_default_na=False)
    merged = out_folder / report_name(args.manifest)
    tag = manifest_tag(args.manifest)
    every = sorted((out_folder / "_reports").glob("clean_report_*.csv"))
    parts = sorted((p for p in every if p.stem.endswith(f"_{tag}")), key=lambda p: p.stat().st_mtime) #oldest first
    if len(parts) < len(every):
        print(f"left out {len(every) - len(parts)} report(s) from another manifest (or made before reports were named after the manifest)")
    if not parts:
        sys.exit(f"no reports for this manifest ({tag}) in {out_folder / '_reports'}")
    report = pd.concat([pd.read_csv(p, keep_default_na=False) for p in parts], ignore_index=True)
    outside = ~report["row"].isin(manifest.index)
    if outside.any():
        print(f"warning: {int(outside.sum())} report row(s) aren't in the manifest (rows {sorted(set(report['row'][outside]))[:10]}), left out")
        report = report[~outside]
    #a rerun task (or the same rows run with another chunk size) reports its rows again; keep the newest try of each row
    n = len(report)
    report = report.sort_values("finished", kind="stable").drop_duplicates("row", keep="last")
    if len(report) < n:
        print(f"dropped {n - len(report)} older report row(s) for blocks that were run again")
    #a task killed before its first block (time limit, out of memory) leaves no rows, so list those as "not run"
    not_run = manifest.drop(index=report["row"], errors="ignore").drop(columns=["path", "missing", "out"], errors="ignore") #no path, so the report has no full paths
    not_run["row"], not_run["status"], not_run["reason"] = not_run.index, "not run", "no report, check the slurm log for this task"
    report = pd.concat([report, not_run], ignore_index=True).sort_values("row").reset_index(drop=True)
    print(report["status"].value_counts().to_string())
    report.to_csv(merged, index=False)
    copy_guide(out_folder)
    print(f"merged {len(parts)} task reports into {merged}")
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
    if not rows_to_do: #an --array range bigger than the manifest: nothing to do, and not an error
        print(f"task {task} is past the end of the manifest ({len(manifest)} rows, chunk {args.chunk}), nothing to do")
        sys.exit(0)
    commit, dirty = code_version()
    if (commit == "unknown" or dirty) and not args.allow_dirty: #every file records the commit, which only means something if the code is committed
        sys.exit(f"clean_dbs is at commit {commit} with uncommitted changes={dirty}; commit (or git stash) first so the cleaned files can be traced "
                 "to the code, or add --allow-dirty for a test run")
    blocks = [manifest["path"][i] for i in rows_to_do]
    #make_manifest.py already worked out the extra copies; older manifests have no duplicate column, so nothing is skipped
    skipped = {manifest["path"][i]:(manifest["duplicate_of"][i], manifest["duplicate"][i]) for i in rows_to_do
               if "duplicate" in manifest and manifest["duplicate"][i] != ""}
    excluded = {manifest["path"][i]:manifest["excluded"][i] for i in rows_to_do if "excluded" in manifest and manifest["excluded"][i] != ""}
    planned = {manifest["path"][i]:manifest["out"][i] for i in rows_to_do if "out" in manifest and manifest["out"][i] != ""}
    #each task writes its own file so parallel tasks never write the same csv at the same time; named after the chunk and manifest
    #so --merge-reports never mixes in reports from another manifest
    report_path = Path(out_folder) / "_reports" / f"clean_report_{task}_chunk{args.chunk}_{manifest_tag(args.manifest)}.csv"
    print(f"task {task}: manifest rows {rows_to_do}")
    try:
        privacy.use_out_folder(out_folder) #the same renames make_manifest.py used, so paths match the manifest's out column
    except ValueError as error:
        sys.exit(str(error))
else:
    if len(args.folders) < 2:
        parser.error("give <raw folder> [<raw folder> ...] <out folder>")
    *raw_folders, out_folder = args.folders
    rows_to_do = []
    try:
        privacy.use_out_folder(out_folder)
    except ValueError as error:
        sys.exit(str(error))
    blocks = sorted({Path(b).resolve(): b for raw in raw_folders for b in find_all_blocks(raw)}.values())
    excluded, planned = {block:held_back(block) for block in blocks}, {}
    excluded = {block:why for block, why in excluded.items() if why}
    #the same recording copied under several number folders would land on one output path, so only one copy is cleaned
    skipped = {block:(source_path(kept), why) for block, (kept, why) in find_duplicates([b for b in blocks if b not in excluded]).items()}
    report_path = Path(out_folder) / "clean_report.csv"

#the targets and case_status in every file come from the case tables, so a missing table would quietly leave them all blank
if not (Path(out_folder) / "_cases" / "cases.csv").is_file() and not args.no_cases:
    sys.exit(f"no _cases/cases.csv in the out folder; copy _cases (cases.csv, blocks.csv) in first, or add --no-cases for a cohort without case tables")
if (Path(out_folder) / "_cases" / "cases.csv").is_file() and not (Path(out_folder) / "_cases" / "dnu_channels.csv").is_file() and not args.no_dnu:
    #without it every file would quietly say "not listed", the same as a recording the lab never marked
    sys.exit("no _cases/dnu_channels.csv in the out folder; copy it in with the other case tables (scripts/make_dnu_table.py builds it), "
             "or add --no-dnu to clean without the lab's do-not-use list")

rows = []
for n, block in enumerate(blocks):
    held = held_back(block, excluded.get(block, "")) #checked again here, also for manifests made before the name check
    if held: #never cleaned, and no name in the report that isn't cleared
        print(f"EXCLUDED {privacy.mask(source_path(block))}: {held}")
        row = {"block":privacy.WITHHELD if "names" in held else Path(block).stem, "status":"excluded", "reason":held}
    elif block in skipped: #listed in the report as duplicate or conflict, with the copy that was cleaned instead
        row = describe_block(block)
        row["status"], row["reason"] = skipped[block][1], f"same output as {skipped[block][0]}"
    else:
        try: #one bad block (missing file, corrupt data) shouldnt stop the rest
            out = output_path(block, out_folder)
            if block in planned and Path(out).relative_to(out_folder).as_posix() != planned[block]:
                raise ValueError("its output path isn't the one in the manifest (renames.csv changed since?), make the manifest again")
            existed = out.exists() #checked before cleaning so we can tell cleaned / recleaned / skipped
            why = ("--force" if args.force else stale_reason(out)) if existed else ""
            print(clean_block(block, out_folder, force=args.force))
            row = describe_block(block)
            #skipped ones were cleaned by this version on an earlier run, so adding new patients only cleans whats new;
            #recleaned ones were there but broken, older than CLEAN_VERSION, or --force
            row["status"] = "recleaned" if existed and why else "skipped" if existed else "cleaned"
            row["reason"] = why
        except Exception as error:
            print(f"FAILED {privacy.mask(source_path(block))}: {short_reason(str(error), block)}")
            try:
                row = describe_block(block)
            except ValueError: #if the stage folder name was the problem, describe_block fails too, so keep just the name
                row = {"block":Path(block).name}
            row["status"] = "failed"
            row["reason"] = short_reason(str(error), block)
    if rows_to_do: #manifest row number lets --merge-reports sort the table and spot rows that never ran
        row["row"] = rows_to_do[n]
    row["finished"] = datetime.now().isoformat(timespec="microseconds") #--merge-reports keeps the newest try of a row
    rows.append(row)
    #saved after every block so a run killed partway still reports the blocks it finished
    report_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(report_path, index=False)

report = pd.DataFrame(rows)
print(report)
#saved next to the cleaned data so it outlives the script; pd.read_csv brings it back as a table
report_path.parent.mkdir(parents=True, exist_ok=True)
report.to_csv(report_path, index=False)
if not args.manifest: #array tasks skip this, --merge-reports does it once at the end
    copy_guide(out_folder)
#nonzero exit makes sacct show FAILED for tasks with a failed block (and a plain run's exit code says so too)
if not report.empty and (report["status"] == "failed").any():
    sys.exit(1)

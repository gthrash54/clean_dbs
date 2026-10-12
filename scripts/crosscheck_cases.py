import argparse
import shutil
from datetime import datetime
from pathlib import Path
import pandas as pd
from clean_dbs.crosscheck import crosscheck, read_crosswalk

#checks every row of the OR sheet tables (<cleaned>/_cases/cases.csv and blocks.csv) against the recordings and adds check columns,
#so only the conflicts need checking by hand against the paper sheets
parser = argparse.ArgumentParser(description="Cross-check cases.csv and blocks.csv against the recordings, adds check_* / status / check_detail columns")
parser.add_argument("cleaned_root", help="top of the cleaned set, the tables are in <cleaned_root>/_cases/")
parser.add_argument("data_root", help="raw data, one folder per study ID (e.g. ~/DBS Data)")
parser.add_argument("--box", help="Box raw data folder as mounted (e.g. ~/Box/UAB DBS Data/UAB DBS Raw Data); needed for stim checks and blocks that aren't local")
parser.add_argument("--crosswalk", help="crosswalk.csv (study_id, block, box_path), local only")
parser.add_argument("--cache", default=str(Path.home() / "clean_dbs_crosscheck_tmp"), help="folder for cached per-block results")
parser.add_argument("--no-stim", action="store_true", help="skip the Box reads, leads and block checks only")
parser.add_argument("--write", action="store_true", help="add the columns to the tables (a backup copy goes in _cases/backup_<time>/ first); without it, only prints the counts")
args = parser.parse_args()

folder = Path(args.cleaned_root).expanduser() / "_cases"
cases = pd.read_csv(folder / "cases.csv", dtype=str, keep_default_na=False)
blocks = pd.read_csv(folder / "blocks.csv", dtype=str, keep_default_na=False)
added = ["check_block", "check_leads_in_place", "check_stim_contacts", "check_target", "check_implant", "status", "check_detail"]
#rerunning replaces the old check columns instead of stacking new ones
cases_in, blocks_in = cases.drop(columns=[c for c in added if c in cases]), blocks.drop(columns=[c for c in added if c in blocks])
box = Path(args.box).expanduser() if args.box else None
cases_out, blocks_out = crosscheck(cases_in, blocks_in, Path(args.data_root).expanduser(), box, read_crosswalk(args.crosswalk), args.cache,
                                   read_stim=not args.no_stim, progress=True)

for col in ["check_block", "check_leads_in_place", "check_stim_contacts", "check_target", "status"]:
    print(f"blocks {col}: {blocks_out[col].value_counts().to_dict()}")
print(f"cases check_implant: {cases_out['check_implant'].value_counts().to_dict()}")
for _, r in blocks_out[blocks_out["status"] == "conflict"].iterrows():
    print(f"  {r['study_id']} {r['block']}: {r['check_detail']}")
for _, r in cases_out[cases_out["status"] == "conflict"].iterrows():
    print(f"  {r['study_id']} implant: {r['check_detail']}")

if args.write:
    backup = folder / f"backup_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    backup.mkdir()
    for name in ("cases.csv", "blocks.csv"):
        shutil.copy2(folder / name, backup / name)
    #the typed columns come back exactly as they were, the check columns go on the end
    cases_out.to_csv(folder / "cases.csv", index=False)
    blocks_out.to_csv(folder / "blocks.csv", index=False)
    print(f"wrote the check columns, old tables copied to {backup.name}/")

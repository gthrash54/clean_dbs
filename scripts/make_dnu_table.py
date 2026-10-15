import sys
import argparse
from pathlib import Path
import pandas as pd
from clean_dbs.dnu import dnu_rows

#builds _cases/dnu_channels.csv from the lab's do-not-use sheet ("UH3 or notes - data to not use.xlsx", sheet datatable_wide, on Box in
#UAB DBS Data/BRAIN Initiative grant beam steering STN GPi/UH3 data). Run it where the case tables live (the desktop), then copy the
#csv to Cheaha with the other case tables. Reading .xlsx needs openpyxl (pip install openpyxl); cleaning itself doesn't
parser = argparse.ArgumentParser(description="Build _cases/dnu_channels.csv from the lab's do-not-use sheet")
parser.add_argument("sheet", help="the lab's .xlsx, e.g. ~/DBS Cleaned/_cases/lab_sheets/UH3 or notes - data to not use.xlsx")
parser.add_argument("cases", help="the _cases folder to write dnu_channels.csv into (blocks.csv there is used to list rows that match no block)")
args = parser.parse_args()

sheet = pd.read_excel(args.sheet, sheet_name="datatable_wide", engine="openpyxl")
rows, unknown = dnu_rows(sheet, f"lab sheet: {Path(args.sheet).name}")
rows = rows.drop_duplicates() #the sheet lists a few recordings twice (uh3001 increment2_depth3)
if unknown:
    print(f"{len(unknown)} names follow no known pattern and were left out: {sorted(set(unknown))}")
cases = Path(args.cases)
rows.to_csv(cases / "dnu_channels.csv", index=False)
listed = rows[["study_id", "stage", "block"]].drop_duplicates()
print(f"{len(listed)} recordings, {int((rows.channel != '').sum())} channels marked do not use -> {cases / 'dnu_channels.csv'}")
if (cases / "blocks.csv").is_file(): #rows that won't match any recording (typos, notes typed in the block column) are listed to fix by hand
    blocks = pd.read_csv(cases / "blocks.csv", dtype=str, keep_default_na=False)
    known = set(zip(blocks.study_id.str.lower(), blocks.get("stage", pd.Series([""] * len(blocks))).str.lower(), blocks.block))
    orphans = [tuple(x) for x in listed.values if tuple(x) not in known and (x[0], "", x[2]) not in known]
    for study_id, stage, block in orphans:
        print(f"  no block in blocks.csv for {study_id} {stage} {block!r}")
    print(f"{len(orphans)} of {len(listed)} recordings match no row of blocks.csv")
sys.exit(0)

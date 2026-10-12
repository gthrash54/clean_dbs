import sys
from pathlib import Path
import pandas as pd
from clean_dbs import find_all_blocks, clean_block, describe_block, output_path

raw_folder = sys.argv[1]
out_folder = sys.argv[2]
rows = []
for block in find_all_blocks(raw_folder):
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

report = pd.DataFrame(rows)
print(report)
#saved next to the cleaned data so it outlives the script; pd.read_csv brings it back as a table
report.to_csv(Path(out_folder) / "clean_report.csv", index=False)

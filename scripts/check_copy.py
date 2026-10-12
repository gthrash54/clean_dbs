import sys
import argparse
from pathlib import Path
from clean_dbs import find_all_blocks, missing_files
from clean_dbs.block_info import source_path
from clean_dbs.read_tdt import tev_problem

#checks a copied raw folder before cleaning: every recording has its signal file, and no .tev was cut short by a stopped copy
#(tdt reads a short .tev without an error and fills the rest with zeros, so this is the only place it shows up)
parser = argparse.ArgumentParser(description="List recordings in a copied raw folder that are missing their signal file or have a cut-short .tev")
parser.add_argument("raw_folder", help="the copied folder, e.g. $USER_DATA/dbs_raw")
args = parser.parse_args()

blocks = find_all_blocks(args.raw_folder)
problems = 0
for block in blocks:
    try:
        missing = missing_files(block)
        problem = f"missing {';'.join(missing)}" if missing else ""
        if not problem and Path(block).is_dir() and any(Path(block).glob("*.tev")):
            problem = tev_problem(block)
    except Exception as error: #a .tsq that won't even parse is a problem too
        problem = f"could not check: {type(error).__name__}: {str(error).replace(str(Path(block)), source_path(block))}" #no full paths
    if problem:
        problems += 1
        print(f"{source_path(block)}: {problem}") #from the study ID down, never the full path
print(f"{len(blocks)} recordings checked, {problems} with problems")
sys.exit(1 if problems else 0)

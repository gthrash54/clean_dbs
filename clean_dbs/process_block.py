from pathlib import Path
import re
from clean_dbs.read_tdt import read_tdt
from clean_dbs.save_h5 import save_h5

def process_block(block_dir, out_root):
    """
    Inputs:
    block_dir is the block to process (ie one item from find_tdt_blocks)
    out_root is the top of the cleaned set where outputs go
  
    """
    block_dir = Path(block_dir)
    out_root = Path(out_root)
    patient = block_dir.parent.parent.name
    #folder names vary ("stage1 surgery", "stage 2 surgery", "stage2 surgery (sedated, not GA)"), so pull out just "stage" + the number
    match = re.search(r"stage\s*(\d)", block_dir.parent.name, re.IGNORECASE)
    if match is None:
        raise ValueError(f"Can't tell the stage from {block_dir.parent.name!r}, expected a folder like 'stage1 surgery' above {block_dir.name}")
    stage = "stage" + match.group(1)
    out_path = out_root / patient / stage / (block_dir.name + ".h5")
    if out_path.exists():
        return out_path
    r = read_tdt(block_dir)
    #cleaning goes here(next)
    save_h5(r, out_path, block_dir)
    return out_path
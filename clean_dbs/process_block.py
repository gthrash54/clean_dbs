from pathlib import Path
import re
from clean_dbs.read_tdt import read_tdt
from clean_dbs.read_bv import read_bv
from clean_dbs.save_h5 import save_h5

def process_block(block_dir, out_root):
    """
    Inputs:
    block_dir is the block to process: a tdt block folder (from find_tdt_blocks) or a .vhdr file (from find_bv_recordings)
    out_root is the top of the cleaned set where outputs go
    Returns: path of the saved file, out_root/<patient>/<stage>/<block>.h5
    """
    block_dir = Path(block_dir)
    out_root = Path(out_root)
    patient = block_dir.parent.parent.name
    #folder names vary ("stage1 surgery", "stage 2 surgery", "stage2 surgery (sedated, not GA)"), so pull out just "stage" + the number
    match = re.search(r"stage\s*(\d)", block_dir.parent.name, re.IGNORECASE)
    if match is None:
        raise ValueError(f"Can't tell the stage from {block_dir.parent.name!r}, expected a folder like 'stage1 surgery' above {block_dir.name}")
    stage = "stage" + match.group(1)
    is_bv = block_dir.suffix.lower() == ".vhdr"
    name = block_dir.stem if is_bv else block_dir.name #stem drops the .vhdr, tdt block folders keep their full name
    out_path = out_root / patient / stage / (name + ".h5")
    if out_path.exists():
        return out_path
    if is_bv: #brainvision is a .vhdr file, tdt is a block folder
        r = read_bv(block_dir)
    else:
        r = read_tdt(block_dir)
    #cleaning goes here(next)
    save_h5(r, out_path, block_dir)
    return out_path
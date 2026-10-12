from pathlib import Path
from clean_dbs.read_tdt import read_tdt
from clean_dbs.read_bv import read_bv
from clean_dbs.save_h5 import save_h5
from clean_dbs.block_info import output_path
from clean_dbs.flags import make_flags

def clean_block(block_dir, out_root):
    """
    Inputs:
    block_dir is the block to process: a tdt block folder (from find_tdt_blocks) or a .vhdr file (from find_bv_recordings)
    out_root is the top of the cleaned set where outputs go
    Returns: path of the saved file, out_root/<patient>/<stage>/<block>.h5
    """
    block_dir = Path(block_dir)
    out_root = Path(out_root)
    out_path = output_path(block_dir, out_root) #works out patient/stage/name, see block_info.py
    if out_path.exists():
        return out_path
    if block_dir.suffix.lower() == ".vhdr": #brainvision is a .vhdr file, tdt is a block folder
        r = read_bv(block_dir)
    else:
        r = read_tdt(block_dir)
    r["flags"] = make_flags(r, block_dir) #saturation, new segments, dnu, stored as notes, the data isn't changed
    save_h5(r, out_path, block_dir)
    return out_path
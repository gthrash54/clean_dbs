from pathlib import Path
from clean_dbs.read_tdt import read_tdt
from clean_dbs.read_bv import read_bv
from clean_dbs.save_h5 import save_h5
from clean_dbs.block_info import output_path
from clean_dbs.flags import make_flags
from clean_dbs.datatypes import stream_type, event_type

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
    for s in r["streams"].values(): #tag each stream lfp / mer / emg / audio / stim_monitor, see datatypes.py
        s["type"], s["subtype"] = stream_type(s["device"])
    is_bv = block_dir.suffix.lower() == ".vhdr"
    for name, e in r["epocs"].items():
        e["type"] = event_type(name, is_bv)
    for name, sc in r["scalars"].items():
        sc["type"] = event_type(name, is_bv)
    #the video stays outside the .h5 (too big, and video tools cant read it from inside); record where it is, Cam1 lines frames up with the data
    videos = sorted(block_dir.glob("*.avi")) if block_dir.is_dir() else []
    r["info"]["video"] = ";".join("/".join(v.parts[-4:]) for v in videos) #from the study ID down, like source
    r["flags"] = make_flags(r, block_dir) #saturation, new segments, dnu, stored as notes, the data isn't changed
    save_h5(r, out_path, block_dir)
    return out_path
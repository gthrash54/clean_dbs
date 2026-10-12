from pathlib import Path

def find_tdt_blocks(folder, include_dnu=False):
    """
    Input:
    Folder: can be the whole inventory, one patient, or one surgery folder
    include_dnu: false skips plock with dnu (do not use) in the name, set true to keep
    Returns: List of block folders
    """
    folder = Path(folder)

    blocks = []
    for file in folder.rglob("*.tsq"):
        if ("DNU" not in file.parent.name) or include_dnu:
            blocks.append(file.parent)

    return sorted(blocks)
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
        if ("dnu" not in file.parent.name.lower()) or include_dnu: #lower() so "_DNU" and "_dnu" both count
            blocks.append(file.parent)

    return sorted(blocks)
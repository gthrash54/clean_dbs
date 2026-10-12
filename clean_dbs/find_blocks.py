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
        if _in_imaging(file, folder):
            continue
        if ("dnu" not in file.parent.name.lower()) or include_dnu: #lower() so "_DNU" and "_dnu" both count
            blocks.append(file.parent)

    return sorted(blocks)

def find_bv_recordings(folder, include_dnu=False):
    """
    Input:
    Folder: can be the whole inventory, one patient, or one surgery folder
    include_dnu: false skips recordings with dnu (do not use) in the name, set true to keep
    Returns: List of .vhdr files (a brainvision recording is one file, not a folder like tdt)
    """
    folder = Path(folder)

    recordings = []
    for file in folder.rglob("*.vhdr"):
        if _in_imaging(file, folder):
            continue
        if ("dnu" not in file.name.lower()) or include_dnu: #checks the file name since theres no block folder
            recordings.append(file)

    return sorted(recordings)


def find_all_blocks(folder, include_dnu=False):
    """
    Input:
    Folder: can be the whole inventory, one patient, or one surgery folder
    include_dnu: false skips anything with dnu (do not use) in the name, set true to keep
    Returns: List of every recording, tdt block folders and brainvision .vhdr files together, sorted
    """
    return sorted(find_tdt_blocks(folder, include_dnu) + find_bv_recordings(folder, include_dnu))


def _in_imaging(file, folder):
    #scans aren't recordings, and some copies on cheaha still have Imaging folders from before the copy filter was fixed
    return any(part.lower() == "imaging" for part in file.relative_to(folder).parts[:-1])

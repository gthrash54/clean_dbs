import re
from pathlib import Path


def describe_block(block):
    """
    Input:
    block: a tdt block folder or a brainvision .vhdr file (one item from find_all_blocks)
    Returns: dictionary with patient, stage, block (name), format ("tdt" or "bv")
    """
    block = Path(block)
    is_bv = block.suffix.lower() == ".vhdr"
    #folder names vary ("stage1 surgery", "stage 2 surgery", "stage2 surgery (sedated, not GA)"), so pull out just "stage" + the number
    match = re.search(r"stage\s*(\d)", block.parent.name, re.IGNORECASE)
    if match is None:
        raise ValueError(f"Can't tell the stage from {block.parent.name!r}, expected a folder like 'stage1 surgery' above {block.name}")
    return {"patient":block.parent.parent.name, "stage":"stage" + match.group(1),
            "block":block.stem if is_bv else block.name, "format":"bv" if is_bv else "tdt"}


def missing_files(block):
    """
    Input:
    block: a tdt block folder or a brainvision .vhdr file
    Returns: list of the signal files that should be there but aren't (empty list = everything is there)
    """
    block = Path(block)
    if block.suffix.lower() == ".vhdr":
        #the .vhdr names its signal file on the "DataFile=" line
        data_file = next((line.split("=", 1)[1].strip() for line in block.read_text(errors="ignore").splitlines() if line.startswith("DataFile=")), block.stem + ".eeg")
        return [] if (block.parent / data_file).is_file() else [data_file]
    #tdt keeps the signal in .tev (or .sev files on some rigs)
    return [] if any(block.glob("*.tev")) or any(block.glob("*.sev")) else [".tev"]

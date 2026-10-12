import pandas as pd
from clean_dbs.find_blocks import find_all_blocks
from clean_dbs.block_info import describe_block, missing_files

def check_blocks(folder):
    """
    Input:
    folder: can be the whole inventory, one patient, or one surgery folder
    Returns: table (pandas DataFrame) with one row per recording: patient, stage, block, format, missing
        missing is a list of signal files that arent there, empty means its ready to clean
    Doesnt read any data so its fast even on the whole inventory
    """
    rows = []
    for block in find_all_blocks(folder):
        row = describe_block(block)
        row["missing"] = missing_files(block)
        rows.append(row)
    return pd.DataFrame(rows)
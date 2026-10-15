import re
from functools import lru_cache
from pathlib import Path
import pandas as pd

#which channels the lab marked "do not use" in a recording (bad or noisy contacts, contacts still in the cannula, bad ECoG or EMG),
#from <cleaned root>/_cases/dnu_channels.csv, which scripts/make_dnu_table.py builds from the lab's sheet
#("UH3 or notes - data to not use.xlsx"). Stored with every channel as channel_dnu; the data itself is never changed

COLUMNS = ["study_id", "stage", "block", "channel", "source"]


@lru_cache(maxsize=8)
def _load(path):
    if not Path(path).is_file():
        return None
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def dnu_channels(out_root, study_id, stage, block, names, channels=(), stream_type=""):
    """
    Which channels of one stream the lab marked do not use in this recording.
    Input:
    out_root: top of the cleaned set (holds _cases/dnu_channels.csv)
    study_id, stage, block: the recording (study ID matched without caring about case; a blank stage in the table matches any)
    names: the stream's channel names (lead1_c3, ecog_a4, ...), or None for a stream without names
    channels: the stream's input numbers. EMG streams without names are matched by them: the lab's sheet calls them muscle_<input>
    stream_type: the stream's type (lfp, emg, ...); only "emg" streams match muscle_<input>
    Returns: (flags, listed): flags is a list with one 0/1 per channel (1 = do not use); listed is the table's source for this
        recording ("lab sheet: ...") when it has any row, even one saying no channel is bad, and "" when the recording isn't
        listed at all (then nothing is known and every flag is 0)
    """
    n = len(names) if names is not None else len(channels)
    flags = [0] * n
    table = _load(str(Path(out_root) / "_cases" / "dnu_channels.csv"))
    if table is None:
        return flags, ""
    rows = table[(table["study_id"].str.strip().str.lower() == str(study_id).strip().lower()) & (table["block"] == block)]
    if stage:
        row_stage = rows["stage"].str.replace(" ", "").str.lower()
        rows = rows[(row_stage == str(stage).replace(" ", "").lower()) | (row_stage == "")]
    if rows.empty:
        return flags, ""
    bad = set(rows["channel"]) - {""}
    for i in range(n):
        name = str(names[i]) if names is not None else ""
        muscle = f"muscle_{channels[i]}" if stream_type == "emg" and i < len(channels) else None
        if name in bad or (muscle is not None and muscle in bad):
            flags[i] = 1
    return flags, rows["source"].iloc[0]


def dnu_rows(sheet, source):
    """
    Turns the lab's do-not-use sheet into dnu_channels.csv rows.
    Input:
    sheet: the sheet as a table (columns id, surgery, trial, contact_dnu1, contact_dnu2, ...), one row per recording
    source: what to write in the source column, e.g. "lab sheet: UH3 or notes - data to not use.xlsx"
    Returns: (rows, unknown): rows is a table with COLUMNS, one row per channel marked do not use, plus a row with a blank channel for
        a recording listed with none; the sheet's names become clean_dbs names (dbs1_8 -> lead1_c8, ecog_a4 and muscle_4 stay).
        unknown lists the names that don't follow any known pattern (left out)
    """
    out, unknown = [], []
    columns = [c for c in sheet.columns if str(c).startswith("contact_dnu")]
    for _, r in sheet.iterrows():
        study_id = re.sub(r"\s+", "", str(r["id"])).lower() #the sheet has "uh 3023"
        stage, block = str(r["surgery"]).strip().lower(), str(r["trial"]).strip()
        names = [str(r[c]).strip() for c in columns if pd.notna(r[c]) and str(r[c]).strip()]
        if not names:
            out.append([study_id, stage, block, "", source])
        for name in names:
            lead = re.fullmatch(r"dbs(\d)_(\d+)", name)
            if lead:
                out.append([study_id, stage, block, f"lead{lead.group(1)}_c{lead.group(2)}", source])
            elif re.fullmatch(r"ecog_[ap]\d+|muscle_\d+", name):
                out.append([study_id, stage, block, name, source])
            else:
                unknown.append(name)
    return pd.DataFrame(out, columns=COLUMNS), unknown

import h5py
import numpy as np
import pandas as pd
from clean_dbs.dnu import dnu_channels, dnu_rows, _load, COLUMNS


def _table(root, rows):
    #a tiny _cases/dnu_channels.csv
    (root / "_cases").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=COLUMNS).to_csv(root / "_cases" / "dnu_channels.csv", index=False)
    _load.cache_clear()


NAMES = ["ecog_a1", "ecog_a4", "lead1_c1", "lead1_c8", "lead2_c3"]


def test_flags_by_name_and_listed_or_not(tmp_path):
    _table(tmp_path, [["uh9999", "stage1", "increment1_depth1", "lead1_c8", "lab sheet: x.xlsx"],
                      ["uh9999", "stage1", "increment1_depth1", "ecog_a4", "lab sheet: x.xlsx"],
                      ["uh9999", "stage1", "rest", "", "lab sheet: x.xlsx"]]) #listed, nothing bad
    flags, listed = dnu_channels(tmp_path, "UH9999", "stage1", "increment1_depth1", NAMES) #study ID case doesn't matter
    assert flags == [0, 1, 0, 1, 0] and listed == "lab sheet: x.xlsx"
    assert dnu_channels(tmp_path, "uh9999", "stage1", "rest", NAMES) == ([0] * 5, "lab sheet: x.xlsx") #checked: none bad
    assert dnu_channels(tmp_path, "uh9999", "stage1", "crosstalk", NAMES) == ([0] * 5, "") #not listed: nothing known
    assert dnu_channels(tmp_path, "uh9999", "stage2", "increment1_depth1", NAMES) == ([0] * 5, "") #other stage


def test_emg_by_input_number(tmp_path):
    #the lab's sheet calls EMG inputs muscle_<n>; TDT EMG streams have input numbers, no names
    _table(tmp_path, [["uh9999", "stage1", "motor_baseline", "muscle_4", "lab sheet: x.xlsx"]])
    flags, _ = dnu_channels(tmp_path, "uh9999", "stage1", "motor_baseline", None, channels=[1, 2, 3, 4], stream_type="emg")
    assert flags == [0, 0, 0, 1]
    flags, _ = dnu_channels(tmp_path, "uh9999", "stage1", "motor_baseline", None, channels=[1, 2, 3, 4], stream_type="lfp")
    assert flags == [0, 0, 0, 0] #only EMG streams match muscle_<n>


def test_no_table(tmp_path):
    _load.cache_clear()
    assert dnu_channels(tmp_path, "uh9999", "stage1", "rest", NAMES) == ([0] * 5, "")


def test_rows_from_the_lab_sheet():
    sheet = pd.DataFrame({"id": ["uh 3023", "uh3001"], "surgery": ["stage1", "stage1"], "trial": ["increment1_depth1", "rest"],
                          "contact_dnu1": ["dbs1_8", np.nan], "contact_dnu2": ["muscle_4", np.nan], "contact_dnu3": ["ecog_p2", np.nan],
                          "contact_dnu4": ["typo!", np.nan]})
    rows, unknown = dnu_rows(sheet, "lab sheet: x.xlsx")
    assert rows.values.tolist() == [["uh3023", "stage1", "increment1_depth1", "lead1_c8", "lab sheet: x.xlsx"],
                                    ["uh3023", "stage1", "increment1_depth1", "muscle_4", "lab sheet: x.xlsx"],
                                    ["uh3023", "stage1", "increment1_depth1", "ecog_p2", "lab sheet: x.xlsx"],
                                    ["uh3001", "stage1", "rest", "", "lab sheet: x.xlsx"]] #listed with none bad
    assert unknown == ["typo!"]


def test_clean_block_writes_channel_dnu(tmp_path):
    from test_final_fixes import recording
    from clean_dbs.clean_block import clean_block
    vhdr = recording(tmp_path / "raw" / "ug9999" / "stage2 surgery")
    out = tmp_path / "out"
    first = clean_block(vhdr, out)
    with h5py.File(first) as f: #no table: every channel 0, and the file says the recording wasn't listed
        names = {s: [x.decode() for x in f["streams"][s]["channel_names"][:]] for s in f["streams"] if "channel_names" in f["streams"][s]}
        assert all(f["streams"][s]["channel_dnu"][:].tolist() == [0] * len(n) for s, n in names.items())
        assert all(f["streams"][s].attrs["dnu_list"] == "none" for s in f["streams"])
    stream, picked = next((s, n[0]) for s, n in names.items() if n)
    _table(out, [["ug9999", "stage2", vhdr.stem, picked, "lab sheet: x.xlsx"]])
    again = clean_block(vhdr, out, force=True)
    with h5py.File(again) as f:
        g = f["streams"][stream]
        assert g["channel_dnu"][:].tolist() == [1] + [0] * (len(names[stream]) - 1) and g["channel_dnu"].dtype == np.int8
        assert g.attrs["dnu_list"] == "lab sheet: x.xlsx"


def test_block_fixes_and_session_wide_rows(tmp_path):
    from clean_dbs.dnu import apply_block_fixes
    rows = pd.DataFrame([["uh9999", "stage1", "clincal_dbs1", "muscle_3", "s"], ["uh9999", "stage1", "gone", "muscle_4", "s"],
                         ["uh9999", "stage1", "a note, not a block", "muscle_4", "s"]], columns=COLUMNS)
    fixes = pd.DataFrame([["uh9999", "stage1", "clincal_dbs1", "clinical_dbs1"], ["uh9999", "stage1", "gone", ""],
                          ["uh9999", "stage1", "a note, not a block", "*"]], columns=["study_id", "stage", "sheet_block", "block"])
    fixed = apply_block_fixes(rows, fixes)
    assert fixed.values.tolist() == [["uh9999", "stage1", "clinical_dbs1", "muscle_3", "s"], ["uh9999", "stage1", "*", "muscle_4", "s"]]
    _table(tmp_path, fixed.values.tolist())
    names = ["muscle_x"] #an EMG stream, matched by input number
    assert dnu_channels(tmp_path, "uh9999", "stage1", "any_recording", None, [4], "emg") == ([1], "s") #"*" reaches every recording
    assert dnu_channels(tmp_path, "uh9999", "stage1", "clinical_dbs1", None, [3, 4], "emg") == ([1, 1], "s")
    assert dnu_channels(tmp_path, "uh9999", "stage2", "any_recording", None, [4], "emg") == ([0], "") #only that stage

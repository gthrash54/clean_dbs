from pathlib import Path
import pytest
from clean_dbs.cohorts import COHORTS, cohort_of, condition_of
from clean_dbs.block_info import describe_block, output_path, find_duplicates, source_path, recording_condition
from clean_dbs.find_blocks import find_all_blocks
from clean_dbs import privacy

#made-up Box-like paths: <number folder>/<study id>/...; numbers and IDs are fake
BOX = Path("box") / "100.0000001"


def test_cohort_of():
    assert cohort_of("uh3021") == "uh" and cohort_of("UG3005") == "ug"
    assert cohort_of("TMSCS001") == "tmscs" and cohort_of("TM001") == "tm" #longest prefix first
    assert cohort_of("ORBV001") == "orbv" and cohort_of("BV001") == "bv" and cohort_of("BPP001") == "bpp"
    assert cohort_of("dlpfc045") == "dlpfc" and cohort_of("SEP_03") == "sep" and cohort_of("Burst012") == "burst"
    #sessions and other folders that only look a bit like IDs
    for name in ("PP1", "ORBV1_Anesth", "AXCPT-009", "DIR009", "LFP", "stage1 surgery", "6 month", "X123"):
        assert cohort_of(name) is None, name
    assert all(p == p.lower() for p in COHORTS)


def test_uh_ug_unchanged():
    #exactly the four keys and the same output path as before cohorts existed, even with the number folder above
    info = describe_block(BOX / "uh9999" / "stage 2 surgery" / "STN_dbs1.vhdr")
    assert info == {"patient":"uh9999", "stage":"stage2", "block":"STN_dbs1", "format":"bv"}
    assert output_path(BOX / "uh9999" / "stage1 surgery" / "clinical_dbs1", "out") == Path("out") / "uh9999" / "stage1" / "clinical_dbs1.h5"
    assert output_path(BOX / "UG9999" / "stage2 surgery (sedated, not GA)" / "LFP1.vhdr", "out") == Path("out") / "UG9999" / "stage2" / "LFP1.h5"
    assert source_path(BOX / "uh9999" / "stage1 surgery" / "clinical_dbs1") == "uh9999/stage1 surgery/clinical_dbs1"
    with pytest.raises(ValueError, match="stage"):
        describe_block(BOX / "UG9999" / "LFP1.vhdr") #UG still needs its stage folder


def test_files_at_id_root():
    info = describe_block(BOX / "DLPFC045" / "rest.vhdr")
    assert info["patient"] == "DLPFC045" and info["session"] == "" and info["stage"] == "" and info["cohort"] == "dlpfc"
    assert info["condition"] == "rest" and info["parent_id"] == ""
    assert output_path(BOX / "DLPFC045" / "rest.vhdr", "out") == Path("out") / "DLPFC045" / "rest.h5"
    assert describe_block(BOX / "dlpfc045" / "rest.vhdr")["patient"] == "DLPFC045" #one spelling per person


def test_session_subfolders():
    info = describe_block(BOX / "DLPFC045" / "Task data" / "AXCPT" / "ieeg" / "task1.vhdr")
    assert info["session"] == "task_data/axcpt/ieeg" and info["condition"] == "task_axcpt"
    assert output_path(BOX / "DLPFC045" / "Task data" / "AXCPT" / "ieeg" / "task1.vhdr", "out") == Path("out") / "DLPFC045" / "task_data" / "axcpt" / "ieeg" / "task1.h5"
    #the same file name in two subfolders lands in two places
    a = output_path(BOX / "BI012" / "LFP" / "LFP1.vhdr", "out")
    b = output_path(BOX / "BI012" / "6 month" / "LFP" / "LFP1.vhdr", "out")
    assert a == Path("out") / "BI012" / "lfp" / "LFP1.h5" and b == Path("out") / "BI012" / "6_month" / "lfp" / "LFP1.h5"
    assert describe_block(BOX / "BI012" / "6 month" / "LFP" / "LFP1.vhdr")["state"] == "followup_6mo"


def test_awake_asleep_and_anesth():
    assert describe_block(BOX / "MC010" / "Asleep" / "1-2.vhdr")["state"] == "asleep"
    assert describe_block(BOX / "MC010" / "1-2_asleep.vhdr")["condition"] == "contact_pair"
    info = describe_block(BOX / "ORBV002" / "ORBV1_Anesth" / "1+2-.vhdr")
    assert info["session"] == "orbv1_anesth" and info["state"] == "asleep" and info["condition"] == "stim_bipolar_config"


def test_nested_ids():
    #an ID inside another ID belongs to the inner one, and the outer one is kept
    info = describe_block(BOX / "PP012" / "CHR003" / "1-2-3.vhdr")
    assert info["patient"] == "CHR003" and info["parent_id"] == "PP012" and info["cohort"] == "chr" and info["session"] == ""
    info = describe_block(BOX / "LFP03" / "CHR003" / "2_mono" / "x.vhdr")
    assert info["patient"] == "CHR003" and info["parent_id"] == "LFP03" and info["session"] == "2_mono"
    #PP1 is a session (one digit), not an ID
    info = describe_block(BOX / "PP012" / "PP1" / "pp1.vhdr")
    assert info["patient"] == "PP012" and info["session"] == "pp1" and info["condition"] == "stim_paired_pulse"
    assert source_path(BOX / "PP012" / "CHR003" / "1-2-3.vhdr") == "PP012/CHR003/1-2-3.vhdr"


def test_stage_folder_in_other_cohort():
    info = describe_block(BOX / "BB004" / "Stage 2 surgery" / "LFP1.vhdr")
    assert info["stage"] == "stage2" and info["session"] == "stage_2_surgery"


def test_straight_in_number_folder():
    info = describe_block(BOX / "BPP012.vhdr")
    assert info["patient"] == "BPP012" and info["cohort"] == "bpp" and info["session"] == ""
    assert output_path(BOX / "BPP012.vhdr", "out") == Path("out") / "BPP012" / "BPP012.h5"
    with pytest.raises(ValueError, match="study ID"):
        describe_block(BOX / "baseline.vhdr")
    with pytest.raises(ValueError, match="study ID"):
        describe_block(BOX / "MOTIF" / "1ma1.vhdr") #a cohort folder without an ID number


def test_tdt_in_other_cohort():
    info = describe_block(BOX / "OR001" / "day1" / "block1")
    assert info["format"] == "tdt" and info["block"] == "block1" and info["session"] == "day1"


def test_condition_of():
    assert condition_of("rest eyes closed") == ("rest_eyes_closed", "")
    assert condition_of("exp1_postdbs")[1] == "post_dbs"
    assert condition_of("itbs") == ("stim_itbs", "")
    assert condition_of("burstpp")[0] == "stim_burst"
    assert condition_of("1_c_130hz_3ma_rev")[0] == "stim_setting"
    assert condition_of("task1", "joy")[0] == "task_joystick"
    assert condition_of("lfp_lesion")[0] == "lesion"
    assert condition_of("something") == ("", "")
    assert recording_condition(BOX / "uh9999" / "stage 2 surgery" / "contact1.vhdr") == ("contact", "")


def test_postop_state():
    #"post op", "postop" and "post_op" folders are on Box (DLPFC) and used to give no state; postop is kept apart from post_dbs
    for folder in ("post op", "Postop", "post_op", "post-op"):
        assert condition_of(f"{folder}/rest") == ("rest", "postop")
        assert privacy.is_clear(folder) #post, op and postop are already known words, so the names check lets them through
    assert condition_of("post op/exp1_postdbs")[1] == "post_dbs;postop" #both apply, neither eats the other
    assert condition_of("exp1_postdbs")[1] == "post_dbs" and condition_of("post dbs")[1] == "post_dbs"
    assert condition_of("impedencetest_postopening")[1] == "" #like preopening, opening isn't op
    #the other two-word states take a space or - too
    assert condition_of("pre op/x")[1] == "preop" and condition_of("x/off meds")[1] == "off_meds" and condition_of("x/on meds")[1] == "on_meds"
    assert condition_of("x/stimulation meds")[1] == "" #"on" at the end of a word isn't on_meds


def write_bv(folder, name, n_bytes):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / f"{name}.vhdr").write_text(f"DataFile={name}.eeg\n")
    (folder / f"{name}.eeg").write_bytes(b"x" * n_bytes)


def test_find_duplicates(tmp_path):
    a, b, c = tmp_path / "100.1", tmp_path / "100.2", tmp_path / "100.3"
    write_bv(a / "DLPFC045", "rest", 10)
    write_bv(b / "dlpfc045", "rest", 10) #same recording, another number folder, other spelling
    write_bv(c / "DLPFC045", "rest", 4) #same name and place, smaller: an incomplete copy, or a different recording
    write_bv(a / "DLPFC045", "itbs", 10)
    write_bv(a / "PP012" / "CHR003", "1-2", 7) #nested copy of the CHR003 folder
    write_bv(b / "CHR003", "1-2", 7)
    skipped = find_duplicates(find_all_blocks(tmp_path))
    assert skipped == {b / "dlpfc045" / "rest.vhdr":(a / "DLPFC045" / "rest.vhdr", "duplicate"),
                       c / "DLPFC045" / "rest.vhdr":(a / "DLPFC045" / "rest.vhdr", "conflict"),
                       b / "CHR003" / "1-2.vhdr":(a / "PP012" / "CHR003" / "1-2.vhdr", "duplicate")}
    #what's left has one recording per output path
    kept = [r for r in find_all_blocks(tmp_path) if r not in skipped]
    assert len({output_path(r, "out") for r in kept}) == len(kept) == 3


def test_manifest_marks_duplicates(tmp_path):
    import subprocess, sys
    import pandas as pd
    raw = tmp_path / "raw"
    write_bv(raw / "100.1" / "DLPFC045", "rest", 10)
    write_bv(raw / "100.2" / "DLPFC045", "rest", 10)
    write_bv(raw / "100.1" / "uh9999" / "stage 2 surgery", "contact1", 5)
    script = Path(__file__).resolve().parent.parent / "scripts" / "make_manifest.py"
    subprocess.run([sys.executable, str(script), str(raw), str(tmp_path / "m.csv")], check=True, capture_output=True)
    m = pd.read_csv(tmp_path / "m.csv", keep_default_na=False)
    assert list(m["duplicate"]) == ["", "", "duplicate"] #sorted by path: 100.1/DLPFC045, 100.1/uh9999, 100.2/DLPFC045
    assert m["duplicate_of"][2] == "DLPFC045/rest.vhdr"
    assert list(m["condition"]) == ["rest", "contact", "rest"] and list(m["cohort"]) == ["dlpfc", "", "dlpfc"]

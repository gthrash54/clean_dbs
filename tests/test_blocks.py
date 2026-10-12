from pathlib import Path
import pytest
from clean_dbs.block_info import describe_block, output_path, missing_files
from clean_dbs.find_blocks import find_tdt_blocks, find_bv_recordings, find_all_blocks


def test_describe_bv_stage_with_space():
    info = describe_block(Path("inventory") / "uh9999" / "stage 2 surgery" / "STN_dbs1.vhdr")
    assert info == {"patient":"uh9999", "stage":"stage2", "block":"STN_dbs1", "format":"bv"}


def test_describe_tdt_stage_with_note():
    info = describe_block(Path("inventory") / "ug9999" / "stage2 surgery (sedated, not GA)" / "clinical_dbs1")
    assert info == {"patient":"ug9999", "stage":"stage2", "block":"clinical_dbs1", "format":"tdt"}


def test_describe_stage1_uppercase_vhdr():
    info = describe_block(Path("uh9999") / "Stage1 Surgery" / "rec.VHDR") #case shouldn't matter
    assert info["stage"] == "stage1" and info["format"] == "bv" and info["block"] == "rec"


def test_describe_bad_folder():
    with pytest.raises(ValueError, match="stage"):
        describe_block(Path("uh9999") / "surgery day" / "clinical_dbs1")


def test_output_path(tmp_path):
    out = output_path(Path("raw") / "uh9999" / "stage 2 surgery" / "STN_dbs1.vhdr", tmp_path)
    assert out == tmp_path / "uh9999" / "stage2" / "STN_dbs1.h5"
    out = output_path(Path("raw") / "ug9999" / "stage2 surgery (sedated, not GA)" / "clinical_dbs1", tmp_path)
    assert out == tmp_path / "ug9999" / "stage2" / "clinical_dbs1.h5"


def make_tree(root):
    #two patients, tdt blocks (folders with a .tsq) and brainvision recordings (.vhdr files), some marked do not use
    tdt = root / "uh9999" / "stage1 surgery"
    for block in ("block1", "block2_DNU", "block3"):
        (tdt / block).mkdir(parents=True)
        (tdt / block / f"{block}.tsq").write_bytes(b"")
    (tdt / "block1" / "block1.tev").write_bytes(b"x")
    (tdt / "notes").mkdir() #a folder without a .tsq isn't a block
    bv = root / "ug9999" / "stage 2 surgery"
    bv.mkdir(parents=True)
    for name in ("rec1", "rec2_dnu"):
        (bv / f"{name}.vhdr").write_text(f"DataFile={name}.eeg\nMarkerFile={name}.vmrk\n")
    (bv / "rec1.eeg").write_bytes(b"x")
    return tdt, bv


def test_find_blocks(tmp_path):
    tdt, bv = make_tree(tmp_path)
    assert find_tdt_blocks(tmp_path) == [tdt / "block1", tdt / "block3"]
    assert find_tdt_blocks(tmp_path, include_dnu=True) == [tdt / "block1", tdt / "block2_DNU", tdt / "block3"]
    assert find_bv_recordings(tmp_path) == [bv / "rec1.vhdr"]
    assert find_bv_recordings(tmp_path, include_dnu=True) == [bv / "rec1.vhdr", bv / "rec2_dnu.vhdr"]
    assert find_all_blocks(tmp_path) == sorted([tdt / "block1", tdt / "block3", bv / "rec1.vhdr"])
    assert find_all_blocks(tmp_path / "ug9999") == [bv / "rec1.vhdr"] #can point at just one patient


def test_missing_files(tmp_path):
    tdt, bv = make_tree(tmp_path)
    assert missing_files(tdt / "block1") == []
    assert missing_files(tdt / "block3") == [".tev"] #only the index file was copied
    assert missing_files(bv / "rec1.vhdr") == []
    assert missing_files(bv / "rec2_dnu.vhdr") == ["rec2_dnu.eeg"]

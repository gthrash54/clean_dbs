from clean_dbs import find_all_blocks


def test_imaging_folders_are_skipped(tmp_path):
    #scans left in an Imaging folder (any capitals, at any depth) are never recordings
    keep = tmp_path / "uh3001" / "stage1 surgery" / "block1"
    keep.mkdir(parents=True)
    (keep / "block1.tsq").touch()
    for skip in ("Imaging/x", "stage1 surgery/imaging/y"):
        d = tmp_path / "uh3001" / skip
        d.mkdir(parents=True)
        (d / "a.tsq").touch()
        (d / "a.vhdr").touch()
    assert find_all_blocks(tmp_path) == [keep]
    assert find_all_blocks(tmp_path / "uh3001" / "stage1 surgery" / "block1") == [keep] #a block folder passed straight in still counts

#cleaning every cohort on Box, not just uh/UG: names that can't reach the outputs, several raw folders, excluded recordings, cohorts
#sharing one out folder, and the copy script. "Zorblax" and "Quuxly" are made-up stand-ins for a person's name; never put a real one here
import os
import re
import shutil
import subprocess
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import pytest
from bv_fake import write_bv, write_mcs
from test_clean_script import run, SCRIPTS
from clean_dbs import privacy
from clean_dbs.cohorts import COHORTS
from clean_dbs.block_info import describe_block, output_path, source_path, names_to_check, has_study_id

BOX = Path("box") / "100.0000001"
FS = 100000.0


@pytest.fixture(autouse=True)
def no_renames():
    #renames are module state, so every test starts and ends without any
    privacy.reset()
    yield
    privacy.reset()


def renames(out, rows, words=()):
    (out / "_cases").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["study_id", "original", "renamed"]).to_csv(out / "_cases" / "renames.csv", index=False)
    if words:
        (out / "_cases" / "allowed_words.txt").write_text("\n".join(words) + "\n")


def bv(folder, name, n=2000):
    return write_bv(folder, name, [("DBS1", "µV"), ("17", "µV")], np.zeros((2, n), dtype=np.float32), 10000.0, [])


def test_words_and_what_passes():
    assert privacy.words_of("HandOCAndBuzz2") == ["hand", "oc", "and", "buzz"]
    assert privacy.words_of("DLPFC045_rest_eyesclosed") == ["rest", "eyesclosed"] #study IDs aren't words
    assert privacy.unknown_words("lfp1_Cam1") == [] #tdt's camera files
    assert privacy.unknown_words("BRAINaim1_events") == [] #read at the case change when the camel case reading doesn't work
    assert privacy.unknown_words("stage2 surgery (sedated, not GA)") == []
    assert privacy.unknown_words("Contents of DLPFC110 folder") == []
    assert privacy.unknown_words("rest_Zorblax2") == ["zorblax"]
    assert privacy.unknown_words("Zorblaxrest") == ["zorblaxrest"] #glued to a known word it's still unknown
    assert privacy.unknown_words("Zörblax") == ["zörblax"] #letters outside a-z stay one word
    assert privacy.unknown_words("rest_Zorblax", allowed={"zorblax"}) == []
    assert privacy.mask("DLPFC045/Zorblax/rest1.vhdr") == "DLPFC045/<withheld>/rest1.vhdr"
    assert privacy.mask("DLPFC045/Misc/Quuxly.vhdr") == "DLPFC045/Misc/<withheld>"
    #no first names hiding in the word list itself: every word is checked against the cohort names by the dry runs (see docs)
    assert all(w == w.lower() and w.isalpha() for w in privacy.WORDS)


def test_renames_reach_every_name(tmp_path):
    out = tmp_path / "out"
    renames(out, [["DLPFC045", "Zorblax", "session_a"], ["DLPFC045", "Quuxly rest", "rest_2"], ["", "Quuxly", "rest_3"]])
    block = BOX / "DLPFC045" / "Zorblax" / "Quuxly rest.vhdr"
    assert privacy.use_out_folder(out) == (3, 0)
    info = describe_block(block)
    assert info["session"] == "session_a" and info["block"] == "rest_2" and info["condition"] == "rest"
    assert output_path(block, "o") == Path("o") / "DLPFC045" / "session_a" / "rest_2.h5"
    assert source_path(block) == "DLPFC045/session_a/rest_2.vhdr"
    #a row with no study ID applies to every folder; a row for another ID doesn't
    assert describe_block(BOX / "BI020" / "Quuxly.vhdr")["block"] == "rest_3"
    assert describe_block(BOX / "BI020" / "Zorblax" / "x.vhdr")["session"] == "zorblax"
    assert privacy.scrub("can't read Quuxly rest.eeg") == "can't read rest_2.eeg"
    privacy.use_out_folder(tmp_path / "elsewhere") #another out folder without renames: nothing renamed
    assert describe_block(block)["block"] == "Quuxly rest"


def test_a_renamed_value_has_to_pass_too(tmp_path):
    renames(tmp_path, [["", "Zorblax", "Quuxly_rest"]])
    with pytest.raises(ValueError, match="doesn't pass"):
        privacy.use_out_folder(tmp_path)
    renames(tmp_path, [["", "Zorblax", "rest/2"]])
    with pytest.raises(ValueError, match="/"):
        privacy.use_out_folder(tmp_path)
    renames(tmp_path, [["", "Zorblax", ""]], words=["QuuxReach  # a task", "SNRx"]) #a row not filled in yet is ignored
    assert privacy.use_out_folder(tmp_path) == (0, 2)
    assert privacy.unknown_words("QuuxReach_1") == [] and privacy.unknown_words("quux") == ["quux"] #a word is allowed whole, never in pieces


def test_names_to_check_covers_folders_files_and_sidecars(tmp_path):
    folder = tmp_path / "raw" / "100.0000001" / "MC9999" / "Zorblax"
    vhdr = bv(folder, "Quuxly_rest")
    (folder / "Quuxly_results.txt").write_text("1\n") #a sidecar name gets into /sidecars
    found = {(n["kind"], n["original"]) for n in names_to_check(vhdr)}
    assert found == {("folder", "Zorblax"), ("recording", "Quuxly_rest"), ("sidecar", "Quuxly_results")}
    assert all(n["study_id"] == "MC9999" for n in names_to_check(vhdr))
    renames(tmp_path / "out", [["MC9999", "Zorblax", "awake"], ["MC9999", "Quuxly_rest", "rest"], ["MC9999", "Quuxly_results", "results"]])
    privacy.use_out_folder(tmp_path / "out")
    assert names_to_check(vhdr) == []


def test_source_path_never_has_the_number_folder():
    assert source_path(BOX / "BPP012.vhdr") == "BPP012.vhdr" #a recording straight in a number folder
    assert source_path(BOX / "MOTIF" / "4and5_10.vhdr") == "MOTIF/4and5_10.vhdr"
    assert source_path(Path("box") / "100.0000001_abcd" / "MOTIF" / "x.vhdr") == "MOTIF/x.vhdr" #number folder with a suffix
    assert source_path(BOX / "PP012" / "CHR003" / "1-2-3.vhdr") == "PP012/CHR003/1-2-3.vhdr" #unchanged
    assert source_path(BOX / "uh9999" / "stage1 surgery" / "clinical_dbs1") == "uh9999/stage1 surgery/clinical_dbs1" #unchanged
    assert not has_study_id(BOX / "MOTIF" / "4and5_10.vhdr") and has_study_id(BOX / "BPP012.vhdr")
    assert not has_study_id(BOX / "DLPFC080_BPP004" / "x.vhdr") #two IDs in one folder name: whose is it?


def make_raw(root, number="100.0000001"):
    #a cleared DLPFC recording, one in a folder named after a person, one in a folder with no study ID
    raw = root / number
    bv(raw / "DLPFC045", "rest1")
    bv(raw / "DLPFC046" / "Zorblax", "Quuxly rest")
    bv(raw / "MOTIF", "4and5_10")
    return raw


def test_manifest_refuses_names_until_renamed(tmp_path):
    raw, out = make_raw(tmp_path / "raw"), tmp_path / "out"
    manifest = out / "manifest_dlpfc.csv"
    p = run(SCRIPTS / "make_manifest.py", raw, manifest)
    assert p.returncode != 0 and not manifest.exists()
    assert "orblax" not in p.stdout + p.stderr and "uuxly" not in p.stdout + p.stderr #never printed
    assert "DLPFC046 2" in p.stdout #by study ID and count only
    todo = pd.read_csv(out / "_cases" / "names_to_check.csv", keep_default_na=False)
    assert set(todo["original"]) == {"Zorblax", "Quuxly rest"} and set(todo["renamed"]) == {""}
    #a person fills in "renamed" and keeps the rows in renames.csv
    todo["renamed"] = todo["original"].map({"Zorblax": "session_a", "Quuxly rest": "rest_2"})
    todo.to_csv(out / "_cases" / "renames.csv", index=False)
    p = run(SCRIPTS / "make_manifest.py", raw, manifest, "--chunk", 2)
    assert p.returncode == 0, p.stderr
    m = pd.read_csv(manifest, keep_default_na=False).set_index("block")
    assert m.loc["rest_2", "out"] == "DLPFC046/session_a/rest_2.h5" and m.loc["rest_2", "session"] == "session_a"
    #MOTIF has no study ID: listed with the reason, not dropped
    assert m.loc["4and5_10", "excluded"].startswith("no study ID") and "MOTIF/4and5_10.vhdr" in m.loc["4and5_10", "excluded"]
    assert "1 excluded" in p.stdout and "--no-cases" in p.stdout #no case table rows for these IDs
    text = m.drop(columns="path").to_csv() #the raw paths have them, the rest of the manifest doesn't
    assert "orblax" not in text and "uuxly" not in text


def test_exclude_unchecked_and_clean(tmp_path):
    raw, out = make_raw(tmp_path / "raw"), tmp_path / "out"
    manifest = out / "manifest_dlpfc.csv"
    p = run(SCRIPTS / "make_manifest.py", raw, manifest, "--exclude-unchecked")
    assert p.returncode == 0 and "orblax" not in p.stdout
    m = pd.read_csv(manifest, keep_default_na=False)
    held = m[m["excluded"].str.startswith("names")]
    assert len(held) == 1 and held["block"].iloc[0] == "<withheld>" and held["out"].iloc[0] == ""
    for task in range(len(m)):
        p = run(SCRIPTS / "clean.py", "--manifest", manifest, "--task", task, "--allow-dirty", "--no-cases", out)
        assert p.returncode == 0 and "orblax" not in p.stdout + p.stderr and "uuxly" not in p.stdout + p.stderr
    p = run(SCRIPTS / "clean.py", "--merge-reports", out, "--manifest", manifest)
    report = pd.read_csv(out / "clean_report_dlpfc.csv", keep_default_na=False) #named after the manifest
    assert sorted(report["status"]) == ["cleaned", "excluded", "excluded"]
    assert "orblax" not in (out / "clean_report_dlpfc.csv").read_text().lower()
    assert [p.relative_to(out).as_posix() for p in out.rglob("*.h5")] == ["DLPFC045/rest1.h5"]


def test_old_manifest_still_checked_by_clean(tmp_path):
    #a manifest made before the name check (no excluded column) can't push a name through: clean.py checks every block itself
    raw, out = make_raw(tmp_path / "raw"), tmp_path / "out"
    vhdr = raw / "DLPFC046" / "Zorblax" / "Quuxly rest.vhdr"
    manifest = tmp_path / "manifest.csv"
    pd.DataFrame([{"path":str(vhdr), "patient":"DLPFC046", "stage":"", "block":"x", "format":"bv", "missing":""}]).to_csv(manifest, index=False)
    p = run(SCRIPTS / "clean.py", "--manifest", manifest, "--task", 0, "--allow-dirty", "--no-cases", out)
    assert p.returncode == 0 and "excluded" in p.stdout and "orblax" not in p.stdout
    assert not list(out.rglob("*.h5"))


def test_renamed_names_never_reach_the_file(tmp_path):
    #recording, folder and stimulator program all named after a person: after renames.csv nothing in the .h5 or the report has the name
    rng = np.random.default_rng(1)
    gaps = rng.integers(70000, 120000, 20)
    ch1 = [(500, 60, -500, 60, 0, g, 0) for g in gaps]
    sync = [(1, 100, 0, 0, 0, 20 + g, 0) for g in gaps]
    folder = tmp_path / "raw" / "100.0000001" / "MC9999" / "Zorblax"
    folder.mkdir(parents=True)
    write_mcs(folder / "Zorblax_results.txt", ch1, sync)
    (folder / "Zorblax_1.1_ttl.txt").write_text("cond={'a' };\ncond_ms=[0.1 ];\n") #a dot in the name: its group name is made from the name
    onsets = 5000 + np.concatenate(([0], np.cumsum(120 + gaps[:-1]))) / 10
    markers = [m for s in onsets for m in (("Stimulus", "S 15", int(round(s))), ("Stimulus", "S 7", int(round(s + 10))))]
    write_bv(folder, "Quuxly_rest", [("DBS1", "µV"), ("17", "µV")], np.zeros((2, int(onsets[-1] + 20000)), dtype=np.float32), FS, markers)
    out = tmp_path / "out"
    renames(out, [["MC9999", "Zorblax", "awake"], ["MC9999", "Quuxly_rest", "rest"], ["MC9999", "Zorblax_results", "stim_results"],
                  ["MC9999", "Zorblax_1.1_ttl", "ttl_1.1"]])
    p = run(SCRIPTS / "clean.py", tmp_path / "raw", out, "--no-cases")
    assert p.returncode == 0, p.stdout + p.stderr
    h5 = out / "MC9999" / "awake" / "rest.h5"
    with h5py.File(h5, "r") as f:
        assert f["info"].attrs["block"] == "rest" and f["info"].attrs["session"] == "awake"
        assert f.attrs["source"] == "MC9999/awake/rest.vhdr"
        assert f["sidecars/stim_results"].attrs["file"] == "stim_results.txt" and f["sidecars/stim_results"].attrs["matched"] == 1
        assert set(v.decode() for v in f["stim/sidecar"][:]) == {"stim_results.txt"}
        assert sorted(f["sidecars"]) == ["stim_results", "ttl_1_1"] and f["sidecars/ttl_1_1"].attrs["file"] == "ttl_1.1.txt"
    for made in [h5, out / "clean_report.csv"]:
        data = made.read_bytes().lower()
        assert b"zorblax" not in data and b"quuxly" not in data, made.name


def test_several_raw_folders_and_cohorts_sharing_an_out_folder(tmp_path):
    out = tmp_path / "out"
    a = tmp_path / "scratch" / "lfp" / "100.0000001"
    b = tmp_path / "scratch" / "pp" / "100.0000002"
    bv(a / "LFP034", "lfp1")
    bv(b / "PP012", "pp1")
    p = run(SCRIPTS / "make_manifest.py", a, b, out / "manifest_lfp.csv")
    assert p.returncode == 0 and "2 recordings found" in p.stdout
    #the PP copy carries LFP034 nested inside it, a different recording that lands on the same file: refused
    bv(b / "PP012" / "LFP034", "lfp1", n=3000)
    p = run(SCRIPTS / "make_manifest.py", b, out / "manifest_pp.csv")
    assert p.returncode != 0 and "LFP034/lfp1.h5" in p.stderr and "manifest_lfp.csv" in p.stderr
    #the same size: the same recording copied twice, cleaned once by the other manifest
    shutil.copy(a / "LFP034" / "lfp1.eeg", b / "PP012" / "LFP034" / "lfp1.eeg")
    p = run(SCRIPTS / "make_manifest.py", b, out / "manifest_pp.csv")
    assert p.returncode == 0, p.stderr
    m = pd.read_csv(out / "manifest_pp.csv", keep_default_na=False).set_index("out")
    assert m.loc["LFP034/lfp1.h5", "duplicate"] == "duplicate" and m.loc["LFP034/lfp1.h5", "duplicate_of"].startswith("manifest_lfp.csv")
    assert m.loc["LFP034/lfp1.h5", "parent_id"] == "PP012"


def test_check_copy_takes_several_folders_and_hides_names(tmp_path):
    raw = make_raw(tmp_path / "raw")
    (raw / "DLPFC046" / "Zorblax" / "Quuxly rest.eeg").unlink()
    other = tmp_path / "raw2" / "100.0000003"
    bv(other / "BI020", "lfp1")
    p = run(SCRIPTS / "check_copy.py", raw, other)
    assert p.returncode == 1 and "4 recordings checked, 1 with problems" in p.stdout
    assert "DLPFC046/<withheld>/<withheld>: missing <withheld>" in p.stdout and "orblax" not in p.stdout


def test_copy_script_cohort_list_matches_cohorts():
    text = (SCRIPTS / "copy_cohorts.sh").read_text()
    listed = re.search(r'^COHORTS="([^"]+)"', text, re.M).group(1).split()
    assert sorted(listed) == sorted(COHORTS)
    assert " sync " not in text.replace("never sync", "") and "rclone sync" not in text and "rclone move" not in text


def fake_box(root):
    #a made-up Box in the shape the real one has: number folders (some with a suffix), study folders, a patient-named folder at the top
    box = root / "UAB DBS Data" / "UAB DBS Raw Data"
    files = ["100.0000001/DLPFC045/rest1.vhdr", "100.0000001/DLPFC045/rest1.eeg", "100.0000001/DLPFC045/Imaging/scan.dcm",
             "100.0000001/DLPFC045/Images/MR AX/1.dcm", "100.0000001/DLPFC045/consent.pdf", "100.0000001/DLPFC045/BPP004/bpp1.vhdr",
             "100.0000002_abcd/dlpfc046/rest1.vhdr", "100.0000003/DLPFC080_BPP004/x.vhdr", "100.0000003/BPP012.vhdr", "100.0000003/BPP012.eeg",
             "100.0000003/notes.txt", "100.0000003/DLPFC048 Zorblax 12/x.vhdr", "100.0000004/MOTIF/4and5_10.vhdr", "100.0000004/Other/x.vhdr", "100.0000005/PP012/LFP034/lfp1.vhdr",
             "Zorblax Quuxly/DLPFC047/rest1.vhdr"]
    for f in files:
        (box / f).parent.mkdir(parents=True, exist_ok=True)
        (box / f).write_text("x")
    return box


@pytest.mark.skipif(shutil.which("rclone") is None, reason="needs rclone")
def test_copy_script_takes_only_study_folders(tmp_path):
    fake_box(tmp_path)
    env = {**os.environ, "REMOTE":"fakebox", "RCLONE_CONFIG_FAKEBOX_TYPE":"local", "DEST":str(tmp_path / "dest"), "RCLONE_CONFIG":str(tmp_path / "none.conf")}
    def copy(*args):
        return subprocess.run(["bash", str(SCRIPTS / "copy_cohorts.sh"), *args], cwd=tmp_path, env=env, capture_output=True, text=True)
    p = copy("--list-only", "dlpfc")
    assert p.returncode == 0 and "dlpfc: 3 folders and 0 loose files" in p.stdout and not (tmp_path / "dest" / "dlpfc" / "100.0000001").exists()
    assert "100.000" not in p.stdout + p.stderr #number folders are never printed
    assert copy("dlpfc", "bpp", "motif").returncode == 0
    got = sorted(f.relative_to(tmp_path / "dest").as_posix() for f in (tmp_path / "dest").rglob("*") if f.is_file() and not f.name.startswith(("_", ".")))
    assert got == ["bpp/100.0000003/BPP012.eeg", "bpp/100.0000003/BPP012.vhdr",
                   "dlpfc/100.0000001/DLPFC045/BPP004/bpp1.vhdr", "dlpfc/100.0000001/DLPFC045/rest1.eeg", "dlpfc/100.0000001/DLPFC045/rest1.vhdr",
                   "dlpfc/100.0000002_abcd/dlpfc046/rest1.vhdr", "dlpfc/100.0000003/DLPFC080_BPP004/x.vhdr",
                   "motif/100.0000004/MOTIF/4and5_10.vhdr"]
    p = copy("--plain-only", "--list-only", "dlpfc")
    assert "dlpfc: 2 folders" in p.stdout #the suffixed number folder left out
    assert copy("uh").returncode != 0 and copy("nosuch").returncode != 0

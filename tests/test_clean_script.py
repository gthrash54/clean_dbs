import os
import subprocess
import sys
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from test_windows import write_tdt
from clean_dbs.save_h5 import stale_reason
from clean_dbs.version import CLEAN_VERSION

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"


def run(*args, env=None):
    #Input: arguments for a script in scripts/. Returns: the finished process (exit code, printed text)
    return subprocess.run([sys.executable, *map(str, args)], capture_output=True, text=True, env={**os.environ, **(env or {})})


def make_raw(tmp_path):
    #one good tdt block, and one missing its .tev (fails), under a made-up number folder
    raw = tmp_path / "raw" / "1234"
    rng = np.random.default_rng(5)
    write_tdt(raw / "uh9999" / "stage1 surgery" / "good", {"ecos":(48828.125, rng.standard_normal((2, 1024 * 10)).astype(np.float32), 1024)})
    bad = write_tdt(raw / "uh9999" / "stage1 surgery" / "bad", {"ecos":(48828.125, np.zeros((2, 1024 * 4), dtype=np.float32), 1024)})
    (bad / "b.tev").unlink()
    return raw


def add_cases(out, status="conflict"):
    (out / "_cases").mkdir(parents=True, exist_ok=True)
    #capital study ID like on Box, lowercase folder: the match has to ignore case
    pd.DataFrame([{"study_id":"UH9999", "implant":"dual", "target1":"STN", "target2":"GPi", "ecog":"", "status":status}]).to_csv(out / "_cases" / "cases.csv", index=False)


def report(out):
    return pd.read_csv(out / "clean_report.csv", keep_default_na=False).set_index("block")


def test_plain_mode_statuses(tmp_path):
    raw, out = make_raw(tmp_path), tmp_path / "out"
    p = run(SCRIPTS / "clean.py", raw, out)
    assert p.returncode != 0 and "_cases" in p.stderr #no case tables: refuses to start
    assert run(SCRIPTS / "clean.py", raw, tmp_path / "out_nocases", "--no-cases").returncode == 1 #runs, but a block failed
    add_cases(out)
    p = run(SCRIPTS / "clean.py", raw, out)
    assert p.returncode == 1 #the block with no .tev failed, so the run says so
    r = report(out)
    assert r.loc["good", "status"] == "cleaned" and r.loc["bad", "status"] == "failed"
    assert str(tmp_path) not in r.loc["bad", "reason"] and "uh9999/stage1 surgery/bad" in r.loc["bad", "reason"] #no full paths
    h5 = out / "uh9999" / "stage1" / "good.h5"
    with h5py.File(h5, "r") as f:
        assert f.attrs["clean_version"] == CLEAN_VERSION
        assert f["info"].attrs["case_status"] == "conflict"
    from clean_dbs import build_index
    idx = build_index(out)
    assert idx.loc[0, "implant"] == "dual" and idx.loc[0, "clean_version"] == CLEAN_VERSION #UH9999 in cases.csv joins uh9999's files
    run(SCRIPTS / "clean.py", raw, out)
    assert report(out).loc["good", "status"] == "skipped"
    with h5py.File(h5, "a") as f: #made by an older version of the code
        f.attrs["clean_version"] = CLEAN_VERSION - 1
    assert "older" in stale_reason(h5)
    run(SCRIPTS / "clean.py", raw, out)
    r = report(out)
    assert r.loc["good", "status"] == "recleaned" and "older" in r.loc["good", "reason"]
    h5.write_bytes(b"not an hdf5 file") #half-copied or corrupt
    run(SCRIPTS / "clean.py", raw, out)
    assert report(out).loc["good", "reason"].startswith("file won't open")
    run(SCRIPTS / "clean.py", raw, out, "--force")
    assert report(out).loc["good", "status"] == "recleaned" and report(out).loc["good", "reason"] == "--force"
    with h5py.File(h5, "a") as f:
        del f.attrs["clean_version"] #made before clean_version existed
    assert stale_reason(h5) == "made before clean_version existed"


def test_manifest_mode_and_merge(tmp_path):
    raw, out = make_raw(tmp_path), tmp_path / "out"
    add_cases(out, status="")
    manifest = tmp_path / "manifest.csv"
    p = run(SCRIPTS / "make_manifest.py", raw, manifest)
    assert p.returncode == 0
    m = pd.read_csv(manifest, keep_default_na=False)
    assert {"signal_gb", "est_peak_gb", "cohort", "condition"} <= set(m.columns) and len(m) == 2
    #no git commit to trace the files to: refuses unless --allow-dirty
    p = run(SCRIPTS / "clean.py", "--manifest", manifest, "--task", 0, out, env={"GIT_DIR":str(tmp_path / "no_git")})
    assert p.returncode != 0 and "allow-dirty" in p.stderr
    p = run(SCRIPTS / "clean.py", "--manifest", manifest, "--task", 5, out) #past the end: nothing to do, not an error
    assert p.returncode == 0 and "past the end" in p.stdout
    for task in (0, 1):
        run(SCRIPTS / "clean.py", "--manifest", manifest, "--task", task, "--allow-dirty", out)
    run(SCRIPTS / "clean.py", "--manifest", manifest, "--task", 0, "--chunk", 2, "--allow-dirty", out) #the same rows again, another chunk size
    names = sorted(p.name for p in (out / "_reports").glob("*.csv"))
    assert len(names) == 3 and all("_chunk" in n for n in names)
    (out / "_reports" / "clean_report_7_chunk1_deadbeef.csv").write_text("row,block,status,reason,finished\n0,x,cleaned,,2000-01-01\n") #another manifest
    assert run(SCRIPTS / "clean.py", "--merge-reports", out).returncode != 0 #needs --manifest
    p = run(SCRIPTS / "clean.py", "--merge-reports", out, "--manifest", manifest)
    assert p.returncode == 0 and "another manifest" in p.stdout and "older report row" in p.stdout
    r = pd.read_csv(out / "clean_report.csv", keep_default_na=False)
    assert r["row"].tolist() == [0, 1] #one row each, the newest try
    assert sorted(r["status"]) == ["failed", "skipped"] #the chunk-2 run came last: good was already there
    with h5py.File(out / "uh9999" / "stage1" / "good.h5", "r") as f:
        assert f["info"].attrs["case_status"] == "unchecked" #row there, cross-check never run


def test_source_bytes_counts_brainvision_files(tmp_path):
    from test_windows import write_bv
    from clean_dbs.save_h5 import _provenance
    vhdr = write_bv(tmp_path / "ug9999" / "stage2 surgery", np.zeros((2, 100)))
    expected = sum((vhdr.parent / n).stat().st_size for n in ("rec.vhdr", "rec.eeg", "rec.vmrk"))
    assert _provenance(vhdr)["source_bytes"] == expected


def test_merge_warns_on_rows_outside_manifest(tmp_path):
    import hashlib
    manifest = tmp_path / "manifest.csv"
    manifest.write_text("path,patient,stage,block,format,missing\n/x/uh9999/stage1 surgery/a,uh9999,stage1,a,tdt,\n")
    tag = hashlib.sha1(manifest.read_bytes()).hexdigest()[:8]
    (tmp_path / "out" / "_reports").mkdir(parents=True)
    (tmp_path / "out" / "_reports" / f"clean_report_0_chunk1_{tag}.csv").write_text(
        "row,block,status,reason,finished\n0,a,cleaned,,2026-01-01T00:00:00\n99,zz,cleaned,,2026-01-01T00:00:00\n")
    p = run(SCRIPTS / "clean.py", "--merge-reports", tmp_path / "out", "--manifest", manifest)
    assert p.returncode == 0 and "aren't in the manifest" in p.stdout
    assert pd.read_csv(tmp_path / "out" / "clean_report.csv")["row"].tolist() == [0]

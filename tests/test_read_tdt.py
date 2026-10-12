from pathlib import Path
import numpy as np
import pytest
from clean_dbs.read_tdt import _store_settings, _to_volts

#made-up StoresListing.txt copying the layout Synapse writes (spaces before the colons, a device then its stores)
LISTING = """Experiment: test_rig
Subject: uh9999
Time: 9:19:18am

Object ID : RZ2(1) - RZn Processor
 Rate     : 48828.1 Hz
 Store ID : Note

Object ID : ecog - Stream Data Storage
 Store ID : ecog
  Format  : Integer-16
  Scale   : Micro
  Rate    : 48828.1 Hz

Object ID : spike - Stream Data Storage
 Store ID : spie
  Format  : Integer-16
  Scale   : Milli

Object ID : eStim1 - Electrical Stim Driver
 Store ID : PeA/
  Source  : Strobe input
"""


def test_store_settings(tmp_path):
    (tmp_path / "StoresListing.txt").write_text(LISTING)
    settings = _store_settings(Path(tmp_path))
    assert settings["ecog"] == {"device":"ecog", "format":"Integer-16", "scale":"Micro"}
    assert settings["spie"] == {"device":"spike", "format":"Integer-16", "scale":"Milli"}
    assert settings["PeA_"] == {"device":"eStim1"} #"PeA/" is renamed the way tdt names it
    assert settings["Note"] == {"device":"RZ2(1)"}


def test_store_settings_missing(tmp_path):
    assert _store_settings(Path(tmp_path)) == {} #older blocks may not have the file


def test_to_volts_integer_stores():
    data = np.array([[100, -200]], dtype=np.int16)
    volts, units = _to_volts(data, {"scale":"Micro"}, "ecog", Path("block1"))
    assert units == "V" and volts.dtype == np.float32
    np.testing.assert_allclose(volts, [[100e-6, -200e-6]], rtol=1e-6)
    volts, _ = _to_volts(data, {"scale":"Milli"}, "spie", Path("block1"))
    np.testing.assert_allclose(volts, [[0.1, -0.2]], rtol=1e-6)


def test_to_volts_float_is_unchanged():
    data = np.array([[0.5, -0.25]], dtype=np.float32)
    for setting in ({}, {"scale":"Unity"}): #float stores are already volts, with or without a listing
        volts, units = _to_volts(data, setting, "ecos", Path("block1"))
        np.testing.assert_array_equal(volts, data)
        assert units == "V"


def test_to_volts_integer_without_scale_fails():
    with pytest.raises(ValueError, match="StoresListing"):
        _to_volts(np.array([[1]], dtype=np.int16), {}, "ecog", Path("block1"))


def small_block(tmp_path, start=1.7e9):
    from test_windows import write_tdt
    rng = np.random.default_rng(4)
    ecos = rng.standard_normal((2, 1024 * 20)).astype(np.float32)
    return write_tdt(tmp_path / "uh9999" / "stage1 surgery" / "block1", {"ecos":(48828.125, ecos, 1024)}, start=start)


def test_truncated_tev_raises(tmp_path):
    from clean_dbs.read_tdt import read_tdt, tdt_windows, tev_problem
    block = small_block(tmp_path)
    assert tev_problem(block) == ""
    tev = block / "b.tev"
    tev.write_bytes(tev.read_bytes()[:-4]) #one sample short, tdt would quietly fill it with zeros
    assert "cut short" in tev_problem(block)
    for read in (read_tdt, tdt_windows):
        with pytest.raises(ValueError, match="cut short"):
            read(block)


def test_start_time_is_chicago_whatever_the_machine(tmp_path, monkeypatch):
    import time
    from clean_dbs.read_tdt import read_tdt
    block = small_block(tmp_path) #1.7e9 s = 2023-11-14 22:13:20 UTC = 16:13:20 in Chicago (CST)
    for tz in ("UTC", "America/Chicago", "Asia/Tokyo"):
        monkeypatch.setenv("TZ", tz)
        time.tzset()
        assert read_tdt(block)["info"]["start"] == "2023-11-14T16:13:20"
    monkeypatch.undo()
    time.tzset()


def test_duration_from_streams_when_no_stop(tmp_path):
    from clean_dbs.read_tdt import read_tdt, tdt_windows
    block = small_block(tmp_path)
    info = read_tdt(block)["info"]
    assert info["duration"] == pytest.approx(1024 * 20 / 48828.125) #longest stream, not the header
    assert info["duration_header"] == pytest.approx(1024 * 19 / 48828.125 + 1) #write_tdt puts the stop 1 s after the last packet
    tsq = block / "b.tsq"
    tsq.write_bytes(tsq.read_bytes()[:-40]) #drop the stop marker, like a block that didn't end cleanly
    with pytest.warns(Warning):
        info = read_tdt(block)["info"]
    assert info["duration"] == pytest.approx(1024 * 20 / 48828.125) and np.isnan(info["duration_header"])
    assert tdt_windows(block)[0]["info"]["duration"] == pytest.approx(1024 * 20 / 48828.125)

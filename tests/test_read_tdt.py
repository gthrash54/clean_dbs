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

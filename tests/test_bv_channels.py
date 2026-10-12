from pathlib import Path
import numpy as np
import pytest
from clean_dbs.bv_channels import classify_channel, channel_types, cohort_of, normalize
from clean_dbs.read_bv import read_bv
from bv_fake import write_bv


@pytest.mark.parametrize("name, kind, group, contact", [
    ("DBS3", ("lfp", "lead"), "lead1", "3"), ("DBS12", ("lfp", "lead"), "lead2", "4"), ("LFP8", ("lfp", "lead"), "lead1", "8"),
    ("DB S8", ("lfp", "lead"), "lead1", "8"), ("DBS2_L", ("lfp", "lead"), "lead_left", "2"), ("DBS1 Ch7", ("lfp", "lead"), "lead1", "1"),
    ("ECOG_A3", ("lfp", "ecog"), "ecog_anterior", "3"), ("ecog_p6", ("lfp", "ecog"), "ecog_posterior", "6"), ("ECOGp2", ("lfp", "ecog"), "ecog_posterior", "2"),
    ("ecog_4a", ("lfp", "ecog"), "ecog_anterior", "4"), ("ECOG5", ("lfp", "ecog"), "ecog", "5"),
    ("Cz", ("eeg", "scalp"), "eeg", ""), ("Y1 AF7", ("eeg", "scalp"), "eeg", ""), ("FCz", ("eeg", "scalp"), "eeg", ""),
    ("Microphone", ("audio", "mic"), "audio", ""), ("mic", ("audio", "mic"), "audio", ""),
    ("AO_MICRO_SPIKE", ("mer", "micro_spike"), "mer", ""), ("MACRO_RAW", ("mer", "macro_raw"), "mer", ""), ("Macro LFP", ("mer", "macro_lfp"), "mer", ""),
    ("NO_MICRO_RAW", ("unused", "empty_mer_slot"), "unused", ""), ("", ("unused", "empty_slot"), "unused", ""),
    ("Accelerometer_X", ("accel", ""), "accel", "x"), ("JOY_Y", ("joystick", ""), "joystick", "y"), ("Dynamometer", ("force", ""), "force", ""),
    ("Aux2", ("aux", ""), "aux", ""), ("photo_sensor", ("aux", "photo_sensor"), "aux", ""),
])
def test_known_names(name, kind, group, contact):
    c = classify_channel(name)
    assert (c["type"], c["subtype"]) == kind and c["group"] == group and c["contact"] == contact


@pytest.mark.parametrize("name, muscle, side", [
    ("FDI_contra", "first dorsal interosseous", "contra"), ("FDI_Conta", "first dorsal interosseous", "contra"),
    ("EMG FCR_Ipsi", "flexor carpi radialis", "ipsi"), ("Contra_FCR", "flexor carpi radialis", "contra"), ("GASTROC LOWER", "gastrocnemius (lower)", ""),
    ("FCU_Flex", "flexor carpi ulnaris (flex)", ""), ("fdi_left", "first dorsal interosseous", "left"), ("OO", "orbicularis oris", ""),
])
def test_emg(name, muscle, side):
    c = classify_channel(name)
    assert c["type"] == "emg" and c["muscle"] == muscle and c["side"] == side


@pytest.mark.parametrize("name", ["Ch12", "17", "Leg", "Mouth", "SCS3", "EP2", "G7", "FRONTAL1 Ch15", "FCI_ipsi", "knee_right"])
def test_generic_names_stay_unknown(name):
    #nothing in these names says what they are, so we never guess
    assert classify_channel(name)["type"] == "unknown"


def test_cohort_override_only_where_listed():
    assert classify_channel("Ch12", "BV")["type"] == "eeg" #64-ch HD-EEG cohort
    assert classify_channel("Ch12", "BI")["type"] == "unknown"
    assert "cohort override" in classify_channel("Ch12", "ORBV")["source"]


def test_cohort_of():
    assert cohort_of(Path("x") / "uh3021" / "stage 2 surgery" / "a.vhdr") == "UH"
    assert cohort_of(Path("x") / "PT001" / "DIR002" / "a.vhdr") == "PT" #subfolder that looks like an ID isn't a cohort
    assert cohort_of(Path("x") / "BI012" / "LFP" / "a.vhdr") == "BI" #a plain LFP folder isn't the LFP cohort
    assert cohort_of(Path("x") / "MOTIF" / "a.vhdr") == "MOTIF"
    assert cohort_of(Path("x") / "notes" / "a.vhdr") == ""
    assert normalize("FCR_contra_BB") == "fcr_contra"


def make_recording(folder, dated=True):
    channels = [("DBS1", "µV"), ("DBS2", "µV"), ("ECOG_A1", "µV"), ("FDI_contra", "µV"), ("mic", "µV"), ("Ch5", "µV"), ("Aux1", "BIP")]
    data = np.arange(len(channels))[:, None] * np.ones((1, 1000)) + 1 #channel k holds k+1 everywhere
    return write_bv(folder, "rec", channels, data, 1000.0, [("Stimulus", "S 15", 100), ("Stimulus", "S 13", 101)], dated=dated)


def test_channel_types_and_units(tmp_path):
    rows = channel_types(make_recording(tmp_path / "uh9999" / "stage 2 surgery"))
    assert [r["stream"] for r in rows] == ["lead", "lead", "ecog", "emg", "audio", "unknown", "aux_BIP"]
    assert rows[-1]["units"] == "BIP"


def test_read_bv_splits_streams(tmp_path):
    r = read_bv(make_recording(tmp_path / "uh9999" / "stage 2 surgery"))
    s = r["streams"]
    assert list(s) == ["lead", "ecog", "emg", "audio", "unknown", "aux_BIP"]
    np.testing.assert_array_equal(s["lead"]["channels"], [1, 2]) #position in the .vhdr, same numbers the impedance checks use
    np.testing.assert_allclose(s["lead"]["data"][:, 0], [1e-6, 2e-6], rtol=1e-5) #uV -> V
    assert s["lead"]["type"] == "lfp" and s["lead"]["subtype"] == "lead" and s["lead"]["units"] == "V"
    assert s["emg"]["channel_muscle"] == ["first dorsal interosseous"] and s["emg"]["channel_side"] == ["contra"]
    #a unit MNE doesn't know comes through unscaled and keeps its own label, never "V"
    assert s["aux_BIP"]["units"] == "BIP" and s["aux_BIP"]["channel_units"] == ["BIP"]
    np.testing.assert_allclose(s["aux_BIP"]["data"][0, 0], 7.0, rtol=1e-5)
    assert r["info"]["cohort"] == "uh" #lowercase, like cohorts.py
    assert read_bv(make_recording(tmp_path / "old"), split=False)["streams"]["dbs"]["data"].shape == (7, 1000) #old single stream still there


def test_read_bv_no_date(tmp_path):
    #no dated New Segment: start falls back to blank instead of crashing
    r = read_bv(make_recording(tmp_path / "nodate", dated=False))
    assert r["info"]["start"] == ""

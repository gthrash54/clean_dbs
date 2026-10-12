from pathlib import Path
from types import SimpleNamespace
import h5py
import numpy as np
from bv_fake import write_bv
from clean_dbs.read_bv import read_bv, bv_windows, ACTICHAMP_RAIL
from clean_dbs.read_tdt import _match_epoc
from clean_dbs.channel_map import channel_table, leads_implanted
from clean_dbs.datatypes import event_type
from clean_dbs.stim_events import stim_events
from clean_dbs.clean_block import clean_block

CHANNELS = [("DBS1", "µV"), ("DBS2", "µV"), ("ECOG_A1", "µV"), ("FDI_contra", "µV"), ("Mic", "µV"), ("Ch7", "µV"), ("Aux1", "BIP")]


def recording(folder, n=150001, actichamp=True):
    #7 channels of 5 types, an odd sample count (windows don't divide it), lead channel 2 railed the whole time
    rng = np.random.default_rng(4)
    data = rng.standard_normal((len(CHANNELS), n)) * 30
    data[1] = -819000.0 #flat at about twice the actiCHamp rail, like a referenced channel whose electrode and reference both rail
    vhdr = write_bv(folder, "rec", CHANNELS, data, 10000.0, [("Stimulus", "S 15", 100), ("Stimulus", "S 13", 101)])
    if actichamp: #the amplifier block real headers carry
        vhdr.write_text(vhdr.read_text(encoding="utf-8") + "\n[Comment]\nAmplifier\nactiCHamp (5001) SN: 1 Rev. 2\n", encoding="utf-8")
    return vhdr


def test_bv_windows_typed_streams_match_read_bv(tmp_path):
    vhdr = recording(tmp_path / "ug9999" / "stage2 surgery")
    whole = read_bv(vhdr)
    meta, windows = bv_windows(vhdr, window_mb=0.5) #small enough for several windows
    parts = {}
    for w in windows:
        for name, (start, x) in w.items():
            assert start == sum(p.shape[1] for p in parts.get(name, []))
            parts.setdefault(name, []).append(x)
    assert list(meta["streams"]) == list(whole["streams"]) == ["lead", "ecog", "emg", "audio", "unknown", "aux_BIP"]
    assert len(parts["lead"]) > 2
    for name, s in whole["streams"].items():
        np.testing.assert_array_equal(np.concatenate(parts[name], axis=1), s["data"])
        assert meta["streams"][name]["n_samples"] == s["data"].shape[1]
        assert {k:v for k, v in meta["streams"][name].items() if k != "n_samples"}.keys() == {k for k in s if k != "data"}
    assert meta["streams"]["lead"]["rail"] == ACTICHAMP_RAIL and np.isnan(meta["streams"]["aux_BIP"]["rail"]) #not volts: no rail


def test_bv_whole_block_rail_flagged_and_window_size_free(tmp_path):
    vhdr = recording(tmp_path / "ug9999" / "stage2 surgery")
    a = clean_block(vhdr, tmp_path / "small", window_mb=0.5)
    b = clean_block(vhdr, tmp_path / "big", window_mb=1000)
    with h5py.File(a, "r") as fa, h5py.File(b, "r") as fb:
        assert fa["flags/saturation_lead/channel"][:].tolist() == [2] #DBS2 flat past the rail for the whole block
        assert fa["streams/lead"].attrs["rail"] == ACTICHAMP_RAIL
        for name in fa["streams"]:
            np.testing.assert_array_equal(fa[f"streams/{name}/data"][:], fb[f"streams/{name}/data"][:])
        assert fa["info"].attrs["start"] == "2024-01-01T12:00:00" #vmrk clock as written, no time zone shift
        assert fa["info"].attrs["cohort"] == "ug"
    #another amplifier writing floats: range unknown, so the flat channel isn't called railed
    c = clean_block(recording(tmp_path / "other" / "ug9998" / "stage2 surgery", actichamp=False), tmp_path / "other_out")
    with h5py.File(c, "r") as f:
        assert len(f["flags/saturation_lead/channel"]) == 0 and np.isnan(f["streams/lead"].attrs["rail"])


def test_leads_implanted():
    assert leads_implanted(None) == (True, True)
    assert leads_implanted({"implant":"single", "target1":"STN", "target2":"", "notes":""}) == (True, False)
    assert leads_implanted({"implant":"dual", "target1":"STN", "target2":"GPi", "notes":"x; GPi lead never implanted: blood in cannula"}) == (True, False)
    assert leads_implanted({"implant":"dual", "target1":"STN", "target2":"GPi", "notes":"the GPi lead was dropped in"}) == (True, True)


def test_channel_table_not_implanted_and_stage2(tmp_path):
    t = channel_table(tmp_path / "uh9999" / "stage1 surgery" / "b", "ecos", "ecogdbs", 28, "STN", "", "", implanted=(True, False))
    assert t["group"][12:20] == ["lead1"] * 8 and t["target"][12] == "STN"
    assert t["group"][20:] == ["lead2_not_implanted"] * 8 and set(t["target"][20:]) == {""}
    t = channel_table(tmp_path / "uh9999" / "stage2 surgery (sedated, not GA)" / "b", "ecos", "ecogdbs", 28, "STN", "GPi", "Ad-Tech 2x6")
    assert set(t["group"][:12]) == {"ecog_anterior_not_placed", "ecog_posterior_not_placed"} #no strips in stage 2
    assert t["group"][12] == "lead1" and t["group"][20] == "lead2" and t["target"][20] == "GPi"


def test_event_types_settings_vs_pulses():
    assert [event_type(n) for n in ("PeA_", "PeB_", "PerB")] == ["stim"] * 3
    assert [event_type(n) for n in ("CnA_", "CnB_", "CntA", "AmpB", "DelC", "ChnA", "eS1p")] == ["stim_params"] * 7
    assert event_type("Cam1") == "video" and event_type("Note") == "other"


def test_stim_events_perb_and_uneven_store():
    t = np.array([1.0, 1.1, 1.2])
    epocs = {"PeA_":{"onset":t, "offset":t, "data":np.full(3, 50.0)}, "PerB":{"onset":t, "offset":t, "data":np.full(3, 25.0)},
             "CntA":{"onset":np.r_[0.0, t], "offset":t, "data":np.full(3, 2.0)}} #one onset more than values, like uh3001 CntA
    result = {"epocs":epocs, "scalars":{}, "info":{}}
    stim = stim_events(result, Path("uh9999/stage1 surgery/b"))
    assert sorted(set(stim["voice"])) == ["A", "B"] and sorted(set(stim["source"])) == ["PeA_", "PerB"]
    assert "CntA: 4 onsets but 3 values" in result["info"]["stim_warning"]


def test_match_epoc_drops_onset_not_in_tsq():
    #tdt hands CntA the 0 s onset of the PeA_ store it rides on; the .tsq knows CntA's own onsets
    e = SimpleNamespace(onset=np.array([0.0, 10.8, 10.9]), offset=np.array([10.85, 10.95]), data=np.array([2.0, 2.0]))
    headers = SimpleNamespace(stores={"CntA":SimpleNamespace(onset=np.array([10.8, 10.9]))})
    warnings = []
    onset, offset, data = _match_epoc("CntA", e, headers, warnings)
    np.testing.assert_array_equal(onset, [10.8, 10.9])
    assert len(offset) == len(data) == 2 and "CntA" in warnings[0]


def test_ug_s13_rule_needs_many_rises(tmp_path):
    #ug3003 erp001: S13 only rises, at the ERP program's spacings -> stim; a single stray S13 stays a marker
    from clean_dbs.markers import read_vmrk, label_bv
    starts = 1000 + np.cumsum(np.tile([3, 8, 20, 500], 40)) #pairs and gaps in samples at 10 kHz
    many = write_bv(tmp_path / "ug9999" / "stage2 surgery", "erp", CHANNELS[:2], np.zeros((2, 90000)), 10000.0, [("Stimulus", "S 13", int(s)) for s in starts])
    one = write_bv(tmp_path / "ug9998" / "stage1 surgery", "lfp", CHANNELS[:2], np.zeros((2, 20000)), 10000.0, [("Stimulus", "S 13", 5000)])
    r = {"epocs":read_vmrk(many)}
    stim = label_bv(r, many)["stim"]
    assert len(stim["time"]) == len(starts) and r["epocs"]["S13"]["type"] == "stim" and set(stim["confirmed"]) == {0}
    r = {"epocs":read_vmrk(one)}
    assert len(label_bv(r, one)["stim"]["time"]) == 0 and r["epocs"]["S13"]["type"] == "marker"

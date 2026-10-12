import numpy as np
from clean_dbs.markers import edges, classify_code, label_bv, read_vmrk
from clean_dbs.stim_events import stim_events
from bv_fake import write_bv, pulse_markers

FS = 100000.0 #100 kHz like most of the cohorts, 1 sample = 0.01 ms


def epocs_of(pairs):
    #{"S15": [times]} -> epocs like read_bv
    return {k:{"onset":np.asarray(v, dtype=float), "offset":np.asarray(v, dtype=float)} for k, v in pairs.items()}


def test_edges_widths():
    #S15 then S13 60 us later = line 2 high for 0.06 ms; the S13 is a falling edge, not a pulse
    ed = edges(epocs_of({"S15":[1.0, 2.0], "S13":[1.00006, 2.00006]}))
    np.testing.assert_array_equal(ed["rise"], [1, 0, 1, 0])
    np.testing.assert_allclose(ed["width_ms"][ed["rise"] == 1], [0.06, 0.06], atol=1e-6)


def test_classify_code():
    t = np.arange(50) * 0.05
    assert classify_code(np.full(50, 0.06), t) == "pulse"
    assert classify_code(np.full(50, 0.24), t) == "pulse" #sdcurve: 240 us pulse still a pulse
    assert classify_code(np.tile([0.1, 0.5, 1.2, 2.2, 3.4], 10), t) == "ttl_code"
    assert classify_code(np.full(50, 2.0), t) == "fixed_ttl" #one fixed width: stim or a task trigger, can't tell
    assert classify_code(np.full(50, 100.0), t) == "long"
    assert classify_code(np.full(50, np.nan), t) == "single"


def rec(tmp_path, cohort_folder, markers, name="rec"):
    data = np.zeros((2, 200000), dtype=np.float32)
    return write_bv(tmp_path / cohort_folder / "stage 2 surgery", name, [("DBS1", "µV"), ("DBS2", "µV")], data, FS, markers)


def test_sdcurve_end_markers_never_extra_pulses(tmp_path):
    #pulse widths 60-240 us: the end marker comes up to 0.24 ms later; the old 0.2 ms window turned 240 us ends into extra pulses
    starts = np.arange(1000, 190000, 5000)
    widths = np.tile([6, 12, 18, 24], len(starts))[:len(starts)]
    markers = []
    for s, w in zip(starts, widths):
        markers += [("Stimulus", "S 15", int(s)), ("Stimulus", "S 13", int(s + w))]
    vhdr = rec(tmp_path, "ug9999", markers)
    r = {"epocs":read_vmrk(vhdr)}
    stim = stim_events(r, vhdr)
    assert len(stim["time"]) == len(starts) and set(stim["source"]) == {"S15"}
    np.testing.assert_allclose(stim["marker_lag_ms"], widths / 100, atol=1e-6)
    assert r["epocs"]["S13"]["meaning"] == "ttl_end" and r["epocs"]["S15"]["type"] == "stim"
    assert stim["confirmed"].all() #UG S15 pulse rule is checked against a stimulator file


def test_100ms_ttl_is_not_stim(tmp_path):
    #ug baseline / mer: S15 held exactly 100 ms at a jittered ~2 Hz = the Buzzy vibration, not DBS
    rng = np.random.default_rng(0)
    starts = 1000 + np.cumsum(rng.integers(40000, 60000, 3))
    vhdr = rec(tmp_path, "ug9999", pulse_markers(starts, 10000))
    r = {"epocs":read_vmrk(vhdr)}
    assert len(stim_events(r, vhdr)["time"]) == 0
    assert r["epocs"]["S15"]["type"] == "task" and r["epocs"]["S15"]["meaning"] == "buzz_start"
    assert r["epocs"]["S15"]["kind"] == "long"


def test_unknown_cohort_fixed_ttl_stays_marker(tmp_path):
    #a 2 ms TTL every 2 s in a cohort with no rule: could be an E-Prime trigger, so never stim
    vhdr = rec(tmp_path, "zz9999", pulse_markers(np.arange(1000, 190000, 20000), 200, "R 15", "R 14"))
    r = {"epocs":read_vmrk(vhdr)}
    assert len(stim_events(r, vhdr)["time"]) == 0
    assert r["epocs"]["R15"]["type"] == "marker" and r["epocs"]["R15"]["meaning"] == "unknown"


def test_generic_pulse_is_unconfirmed_stim(tmp_path):
    vhdr = rec(tmp_path, "zz9999", pulse_markers(np.arange(1000, 190000, 2000), 6, "S 7", "S 5"))
    r = {"epocs":read_vmrk(vhdr)}
    stim = stim_events(r, vhdr)
    assert len(stim["time"]) == len(np.arange(1000, 190000, 2000)) and not stim["confirmed"].any()
    assert "generic" in r["epocs"]["S7"]["evidence"]


def test_button_and_new_segment(tmp_path):
    vhdr = rec(tmp_path, "zz9999", [("Comment", "BtnPressed", 500), ("Comment", "BtnReleased", 900)])
    r = {"epocs":read_vmrk(vhdr)}
    label_bv(r, vhdr)
    assert r["epocs"]["BtnPressed"]["meaning"] == "recorder_button_press" and r["epocs"]["NewSegment"]["type"] == "other"

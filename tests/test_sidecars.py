import h5py
import numpy as np
from clean_dbs import sidecars
from clean_dbs.clean_block import clean_block
from clean_dbs.markers import read_vmrk
from clean_dbs.stim_events import stim_events
from bv_fake import write_bv, write_mcs

FS = 100000.0


def program(n=40, seed=1):
    #MC011-style program: biphasic pulse per row with a jittered gap, plus a sync TTL whose width codes the condition
    rng = np.random.default_rng(seed)
    amps = rng.choice([0, 500, 1000, 1500], n) * rng.choice([-1, 1], n)
    gaps = rng.integers(70000, 120000, n)
    ttl = (np.abs(amps) // 500 + 1) * 100 #100-400 us
    ch1 = [(a, 60, -a, 60, 0, g, 0) for a, g in zip(amps, gaps)]
    sync = [(1, w, 0, 0, 0, 120 + g - w, 0) for w, g in zip(ttl, gaps)]
    return ch1, sync, amps, ttl, gaps


def test_read_mcs_and_table(tmp_path):
    ch1, sync, amps, ttl, gaps = program()
    m = sidecars.read_mcs(write_mcs(tmp_path / "results.txt", ch1, sync))
    assert m["mode"] == "current" and m["units"] == "uA" and m["n_channels"] == 8
    t = sidecars.mcs_table(m)
    assert len(t["onset"]) == 40 #one row per sync TTL, 0 uA rows kept (sham conditions still send their TTL)
    np.testing.assert_allclose(t["amplitude"], amps)
    np.testing.assert_allclose(t["ttl_ms"], ttl / 1000)
    np.testing.assert_allclose(np.diff(t["onset"]), (120 + gaps[:-1]) / 1e6)


def test_repeat_column():
    m = {"mode":"voltage", "units":"mV", "n_channels":1, "format":3,
         "channels":{1:np.array([[3000, 60, 0, 0, 0, 99940, 3]], dtype=float)}}
    e = sidecars.mcs_events(m)[1]
    np.testing.assert_allclose(e["onset"], [0, 0.1, 0.2]) #the row is played 3 times


def test_read_ttl_txt(tmp_path):
    p = tmp_path / "TTL.txt"
    p.write_text("voltagesp={'ch1_0V_0.06ms' 'ch1_0.5V_0.06ms' };\npulesp_ms=[0.1 0.2 ];\nvoltagesn={'-ch1_0V_0.06ms' };\npulesn_ms=[0.3 ];\n")
    t = sidecars.read_ttl_txt(p)
    assert list(t["condition"]) == ["ch1_0V_0.06ms", "ch1_0.5V_0.06ms", "-ch1_0V_0.06ms"]
    np.testing.assert_allclose(t["ttl_ms"], [0.1, 0.2, 0.3])


def test_match_times_offset_missing_drift():
    rng = np.random.default_rng(3)
    events = np.cumsum(rng.uniform(0.05, 0.15, 300))
    markers = 12.0 + events[20:250] * (1 + 5e-5) #recording starts late, stops early, clocks differ by 0.005%
    markers = np.delete(markers, [50, 51, 120]) #a few markers missing
    m = sidecars.match_times(markers, events)
    assert m["matched"] == 1 and m["n_pairs"] == len(markers)
    assert m["pairs"][0, 1] == 20 and m["max_err_ms"] < 0.1
    assert sidecars.match_times(markers, rng.uniform(0, 30, 300))["matched"] == 0 #random times don't match


def test_eprime(tmp_path):
    text = ("*** Header Start ***\r\nSubject: 91\r\nSessionDate: 01-01-2024\r\n*** Header End ***\r\n"
            "\t*** LogFrame Start ***\r\n\tProcedure: TrialProc\r\n\tCondition: GO\r\n\tStimulusDisplay.RT: 433\r\n\tSubject: 91\r\n\t*** LogFrame End ***\r\n"
            "\t*** LogFrame Start ***\r\n\tProcedure: TrialProc\r\n\tCondition: NOGO\r\n\tStimulusDisplay.RT: 0\r\n\t*** LogFrame End ***\r\n")
    p = tmp_path / "Go_NoGo_v5-1-1.txt"
    p.write_bytes(text.encode("utf-16"))
    kind, t = sidecars.read_sidecar(p)
    assert kind == "eprime" and list(t["Condition"]) == ["GO", "NOGO"] and "Subject" not in t


def test_sidecar_lines_up_and_is_saved(tmp_path):
    #a recording whose S15 markers are the sync outputs of results.txt next to it: amplitudes and the condition get copied per pulse
    ch1, sync, amps, ttl, gaps = program()
    folder = tmp_path / "raw" / "mc9999" / "stage 2 surgery"
    folder.mkdir(parents=True)
    write_mcs(folder / "results.txt", ch1, sync)
    (folder / "TTL.txt").write_text("cond={'a' 'b' 'c' 'd' };\ncond_ms=[0.1 0.2 0.3 0.4 ];\n")
    onsets = 5000 + np.concatenate(([0], np.cumsum(120 + gaps[:-1]))) / 10 #us -> samples at 100 kHz
    markers = []
    for s, w in zip(onsets[3:], ttl[3:]): #the recording started after the first 3 trials
        markers += [("Stimulus", "S 15", int(round(s))), ("Stimulus", "S 7", int(round(s + w / 10)))]
    vhdr = write_bv(folder, "0-1", [("DBS1", "µV"), ("17", "µV")], np.zeros((2, int(onsets[-1] + 20000)), dtype=np.float32), FS, markers)
    stim = stim_events({"epocs":read_vmrk(vhdr)}, vhdr)
    assert len(stim["time"]) == 37 and stim["confirmed"].all()
    np.testing.assert_allclose(stim["amplitude_ua"], amps[3:])
    assert list(stim["condition"]) == [["a", "b", "c", "d"][int(w // 100) - 1] for w in ttl[3:]]
    assert stim["marker_lag_meaning"][0] == "ttl_code"
    out = clean_block(vhdr, tmp_path / "out")
    with h5py.File(out, "r") as f:
        assert f["sidecars/results"].attrs["matched"] == 1 and f["sidecars/results"].attrs["kind"] == "mcs"
        assert f["sidecars/TTL"].attrs["kind"] == "ttl_txt"
        np.testing.assert_allclose(f["stim/amplitude_ua"][:], amps[3:])
        assert f["epocs/S15"].attrs["type"] == "stim" and f["epocs/S7"].attrs["meaning"] == "ttl_end"
        assert f["streams/lead"].attrs["type"] == "lfp" and f["streams/unknown"].attrs["type"] == "unknown"
        assert f["streams/lead/channel_group"][0].decode() == "lead1"
        assert "impedance_brainvision" not in f["flags"] #no impedance block in this fake header

#fixes from the trial run on real recordings of the non-uh/UG cohorts (2026-10-12); paths are made up in the shape Box uses
from pathlib import Path
from clean_dbs.cohorts import condition_of
from clean_dbs.block_info import describe_block, output_path

BOX = Path("box") / "100.0000001"


def test_folder_label_when_the_name_says_nothing():
    #the docstring's rule: a folder wins when the name is blank, generic folder labels included
    assert describe_block(BOX / "SENSORY002" / "Initial Peripheral Stim" / "Bicep.vhdr")["condition"] == "stim"
    assert describe_block(BOX / "PNS002" / "Stim" / "PNS1.vhdr")["condition"] == "stim"
    assert describe_block(BOX / "BI020" / "LFP" / "2.vhdr")["condition"] == "lfp"
    #a name that says something still beats a generic folder
    assert describe_block(BOX / "BI029" / "LFP" / "MER.vhdr")["condition"] == "mer"
    assert describe_block(BOX / "DLPFC045" / "Misc" / "x1.vhdr")["condition"] == ""


def test_study_id_in_the_name_is_not_a_condition():
    #BPP032.vhdr names the person; "BPP" isn't a paired pulse recording by itself
    assert describe_block(BOX / "BPP032" / "BPP032.vhdr")["condition"] == ""
    assert describe_block(BOX / "BPP032" / "BPP032_stimRight.vhdr")["condition"] == "stim"
    assert describe_block(BOX / "BPP012.vhdr")["condition"] == "" #straight in a number folder
    assert describe_block(BOX / "BI020" / "LFP" / "BI020.vhdr")["condition"] == "lfp"
    #run numbers that aren't IDs stay
    assert describe_block(BOX / "PP012" / "PP1.vhdr")["condition"] == "stim_paired_pulse"
    assert describe_block(BOX / "DLPFC045" / "pairedpulsesDLPFC045.vhdr")["condition"] == "stim_paired_pulse"


def test_same_id_nested_in_itself():
    path = BOX / "DLPFC110" / "Contents of DLPFC110 folder" / "DLPFC110" / "practice.vhdr"
    info = describe_block(path)
    assert info["patient"] == "DLPFC110" and info["parent_id"] == ""
    assert info["session"] == "contents_of_dlpfc110_folder/dlpfc110"
    #so it can't land on the same file as a practice.vhdr at the ID root
    assert output_path(path, "out") != output_path(BOX / "DLPFC110" / "practice.vhdr", "out")
    #a different ID inside still wins, with the outer one as parent
    info = describe_block(BOX / "PP012" / "PP012" / "CHR003" / "x.vhdr")
    assert info["patient"] == "CHR003" and info["parent_id"] == "PP012"


def test_negated_task_and_prestim():
    assert condition_of("itbs no task")[0] == "stim_itbs"
    assert condition_of("itbs_notask")[0] == "stim_itbs"
    assert condition_of("task1_itbs")[0] == "task"
    assert condition_of("restingstateprestimulation")[0] == "rest"
    assert condition_of("corticalstim")[0] == "stim"


def test_impedance_checks_are_not_the_task():
    info = describe_block(BOX / "DLPFC045" / "Task data" / "AXCPT" / "ieeg" / "impedencetest_preopening impedencetab.vhdr")
    assert info["condition"] == "impedance_check"
    assert info["state"] == "" #"preopening" isn't preop
    assert condition_of("Preop/practice1")[1] == "preop"


def test_typos_and_named_settings():
    assert condition_of("12 month/comtact1")[0] == "contact"
    assert condition_of("best_psuedo_27_60_130")[0] == "stim_setting"
    assert condition_of("1mA_DLPFCatMEPthrehold_postDBS_1_2")[0] == "stim_mep_threshold"
    assert condition_of("joystick")[0] == "task_joystick"
    assert condition_of("preop/61offmeds")[1] == "preop;off_meds"
    assert condition_of("preop/61onmeds")[1] == "preop;on_meds"


#--- sidecars, markers and reading, from the same trial run ---
import numpy as np
from clean_dbs import sidecars
from clean_dbs.markers import label_bv
from clean_dbs.read_bv import read_bv
from bv_fake import write_bv, write_mcs


def test_ttl_list_written_as_rows(tmp_path):
    #BI016 STIM/TTL.txt writes both lists as MATLAB rows separated by ";"
    p = tmp_path / "TTL.txt"
    p.write_text("voltagesp={'a' 'b' ; \r\n'c' 'd' ; \r\n};\r\npulesp_ms=[0.2 0.4 ; \r\n0.6 0.8 ; \r\n];\r\n")
    t = sidecars.read_ttl_txt(p)
    assert list(t["condition"]) == ["a", "b", "c", "d"]
    np.testing.assert_allclose(t["ttl_ms"], [0.2, 0.4, 0.6, 0.8])


def test_program_found_by_its_header_not_its_name(tmp_path):
    #CHR and PT name their MC_Stimulus programs after the setting ("0-1_80us_M.txt", "pulse_directionality_ch1.txt")
    write_mcs(tmp_path / "0-1_80us_M.txt", [(100, 60, -100, 60, 0, 1000, 0)], [(1, 100, 0, 0, 0, 1020, 0)])
    (tmp_path / "notes about the day.txt").write_text("not a program")
    (tmp_path / "x.vhdr").write_text("")
    names = [f.name for f in sidecars.find_sidecars(tmp_path / "x.vhdr")]
    assert names == ["0-1_80us_M.txt"]


def _ttl_stim(widths_ms):
    n = len(widths_ms)
    return {"time":np.arange(n, dtype=float), "marker_lag_ms":np.asarray(widths_ms, dtype=float), "source":np.full(n, "S7")}


def test_two_ttl_tables_that_disagree_leave_the_condition_blank(tmp_path):
    (tmp_path / "x.vhdr").write_text("")
    (tmp_path / "A_TTL.txt").write_text("c={'a1' 'a2' };\nc_ms=[0.1 0.2 ];\n")
    (tmp_path / "B_TTL.txt").write_text("c={'b1' 'b2' };\nc_ms=[0.1 0.2 ];\n")
    stim = _ttl_stim([0.1, 0.2, 0.1, 0.2])
    out = sidecars.attach_sidecars(stim, tmp_path / "x.vhdr")
    assert "condition" not in stim
    assert out["A_TTL"]["matched"] == 0 and out["A_TTL"]["match"]["ambiguous"] == 1
    #tables that say the same thing are fine
    (tmp_path / "B_TTL.txt").write_text("c={'a1' 'a2' };\nc_ms=[0.1 0.2 ];\n")
    stim = _ttl_stim([0.1, 0.2, 0.1, 0.2])
    sidecars.attach_sidecars(stim, tmp_path / "x.vhdr")
    assert list(stim["condition"]) == ["a1", "a2", "a1", "a2"]


def test_ttl_table_named_after_the_matched_program_wins(tmp_path):
    rng = np.random.default_rng(5)
    gaps = rng.integers(70000, 120000, 30) #us
    widths = rng.choice([100, 200], 30) #us
    write_mcs(tmp_path / "P_80us.txt", [(500, 80, -500, 80, 0, g, 0) for g in gaps], [(1, w, 0, 0, 0, 160 + g - w, 0) for w, g in zip(widths, gaps)])
    (tmp_path / "P_80us_TTL.txt").write_text("c={'p1' 'p2' };\nc_ms=[0.1 0.2 ];\n")
    (tmp_path / "Q_20us_TTL.txt").write_text("c={'q1' 'q2' };\nc_ms=[0.1 0.2 ];\n")
    (tmp_path / "x.vhdr").write_text("")
    times = 3.0 + np.concatenate(([0], np.cumsum(160 + gaps[:-1]))) / 1e6
    stim = {"time":times, "marker_lag_ms":widths / 1000, "source":np.full(30, "S7")}
    out = sidecars.attach_sidecars(stim, tmp_path / "x.vhdr")
    assert out["P_80us"]["matched"] == 1
    assert list(stim["condition"]) == [{100:"p1", 200:"p2"}[w] for w in widths]
    assert out["Q_20us_TTL"]["matched"] == 0 and out["P_80us_TTL"]["matched"] == 1


def test_a_few_chance_pairs_are_not_a_match():
    rng = np.random.default_rng(7)
    events = np.cumsum(rng.uniform(0.05, 0.15, 300))
    markers = np.concatenate((events[:4] + 2.0, 100 + np.cumsum(rng.uniform(0.003, 0.004, 500))))
    m = sidecars.match_times(markers, events)
    assert m["n_pairs"] >= 4 and m["matched"] == 0


def test_program_that_follows_one_code_only(tmp_path):
    #TM002 Burst1: the program's 20 sync outputs are the S15 burst onsets, the S14 pulses inside each burst aren't in it
    rng = np.random.default_rng(9)
    gaps = rng.integers(200000, 300000, 20) #us between bursts
    write_mcs(tmp_path / "results.txt", [(1000, 60, -1000, 60, 0, g, 0) for g in gaps], [(1, 100, 0, 0, 0, 20 + g, 0) for g in gaps])
    (tmp_path / "x.vhdr").write_text("")
    onsets = 1.0 + np.concatenate(([0], np.cumsum(120 + gaps[:-1]))) / 1e6
    inside = (onsets[:, None] + np.arange(1, 9) * 0.008).ravel()
    t = np.concatenate((onsets, inside))
    src = np.array(["S15"] * len(onsets) + ["S14"] * len(inside))
    order = np.argsort(t)
    stim = {"time":t[order], "source":src[order], "marker_lag_ms":np.full(len(t), 0.06)}
    out = sidecars.attach_sidecars(stim, tmp_path / "x.vhdr")
    assert out["results"]["matched"] == 1 and out["results"]["match"]["source"] == "S15"
    burst = stim["source"] == "S15"
    np.testing.assert_allclose(stim["amplitude_ua"][burst], 1000)
    assert np.isnan(stim["amplitude_ua"][~burst]).all() and (stim["sidecar"][burst] == "results.txt").all()


def test_falling_edge_names_the_code_on_its_own_port(tmp_path):
    #R15 -> R14 is one line going low on the R port, even when an S marker comes in between
    epocs = {"R15":{"onset":np.array([1.0, 2.0, 3.0]), "offset":np.array([1.0, 2.0, 3.0])},
             "S3":{"onset":np.array([1.05, 2.05, 3.05]), "offset":np.array([1.05, 2.05, 3.05])},
             "R14":{"onset":np.array([1.1, 2.1, 3.1]), "offset":np.array([1.1, 2.1, 3.1])}}
    label_bv({"epocs":epocs}, tmp_path / "none.vhdr", cohort="")
    assert epocs["R14"]["meaning"] == "ttl_end" and "raised by R15" in epocs["R14"]["evidence"]


def test_header_with_a_repeated_channel_name(tmp_path):
    #BI002 12 Month Baseline names Ch20 and Ch25 both LFP1; MNE can't open that on its own
    data = np.vstack([np.full(500, 1.0), np.full(500, 2.0), np.full(500, 3.0)]).astype(np.float32)
    vhdr = write_bv(tmp_path / "BI999", "Baseline", [("LFP1", "µV"), ("Ch2", "µV"), ("LFP1", "µV")], data, 1000.0, [("Stimulus", "S 1", 100)])
    r = read_bv(vhdr)
    lead = r["streams"]["lead"]
    assert list(lead["channels"]) == [1, 3] and lead["channel_names"] == ["LFP1", "LFP1"]
    np.testing.assert_allclose(lead["data"][:, 0], [1e-6, 3e-6])
    assert all("2 channels" in s for s in lead["channel_source"])
    assert "S1" in r["epocs"] #markers still read through the copied header


def test_binary_mc_stimulus_file_is_listed_not_dropped(tmp_path):
    #MC012 results.txt is MC_Stimulus's own binary save; the README says every sidecar found is kept, so it shows up as unreadable
    p = tmp_path / "results.txt"
    p.write_bytes(b"\x1eMC_Stimulus II Version 3.4.4.00\x8d\x03\x02\x00\x00\x00CGXGridParam\x00\x01")
    kind, table = sidecars.read_sidecar(p)
    assert kind == "mcs_binary" and "ASCII" in table["note"][0]
    (tmp_path / "x.vhdr").write_text("")
    out = sidecars.attach_sidecars({"time":np.zeros(0)}, tmp_path / "x.vhdr")
    assert out["results"]["kind"] == "mcs_binary" and out["results"]["matched"] == 0


def test_code_carrying_a_few_of_the_programs_syncs(tmp_path):
    #DLPFC047 itbs: 40 S15 among 1,103 S7 rises are all syncs of results_itbs.dat, but 40 sparse rises can't line up on their own
    from clean_dbs.markers import read_vmrk
    rng = np.random.default_rng(11)
    gaps = rng.integers(40000, 60000, 600) #us
    folder = tmp_path / "DLPFC999"
    folder.mkdir()
    write_mcs(folder / "results_itbs.dat", [(800, 60, -800, 60, 0, g, 0) for g in gaps], [(1, 2000, 0, 0, 0, 120 + g - 2000, 0) for g in gaps])
    onsets = 1000 + np.round(np.concatenate(([0], np.cumsum(120 + gaps[:-1]))) / 100).astype(int) #samples at 10 kHz
    markers = []
    for i, s in enumerate(onsets):
        markers += [("Stimulus", "S 15" if i % 100 == 50 else "S  7", int(s)), ("Stimulus", "S  3", int(s + 20))]
    vhdr = write_bv(folder, "itbs", [("DBS1", "µV")], np.zeros((1, int(onsets[-1] + 2000)), dtype=np.float32), 10000.0, markers)
    epocs = read_vmrk(vhdr)
    lab = label_bv({"epocs":epocs}, vhdr)
    assert epocs["S7"]["type"] == "stim" and epocs["S15"]["type"] == "stim" and "together with S7" in epocs["S15"]["evidence"]
    assert len(lab["stim"]["time"]) == 600 and (lab["stim"]["sidecar"] == "results_itbs.dat").all()

import warnings
import h5py
import numpy as np
import pandas as pd
from bv_fake import write_bv
from clean_dbs.clean_block import clean_block
from clean_dbs.coords import _load
from clean_dbs.datatypes import load_type

#fixes from the 2026-10-12 review: missing contacts.csv (M6), coordinates on non-lead streams (L6), load_type labels (L7)
CHANNELS = [("DBS1", "µV"), ("DBS2", "µV"), ("ECOG_A1", "µV"), ("FDI_contra", "µV")]


def recording(folder):
    #a made-up brainvision recording with lead, ecog and emg channels, 2 s at 1 kHz
    data = np.random.default_rng(5).standard_normal((len(CHANNELS), 2000)) * 30
    return write_bv(folder, "rec", CHANNELS, data, 1000.0, [])


def contacts_table(out):
    #lead 1 of the fake patient at its final position only, x = contact number
    (out / "_cases").mkdir(parents=True)
    rows = [{"study_id":"uh9999", "lead":"1", "contact":str(c), "depth_index":"", "x":str(c), "y":"0", "z":"0", "space":"ACPC test"}
            for c in range(1, 9)]
    pd.DataFrame(rows).to_csv(out / "_cases" / "contacts.csv", index=False)
    _load.cache_clear()


def test_missing_contacts_csv_warns_once_and_is_kept_in_the_file(tmp_path):
    raw = tmp_path / "raw"
    blocks = [recording(raw / "uh9999" / "stage1 surgery"), recording(raw / "uh9999" / "stage2 surgery")]
    ug = recording(raw / "ug9999" / "stage2 surgery")
    _load.cache_clear()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        files = [clean_block(b, tmp_path / "out") for b in blocks]
        clean_block(ug, tmp_path / "out_ug") #UG has no recons, so no contacts.csv is fine there
    said = [w for w in caught if "contacts.csv" in str(w.message)]
    assert len(said) == 1 #once per run, not once per block
    for path in files:
        with h5py.File(path, "r") as f:
            assert f["info"].attrs["coord_status"] == "no_contacts_csv" #tells it apart from a patient with no recon
            assert "channel_x" not in f["streams/lead"]


def test_only_lead_streams_get_coordinates(tmp_path):
    vhdr = recording(tmp_path / "raw" / "uh9999" / "stage2 surgery")
    contacts_table(tmp_path / "out")
    with h5py.File(clean_block(vhdr, tmp_path / "out"), "r") as f:
        assert f["info"].attrs["coord_status"] == "recon"
        assert f["streams/lead/channel_x"][:].tolist() == [1.0, 2.0] and f["streams/lead"].attrs["coord_space"] == "ACPC test"
        for name in ("ecog", "emg"): #these used to get all-NaN channel_x/y/z and a coord_space
            assert "channel_x" not in f[f"streams/{name}"] and "coord_space" not in f[f"streams/{name}"].attrs


def test_load_type_gives_channel_labels_and_coordinates(tmp_path):
    vhdr = recording(tmp_path / "raw" / "uh9999" / "stage2 surgery")
    contacts_table(tmp_path / "out")
    lfp = load_type(clean_block(vhdr, tmp_path / "out"), "lfp")
    lead, ecog = lfp["lead"], lfp["ecog"]
    assert lead["data"].shape == (2, 2000) and lead["channels"].tolist() == [1, 2] and lead["fs"] == 1000.0 #what it gave before
    assert lead["channel_names"] == ["dbs1", "dbs2"] and lead["channel_group"] == ["lead1", "lead1"] and lead["channel_contact"] == ["1", "2"]
    assert lead["channel_x"].dtype == np.float64 and lead["channel_x"].tolist() == [1.0, 2.0]
    assert lead["channel_position"] == ["final", "final"] and lead["coord_space"] == "ACPC test"
    assert ecog["channel_names"] == ["ecog_a1"] and "channel_x" not in ecog
    assert set(load_type(clean_block(vhdr, tmp_path / "out"), "emg")) == {"emg"}

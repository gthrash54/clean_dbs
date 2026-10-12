import h5py
import numpy as np
import pytest
from clean_dbs.save_h5 import save_h5


def fake_result():
    """
    Input: nothing
    Returns: a small made-up reader result (streams, epocs, scalars, info, flags) shaped like read_tdt's output
    """
    return {
        "streams":{"ecos":{"data":np.arange(12, dtype=np.float32).reshape(3, 4), "fs":48828.125, "start_time":0.25,
                           "channels":np.array([1, 2, 3]), "device":"ecogdbs", "units":"V", "type":"lfp", "subtype":"ecog+lead"}},
        "epocs":{"Cam1":{"onset":np.array([0.1, 0.2]), "offset":np.array([0.15, 0.25]), "data":np.array([1.0, 2.0]), "type":"video"}},
        "scalars":{"eS1p":{"ts":np.array([0.5, 0.6, 0.7]), "data":np.arange(6.0).reshape(2, 3), "type":"stim_params"}},
        "info":{"block":"block1", "duration":12.5, "has_video":True},
        "flags":{"saturation_ecos":{"channel":np.array([2]), "start":np.array([1.0]), "end":np.array([1.5])}, "dnu":np.int8(0)},
    }


def make_source(tmp_path):
    #a fake raw block folder, deep enough that provenance keeps only the last three parts
    source = tmp_path / "raw" / "uh9999" / "stage1 surgery" / "block1"
    source.mkdir(parents=True)
    (source / "block1.tev").write_bytes(b"x" * 100)
    (source / "block1.tsq").write_bytes(b"x" * 20)
    return source


def test_round_trip(tmp_path):
    out = tmp_path / "out" / "uh9999" / "stage1" / "block1.h5"
    assert save_h5(fake_result(), out, make_source(tmp_path)) == out
    with h5py.File(out, "r") as f:
        g = f["streams/ecos"]
        assert g["data"].dtype == np.float32
        np.testing.assert_array_equal(g["data"][:], np.arange(12).reshape(3, 4)) #channels x samples, not flipped
        np.testing.assert_array_equal(g["channels"][:], [1, 2, 3])
        assert g.attrs["fs"] == 48828.125 and g.attrs["start_time"] == 0.25
        assert (g.attrs["device"], g.attrs["units"], g.attrs["type"], g.attrs["subtype"]) == ("ecogdbs", "V", "lfp", "ecog+lead")
        assert f["epocs/Cam1"].attrs["type"] == "video" #text goes on as an attribute, arrays as datasets
        np.testing.assert_array_equal(f["epocs/Cam1/data"][:], [1.0, 2.0])
        np.testing.assert_array_equal(f["scalars/eS1p/data"][:], np.arange(6.0).reshape(2, 3))
        assert f["scalars/eS1p"].attrs["type"] == "stim_params"
        assert f["info"].attrs["block"] == "block1" and f["info"].attrs["duration"] == 12.5
        assert f["info"].attrs["has_video"] == 1 #bools are stored as 0/1 so MATLAB doesn't get the text 'TRUE'
        np.testing.assert_array_equal(f["flags/saturation_ecos/start"][:], [1.0])
        assert f["flags"].attrs["dnu"] == 0
        #provenance: source from the study ID down, never the full path on this computer
        assert f.attrs["source"] == "uh9999/stage1 surgery/block1"
        assert f.attrs["source_bytes"] == 120
        assert "created" in f.attrs and "version_h5py" in f.attrs and "clean_dbs_commit" in f.attrs


def test_part_file_renamed(tmp_path):
    out = tmp_path / "block1.h5"
    save_h5(fake_result(), out, make_source(tmp_path))
    assert out.is_file()
    assert not (tmp_path / "block1.h5.part").exists() #renamed once finished


def test_crash_leaves_only_part(tmp_path):
    result = fake_result()
    del result["info"] #makes save_h5 fail partway through writing
    out = tmp_path / "block1.h5"
    with pytest.raises(KeyError):
        save_h5(result, out, make_source(tmp_path))
    assert not out.exists() #a crashed job never leaves a finished-looking .h5
    assert (tmp_path / "block1.h5.part").exists()


def test_empty_stream(tmp_path):
    #a stream with no samples (e.g. a short window) shouldn't break the compressed write
    result = fake_result()
    result["streams"]["ecos"]["data"] = np.zeros((3, 0), dtype=np.float32)
    out = save_h5(result, tmp_path / "empty.h5", make_source(tmp_path))
    with h5py.File(out, "r") as f:
        assert f["streams/ecos/data"].shape == (3, 0)

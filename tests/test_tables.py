import pandas as pd
from clean_dbs.tables import block_rows
from clean_dbs.coords import _leads_in_place, _load


def _blocks(stage_col=True):
    rows = [{"study_id":"ug3007", "block":"erps01", "stage":"stage1", "leads_in_place":"1"},
            {"study_id":"ug3007", "block":"erps01", "stage":"stage2", "leads_in_place":"1"},
            {"study_id":"uh3021", "block":"increment1_depth1", "stage":"", "leads_in_place":"0"}]
    t = pd.DataFrame(rows)
    return t if stage_col else t.drop(columns="stage")


def test_same_block_name_in_two_stages_stays_apart():
    t = _blocks()
    assert len(block_rows(t, "UG3007", "erps01", "stage2")) == 1 #case-insensitive ID, one row per stage
    assert len(block_rows(t, "ug3007", "erps01")) == 2 #no stage given: ambiguous, callers leave it out
    assert len(block_rows(t, "uh3021", "increment1_depth1", "stage1")) == 1 #a blank stage matches any stage


def test_old_tables_without_a_stage_column_still_work():
    assert len(block_rows(_blocks(stage_col=False), "ug3007", "erps01", "stage1")) == 2


def test_leads_in_place_uses_the_stage(tmp_path):
    (tmp_path / "_cases").mkdir()
    t = _blocks()
    t.loc[1, "leads_in_place"] = "2"
    t.to_csv(tmp_path / "_cases" / "blocks.csv", index=False)
    _load.cache_clear()
    assert _leads_in_place(tmp_path, "ug3007", "erps01", "stage1") == 1
    assert _leads_in_place(tmp_path, "ug3007", "erps01", "stage2") == 2
    assert _leads_in_place(tmp_path, "ug3007", "erps01") == 0 #ambiguous without a stage: unknown, not a guess

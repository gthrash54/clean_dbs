import numpy as np
import pandas as pd
from clean_dbs.coords import contact_coords, _load


def _tables(root):
    #a tiny contacts.csv (lead 1 at depths 1-5, lead 2 final only) and blocks.csv for one fake patient
    cases = root / "_cases"
    cases.mkdir(parents=True)
    rows = []
    for depth in ("1", "2", "3", "4", "5"):
        for c in range(1, 9):
            rows.append({"study_id":"uh9999", "lead":"1", "contact":str(c), "depth_index":depth, "x":str(c), "y":depth, "z":"0", "space":"ACPC test"})
    for c in range(1, 9):
        rows.append({"study_id":"uh9999", "lead":"2", "contact":str(c), "depth_index":"", "x":str(-c), "y":"9", "z":"0", "space":"ACPC test"})
    pd.DataFrame(rows).to_csv(cases / "contacts.csv", index=False)
    pd.DataFrame([{"study_id":"uh9999", "block":"clinical_dbs1", "leads_in_place":"1"},
                  {"study_id":"uh9999", "block":"crosstalk", "leads_in_place":"2"},
                  {"study_id":"uh9999", "block":"mer1_traj1_v1", "leads_in_place":"0"}]).to_csv(cases / "blocks.csv", index=False)
    _load.cache_clear()


GROUPS = ["ecog_anterior"] * 2 + ["lead1"] * 8 + ["lead2"] * 8
CONTACTS = ["1", "2"] + [str(c) for c in range(1, 9)] * 2


def test_increment_depth_and_earlier_lead(tmp_path):
    _tables(tmp_path)
    x, y, z, pos, space = contact_coords(tmp_path, "uh9999", "increment1_depth3", GROUPS, CONTACTS)
    assert space == "ACPC test"
    assert np.isnan(x[0]) and pos[0] == "" #ecog channels get no coordinates
    assert pos[2] == "depth3" and y[2] == 3.0 and x[2] == 1.0 #lead 1 contact 1 at depth 3
    assert all(p == "" for p in pos[10:]) #lead 2 not in the brain yet during increment1
    x, y, z, pos, space = contact_coords(tmp_path, "uh9999", "increment2_depth2", GROUPS, CONTACTS)
    assert pos[2] == "final" and y[2] == 5.0 #lead 1 already at its final (imaged) position
    assert pos[10] == "" #lead 2 has no depth 2 row in the table, so nothing


def test_leads_in_place_from_blocks_table(tmp_path):
    _tables(tmp_path)
    _, _, _, pos, _ = contact_coords(tmp_path, "UH9999", "clinical_dbs1", GROUPS, CONTACTS) #case-insensitive id
    assert pos[2] == "final" and all(p == "" for p in pos[10:])
    _, y, _, pos, _ = contact_coords(tmp_path, "uh9999", "crosstalk", GROUPS, CONTACTS)
    assert pos[10] == "final" and y[10] == 9.0 #lead 2 final from its blank-depth rows
    _, _, _, pos, _ = contact_coords(tmp_path, "uh9999", "mer1_traj1_v1", GROUPS, CONTACTS)
    assert all(p == "" for p in pos) #no lead in the brain during MER
    _, _, _, pos, _ = contact_coords(tmp_path, "uh9999", "STN_dbs1", GROUPS, CONTACTS, stage="stage2")
    assert pos[2] == "final" and pos[10] == "final" #stage 2: everything implanted at its final position


def test_no_tables_or_patient(tmp_path):
    _load.cache_clear()
    x, _, _, pos, space = contact_coords(tmp_path, "uh9999", "crosstalk", GROUPS, CONTACTS)
    assert space == "" and all(np.isnan(v) for v in x)
    _tables(tmp_path)
    _, _, _, pos, space = contact_coords(tmp_path, "uh0000", "crosstalk", GROUPS, CONTACTS)
    assert space == "" and all(p == "" for p in pos)


def test_depth_name_without_underscore():
    #uh3002 names its inching blocks increment1depth3; they used to fall through to the final position
    from clean_dbs.coords import _lead_positions
    assert _lead_positions("unused", "uh3002", "increment1depth3", "stage1") == {1:"depth3"}
    assert _lead_positions("unused", "uh3002", "Increment2_Depth4", "stage1") == {1:"final", 2:"depth4"}

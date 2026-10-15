import re
from functools import lru_cache
from pathlib import Path
import numpy as np
import pandas as pd
from .tables import block_rows

#contact coordinates come from <cleaned root>/_cases/contacts.csv, built from the lab's Brainlab lead reconstructions
#(UH3 data/uh3 stl recons: object_centers_acpc.csv). Space is Brainlab ACPC mm in the patient's own hemisphere (+x right,
#+y anterior, +z superior). Only the final lead position is imaged (the lab calls it depth5); the other inching depths are that
#position slid along the lead's trajectory by the depth table's steps. Only uh cases have recons so far (18 patients).

#cohorts with recons in contacts.csv. Cleaning one of them without the table warns (clean_block), because every file would come out
#with no coordinates and look just like a patient who has no recon
RECON_COHORTS = ("uh",)
#a channel_group that is a lead in the brain (lead1, lead2); lead2_not_implanted and lead2_not_connected aren't
LEAD_GROUP = re.compile(r"lead(\d)")


def contact_coords(out_root, study_id, block, groups, contacts, stage=""):
    """
    Where each channel's lead contact was in the brain during this recording.
    Input:
    out_root: top of the cleaned set (holds _cases/contacts.csv and _cases/blocks.csv)
    study_id, block: the recording (e.g. "uh3021", "increment1_depth3")
    groups, contacts: per-channel channel_group / channel_contact from channel_map (lead1/lead2 and "1".."8"; others are skipped)
    stage: "stage1", "stage2", ... (stage 2 recordings have every implanted lead at its final position)
    Returns: (x, y, z, position, space) where x/y/z are float lists (NaN = not a lead contact, lead not in the brain yet, or no recon),
             position is a per-channel text list ("depth3", "final", or ""), and space describes the coordinate space.
             space is "" when there are no coordinates to write: no contacts.csv, no recon for this patient, or no lead channel in this
             stream (an ecog or emg stream gets no coordinates at all, not a column of NaN)
    """
    table = _load(str(Path(out_root) / "_cases" / "contacts.csv"))
    n = len(groups)
    x, y, z, pos = [np.nan] * n, [np.nan] * n, [np.nan] * n, [""] * n
    if table is None or not any(LEAD_GROUP.fullmatch(str(g)) for g in groups):
        return x, y, z, pos, ""
    rows = table[table["study_id"].str.lower() == study_id.lower()]
    if rows.empty:
        return x, y, z, pos, ""
    in_brain = _lead_positions(out_root, study_id, block, stage)
    for i, (group, contact) in enumerate(zip(groups, contacts)):
        m = LEAD_GROUP.fullmatch(str(group))
        if not m or not str(contact).isdigit():
            continue
        where = in_brain.get(int(m.group(1)), "")
        if not where:
            continue
        lead_rows = rows[(rows["lead"] == m.group(1)) & (rows["contact"] == str(contact))]
        if where == "final": #the imaged position; patients without inching depths only have a blank depth
            pick = lead_rows[lead_rows["depth_index"] == "5"]
            pick = pick if len(pick) else lead_rows[lead_rows["depth_index"] == ""]
        else:
            pick = lead_rows[lead_rows["depth_index"] == where.replace("depth", "")]
        if len(pick):
            x[i], y[i], z[i] = (float(pick.iloc[0][k]) for k in ("x", "y", "z"))
            pos[i] = where
    space = rows["space"].iloc[0] if "space" in rows else ""
    return x, y, z, pos, space


def coord_status(out_root, study_id):
    """
    Whether a file could get contact coordinates, saved as info coord_status so a missing table can't pass for a patient with no recon.
    Input: out_root (top of the cleaned set), study_id (e.g. "uh3021")
    Returns: "recon" (_cases/contacts.csv has this patient), "no_recon" (the table is there but this patient isn't in it),
             or "no_contacts_csv" (there was no _cases/contacts.csv in the out folder, so nothing was looked up)
    """
    table = _load(str(Path(out_root) / "_cases" / "contacts.csv"))
    if table is None:
        return "no_contacts_csv"
    return "recon" if (table["study_id"].str.lower() == str(study_id).lower()).any() else "no_recon" #matched like contact_coords does


def _lead_positions(out_root, study_id, block, stage):
    #{lead number: "depthN" or "final"} for the leads that were in the brain during this recording
    m = re.match(r"increment(\d)_?depth(\d+)", block, re.IGNORECASE) #uh3002 writes increment1depth3, no underscore
    if m: #inching: this lead is at depth N, any earlier lead is already at its final position
        lead, depth = int(m.group(1)), m.group(2)
        return {**{k:"final" for k in range(1, lead)}, lead:f"depth{depth}"}
    if stage and stage != "stage1": #stage 2: everything implanted is at its final position
        return {1:"final", 2:"final"}
    leads = _leads_in_place(out_root, study_id, block, stage)
    return {k:"final" for k in range(1, leads + 1)} #0 = MER/baseline before any lead went in


def _leads_in_place(out_root, study_id, block, stage=""):
    #leads in place for this block from _cases/blocks.csv (OR sheet, cross-checked against impedance); 0 if unknown
    blocks = _load(str(Path(out_root) / "_cases" / "blocks.csv"))
    if blocks is None:
        return 0
    row = block_rows(blocks, study_id, block, stage)
    value = row["leads_in_place"].iloc[0] if len(row) == 1 else "" #none or ambiguous: unknown
    return int(float(value)) if value.strip() else 0


@lru_cache(maxsize=8)
def _load(path):
    #read a case table once per run (every recording of a patient uses the same tables)
    return pd.read_csv(path, dtype=str, keep_default_na=False) if Path(path).is_file() else None

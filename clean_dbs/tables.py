#one place for matching a recording to its row in _cases/blocks.csv, so every caller treats case and stage the same way


def block_rows(blocks, study_id, block, stage=""):
    """
    Input:
    blocks: the blocks.csv table (pandas, read as text)
    study_id, block: the recording (study ID matched without caring about case: Box has UG3003 where the tables have ug3003)
    stage: "stage1", "stage2", ... When the table has a stage column, only rows for that stage (or with the stage left blank) match,
        so a stage 1 and a stage 2 block with the same name (ug3007 erps01) don't mix. "" = don't filter by stage
    Returns: the matching rows (more than one = ambiguous; callers leave those out)
    """
    match = blocks[(blocks["study_id"].str.strip().str.lower() == str(study_id).strip().lower()) & (blocks["block"] == block)]
    if stage and "stage" in match:
        row_stage = match["stage"].str.replace(" ", "").str.lower()
        match = match[(row_stage == stage.replace(" ", "").lower()) | (row_stage == "")]
    return match

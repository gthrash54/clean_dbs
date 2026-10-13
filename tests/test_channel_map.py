from clean_dbs.channel_map import channel_table, leads_connected, leads_implanted


def test_lead_never_connected_is_labeled():
    #uh3016: GPi lead dropped in after bleeding, never plugged into the amplifier
    case = {"implant":"dual", "target1":"STN", "target2":"GPi", "notes":"GPi lead never connected (impedance shows nothing plugged in)"}
    assert leads_implanted(case) == (True, True)
    assert leads_connected(case) == (True, False)
    t = channel_table("uh3016/stage1 surgery/increment1_depth1", "ecos", "ecogdbs", 28, "STN", "GPi", "", leads_implanted(case), leads_connected(case))
    lead2 = [i for i, n in enumerate(t["name"]) if n.startswith("lead2")]
    assert {t["group"][i] for i in lead2} == {"lead2_not_connected"} and {t["target"][i] for i in lead2} == {""}
    assert t["target"][t["name"].index("lead1_c1")] == "STN"


def test_no_note_means_connected():
    assert leads_connected({"target1":"STN", "target2":"GPi", "notes":"GPi lead dropped in"}) == (True, True)
    assert leads_connected(None) == (True, True)

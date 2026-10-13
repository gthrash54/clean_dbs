import numpy as np
from clean_dbs.crosscheck import data_contacts, flat_contacts, parse_contacts, switched, _check_stim

#std (uV) of each lead's 8 contacts in real stim windows, copied from the cross-check cache (contact 1 first)
CASE = {"target1":"STN", "target2":"GPi"}


def evidence(*windows, channel=-1):
    #block_evidence-style dictionary: one stim window per argument, each {lead: 8 stds}
    return {"stim_read":True, "trains":[{"voice":"A", "time":float(i), "channel":channel, "std":{str(l):s for l, s in w.items()}}
                                        for i, w in enumerate(windows)]}


def test_lone_ring_8_found_next_to_an_unplugged_contact_on_the_other_lead():
    #uh3005 sedation_pre_pp_trial2, sheet "STN 234; GPi 8": lead 2 contact 8 reads like lead 1's unplugged 2-4
    w = {1:[18403, 5384, 5660, 5766, 19736, 17909, 18126, 11550], 2:[8938, 12642, 12760, 12656, 21160, 21396, 20967, 4992]}
    assert flat_contacts(w[2]) == [] #the per-lead rules alone miss it
    contacts, _ = data_contacts(evidence(w, w))
    assert contacts == {(1, 2), (1, 3), (1, 4), (2, 8)}


def test_lone_ring_1_found_when_it_drops_for_one_train():
    #uh3008 increment1_depth1_v2, sheet "1-2-3-4": contact 1 at 0.12 of its level in train 1, ~0.97 afterwards
    windows = [{1:[1043, 8802, 8784, 8755, 8748, 8755, 8762, 8755], 2:[1127, 1330, 1315, 1164, 1288, 1321, 1334, 1338]},
               {1:[8760, 1086, 9050, 8946, 9125, 9116, 9116, 9145], 2:[1162, 1395, 1369, 1214, 1344, 1390, 1387, 1392]},
               {1:[8746, 9048, 1150, 9480, 9122, 9135, 9115, 9138], 2:[1166, 1382, 1370, 1214, 1339, 1372, 1392, 1388]},
               {1:[8770, 9015, 9562, 1123, 9209, 9160, 9143, 9145], 2:[1168, 1382, 1368, 1214, 1337, 1371, 1385, 1387]}]
    contacts, order = data_contacts(evidence(*windows))
    assert contacts == {(1, 1), (1, 2), (1, 3), (1, 4)}
    assert order == [(1, 1), (1, 2), (1, 3), (1, 4)]


def test_low_plugged_rings_are_not_called_flat():
    #uh3016 increment1_depth1: plugged contact 8 far from the stim reads 0.15 of its level in every train
    far8 = [{1:[27955, 7128, 29915, 30070, 27998, 27666, 27741, 4069], 2:[7361, 7488, 7443, 7412, 7430, 7479, 7472, 7458]},
            {1:[27762, 29886, 6666, 29625, 25612, 27876, 26294, 4028], 2:[6797, 6898, 6884, 6842, 6849, 6892, 6901, 6885]},
            {1:[28052, 29844, 29443, 6751, 27680, 27310, 28088, 4055], 2:[6925, 7028, 7012, 6971, 6979, 7024, 7028, 7014]}]
    assert data_contacts(evidence(*far8))[0] == {(1, 2), (1, 3), (1, 4)}
    #uh3006 increment1_depth1_v2: contact 1 sits low in every train (it never drops for one train only)
    low1 = [{1:[5208, 4653, 37661, 46988, 14760, 14732, 15305, 17866], 2:[4658, 4738, 4731, 4676, 4727, 4741, 4744, 4742]},
            {1:[5011, 52455, 4708, 48495, 13808, 12878, 14272, 9367], 2:[4644, 4706, 4705, 4654, 4699, 4707, 4718, 4712]},
            {1:[6156, 71200, 69121, 6215, 21257, 17698, 18719, 16995], 2:[6167, 6224, 6221, 6164, 6213, 6224, 6231, 6228]}]
    assert data_contacts(evidence(*low1))[0] == {(1, 2), (1, 3), (1, 4)}


def test_mer_stim_through_lead_1_still_gives_2348():
    #uh3012 mer2_traj1: DBS1 234 into stim c1, contact 8 into stim c2 (the sheet's note)
    w = {1:[11191, 192, 82, 51, 14230, 11508, 12301, 383]}
    assert data_contacts(evidence(w, w))[0] == {(1, 2), (1, 3), (1, 4), (1, 8)}
    assert parse_contacts("DBS1 234 to stim c1; DBS1 8 to stim c2", "GPi", CASE)[0] == {(1, 2), (1, 3), (1, 4), (1, 8)}


def test_switch_experiment_flat_channels_are_not_evidence():
    #uh3002/uh3003 uh3_1dbs_directional: the stim log names stimulator channel 2, 3, 4, a switch routes the stim
    w = {1:[74227, 49778, 75043, 75696, 68795, 78386, 75786, 78851]}
    ev = evidence(w, channel=2)
    assert switched(ev) and data_contacts(ev) is None
    detail = []
    row = {"stim_contacts":"", "target":"GPi", "block":"increment1depth1"}
    assert _check_stim(row, ev, {"target1":"GPi"}, detail) == ("sheet_only", "sheet_only")
    assert "switch" in detail[0]
    assert not switched(evidence(w)) #noSwitch experiments log channel -1


def test_railed_or_misscaled_lead_is_skipped():
    #uh3002 dbs1 store: exact zeros next to 1e9 uV, nothing can be read from it
    w = {1:[1622331553, 0, 3025471628, 1363040679, 1389820218, 1930131274, 1904933370, 2133070139]}
    assert data_contacts(evidence(w)) is None
    w2 = {1:w[1], 2:[25954, 7163, 36601, 36694, 25240, 23981, 23530, 13713]} #the readable lead still counts
    assert data_contacts(evidence(w2))[0] == {(2, 2)}

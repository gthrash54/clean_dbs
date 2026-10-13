import contextlib
import csv
import io
import json
import re
import time
from pathlib import Path
import numpy as np
import pandas as pd
import tdt
from .block_info import describe_block
from .find_blocks import find_all_blocks
from .channel_map import connected_contacts, _pz5_stores
from .stim_events import stim_events
from .read_tdt import _events

#what the OR sheets (cases.csv / blocks.csv, typed from the handwritten sheets) are checked against, and how:
#  block exists: a tdt folder or .vhdr with that name under the data root, or on Box (uh3020/uh3024 aren't staged locally)
#  leads_in_place: impedance checks (PZ5 csvs / .vhdr) say which 8-channel lead groups were plugged in, counted up through the surgery
#  stim_contacts: in the noSwitch experiments the stimulated contact is unplugged from the PZ5, so its channel goes flat during stim
#    (not in the old uh3_1dbs_directional setup, where a switch routes the stim: those blocks are sheet_only, see switched)
#  target: the lead that was stimulated, lead 1 = cases.csv target1, lead 2 = target2
#every check gives match / conflict / sheet_only (nothing in the data to check it with) / no_data (recording not found)

#uh3 dual-lead ecos (28 ch): 1-12 ecog strips, 13-20 lead 1, 21-28 lead 2; mer blocks have 20 ch (no lead 2)
#within a lead: contact 1 = distal ring, 2-4 = first segmented level, 5-7 = second, 8 = proximal ring
#old uh3002/uh3003 setups record the lead on its own 8 ch store (dbs1, LFP1), ch k = contact k
LEAD_STORES = ("ecos", "dbs1", "LFP1")
#a lead contact is plugged in when its impedance is a real electrode value: lead contacts read 0.5-10 kOhm, while an empty input reads
#the 536.87 ceiling, 350-420 kOhm through an unconnected adapter (uh3012/13 baselines) or exactly 0 (uh3008 ecogimpedance)
LEAD_KOHM = (0.05, 100.0)
MIN_CONNECTED = 4 #a lead group counts as plugged in when at least 4 of its 8 contacts are
#an unplugged contact isn't truly flat during stim (it picks up some of the artifact), it's just much quieter than where it would be.
#tuned on uh3021, uh3005, uh3008 against sheets that are clear (increments "2-3-4", crosstalk "DBS1 234; DBS2 234", induction blocks):
#  one segment off: under 0.35 of the median of the other two on level 2-4 (flat 0.04-0.30 of them, plugged 0.8-1.2,
#    a slightly weak channel on a lead that isn't stimulating 0.42),
#    under 0.3 of them on level 5-7 (plugged level 5-7 segments next to a stimulated level read down to 0.33)
#  "234" all off: level 2 mean under 0.4 of level 3 mean (flat 0.03-0.30, plugged 0.9 and up; level 2 is never far from the stim otherwise)
#  any segment under 10% of the lead's loudest contact
#  rings 1 and 8 have no siblings and a plugged ring 8 can read under 10% of the max, so a ring only counts when a whole level is off
#    and the ring is no louder than 1.25x those segments (flat ring 0.9-1.2x, plugged ring 1.5x and up) or under 5% of the lead's max
#    (flat 0.025, plugged down to 0.099)
#  nothing is called flat on a lead whose loudest contact is under 2000 uV: that lead isn't carrying the stim (other lead, or no stim)
LEVEL_RATIO = {(2, 3, 4):0.35, (5, 6, 7):0.3}
GANG_RATIO, FLAT_RATIO, RING_RATIO, RING_MAX_RATIO = 0.4, 0.1, 1.25, 0.05
MIN_ARTIFACT_UV = 2000
#a ring unplugged on its own (contact 8 on the GPe in the uh3005 sedation runs, contact 1 in the uh3008 "1-2-3-4" run), see _lone_rings:
#  under 0.4 of its neighboring level (unplugged 0.12-0.25; a plugged ring far from the stim reads down to 0.15, so this alone isn't enough)
#  and within 0.8-1.25x of a contact that's off the amp on the other lead (unplugged 0.86-0.93x), or under 0.3x its own ratio in the
#  block's other windows (unplugged 0.12 vs 0.97). Checked on the 2026-10-12 cache (526 blocks): only those 4 blocks have a candidate at all
RING_ALONE_RATIO, TWIN_RANGE, RING_DROP = 0.4, (0.8, 1.25), 0.3
MAX_UV = 1e6 #lead std over 1 V = mis-scaled or railed store, see _readable
WINDOW_S = 0.5 #seconds read per stim window
MAX_WINDOWS = 6 #stim windows read per block, spread over its trains (each is a slow random read of the .tev on Box)
BOX_PAUSE_S = 0.3 #pause between reads on Box
_CASE_FOLDERS = {} #study ID -> its folder on the Box mount, filled as needed


def lead_groups(store, n_channels):
    """
    Input:
    store: lfp store name (ecos, dbs1, LFP1) or "dbs" for brainvision
    n_channels: channels in that store
    Returns: {lead number: list of channel numbers, contact k = k-th in the list}, empty if the layout isn't known
    """
    if store == "ecos" and n_channels in (20, 28):
        return {1:list(range(13, 21)), 2:list(range(21, 29))} if n_channels == 28 else {1:list(range(13, 21))}
    if store in ("dbs1", "LFP1") and n_channels == 8:
        return {1:list(range(1, 9))}
    if store == "dbs" and n_channels in (8, 9, 15, 16): #brainvision DBS1-8 = lead 1, DBS9-16 = lead 2 (passed the count of DBS channels)
        return {1:list(range(1, 9)), 2:list(range(9, 17))} if n_channels >= 15 else {1:list(range(1, 9))}
    return {}


def find_block(data_root, study_id, block, stage=""):
    """
    Input:
    data_root: folder with one folder per study ID (e.g. ~/DBS Data)
    study_id, block: as typed in blocks.csv (block name matched without caring about case)
    stage: the row's stage column ("stage1", "stage2"); "" = either
    Returns: the tdt block folder or .vhdr, stage1 first if both stages have one with that name and no stage is given; None if not there
    """
    found = [b for b in find_all_blocks(Path(data_root) / study_id, include_dnu=True) if describe_block(b)["block"].lower() == block.lower()]
    if stage:
        found = [b for b in found if describe_block(b)["stage"] == stage.replace(" ", "").lower()]
    found = sorted(found, key=lambda b: describe_block(b)["stage"]) #the sheets are stage 1 sheets, so stage1 wins a tie
    return found[0] if found else None


def read_crosswalk(path):
    """
    Input:
    path: crosswalk.csv (study_id, block, box_path), local only, the box paths start with a site number that never leaves this machine
    Returns: {(study_id, block lowercase): box_path}, duplicates dropped
    """
    if path is None or not Path(path).is_file():
        return {}
    with open(path, newline="") as f:
        return {(r["study_id"], r["block"].lower()):r["box_path"] for r in csv.DictReader(f)}


def find_box_block(box_root, crosswalk, study_id, block):
    """
    Input:
    box_root: the Box raw data folder as mounted (e.g. ~/Box/UAB DBS Data/UAB DBS Raw Data)
    crosswalk: from read_crosswalk
    study_id, block: as in blocks.csv
    Returns: path of that block on the mount, or None. Cases missing from the crosswalk are looked up by their study ID folder.
    """
    if box_root is None:
        return None
    box_root = Path(box_root)
    path = crosswalk.get((study_id, block.lower()))
    if path is not None and (box_root / path).exists():
        return box_root / path
    if study_id not in _CASE_FOLDERS: #site folders hold one folder per study ID; looked up once per case
        _CASE_FOLDERS[study_id] = next((c for top in sorted(p for p in box_root.iterdir() if p.is_dir())
                                        for c in (top / study_id, top / study_id.upper()) if c.is_dir()), None)
    case = _CASE_FOLDERS[study_id]
    if case is not None:
        for stage in sorted(case.iterdir()):
            if stage.is_dir() and re.search(r"stage\s*1", stage.name, re.I):
                hit = next((b for b in stage.iterdir() if b.name.lower() == block.lower() and any(b.glob("*.tsq"))), None)
                if hit is not None:
                    return hit
    return None


def block_evidence(study_id, block, data_root, box_root=None, crosswalk=None, cache=None, read_stim=True, stage=""):
    """
    Collects what the recording says about one blocks.csv row. Cached as json (one file per row), so reruns are fast.
    Input:
    study_id, block: the row
    data_root: local data (one folder per study ID)
    box_root: Box raw data mount, used for blocks that aren't local and for the stim windows (.tev isn't staged locally)
    crosswalk: from read_crosswalk
    cache: folder for the json files (None = no cache)
    read_stim: False skips the Box reads (leads only)
    Returns: dictionary with found (local, box or ""), format, stage, experiment, start_s (clock time of day in s, only for ordering),
      leads: {lead: contacts connected at the first impedance check} (missing = no check), leads_recorded: lead groups the streams can hold,
      trains: one per stim train (or window in it): voice, t0, t1, time (window middle), amplitude_ua, channel, std ({lead: std in uV of its 8 contacts})
    """
    #only the slow Box stim reads are cached; the local files are read again every time (cheap), so rule changes apply on rerun
    key = re.sub(r"[^A-Za-z0-9_.-]", "_", f"{study_id}__{block}" + (f"__{stage}" if stage not in ("", "stage1") else "")) + ".json" #stage1 keeps the old cache names
    cached = json.loads((Path(cache) / key).read_text()) if cache is not None and (Path(cache) / key).is_file() else {}
    path = find_block(data_root, study_id, block, stage)
    where = "local" if path is not None else ""
    box = find_box_block(box_root, crosswalk or {}, study_id, block) if box_root is not None else None
    if path is None and box is not None:
        path, where = box, "box"
    ev = {"study_id":study_id, "block":block, "found":where, "stim_read":False}
    if path is None:
        return _save(ev, cache, key)
    ev["format"] = "bv" if path.suffix.lower() == ".vhdr" else "tdt"
    ev["stage"] = "stage1" if where == "box" else describe_block(path)["stage"]
    ev.update(_bv_leads(path) if ev["format"] == "bv" else _tdt_leads(path))
    if cached.get("stim_read"):
        ev["trains"], ev["stim_read"] = cached["trains"], True
    elif ev["format"] == "tdt" and read_stim and box is not None:
        try:
            ev["trains"] = _stim_windows(box, ev.get("lead_store"), ev.get("lead_channels", 0))
            ev["stim_read"] = True
        except Exception as error: #a block tdt can't read (or a Box hiccup) is left unchecked, the next run tries again
            ev["stim_error"] = type(error).__name__
    return _save(ev, cache, key)


def _save(ev, cache, key):
    if cache is not None:
        Path(cache).mkdir(parents=True, exist_ok=True)
        (Path(cache) / key).write_text(json.dumps(ev))
    return ev


def _clock(text):
    #"9:52:29am" -> seconds since midnight (no date kept)
    m = re.search(r"(\d+):(\d\d):(\d\d)\s*([ap]m)", text, re.I)
    if m is None:
        return None
    h = int(m.group(1)) % 12 + (12 if m.group(4).lower() == "pm" else 0)
    return h * 3600 + int(m.group(2)) * 60 + int(m.group(3))


def _tdt_leads(block):
    out = {"experiment":"", "start_s":None, "leads":{}, "leads_recorded":0, "ecog_connected":None}
    listing = block / "StoresListing.txt"
    text = listing.read_text(errors="ignore") if listing.is_file() else ""
    m = re.search(r"Experiment:\s*(.*)", text)
    out["experiment"] = m.group(1).strip() if m else ""
    notes = block / "Notes.txt"
    if notes.is_file():
        m = re.search(r"Start:\s*(\S+)", notes.read_text(errors="ignore"))
        out["start_s"] = _clock(m.group(1)) if m else None
    stores = {name:n for name, n in _pz5_stores(block).values()}
    store = next((s for s in LEAD_STORES if f"Store ID : {s}" in text), None) #dbs1/LFP1 aren't whole sub-amps, so look in the listing
    n = stores.get(store, 8 if store in ("dbs1", "LFP1") else 0)
    out["lead_store"], out["lead_channels"] = store, n
    groups = lead_groups(store, n)
    out["leads_recorded"] = len(groups)
    checks = connected_contacts(block).get(store, None) if store else None
    if checks is not None and len(checks["channel"]):
        #first check that measured the lead channels (later passes often only measure the reference)
        times = [t for t in np.unique(checks["time"]) if np.isin(checks["channel"][checks["time"] == t], sum(groups.values(), [])).any()]
        if times:
            pick = checks["time"] == times[0]
            on = dict(zip(checks["channel"][pick].tolist(), _plugged(checks["impedance"][pick])))
            out["leads"] = {str(lead):[k + 1 for k, ch in enumerate(chans) if on.get(ch) == 1] for lead, chans in groups.items()}
            if store == "ecos":
                out["ecog_connected"] = sum(on.get(ch, 0) for ch in range(1, 13))
    return out


def _bv_leads(vhdr):
    out = {"experiment":"", "start_s":None, "leads":{}, "leads_recorded":0, "ecog_connected":None, "lead_store":"dbs", "lead_channels":0}
    names = [n.strip().upper() for n in re.findall(r"^Ch\d+=([^,]*),", vhdr.read_text(errors="ignore"), re.M)]
    dbs = {int(n[3:]):i + 1 for i, n in enumerate(names) if re.fullmatch(r"DBS\d+", n)} #DBS contact number -> channel number
    out["lead_channels"] = len(dbs)
    groups = lead_groups("dbs", len(dbs))
    out["leads_recorded"] = len(groups)
    vmrk = vhdr.with_suffix(".vmrk")
    if vmrk.is_file():
        m = re.search(r"New Segment,[^\n]*,(\d{8})(\d\d)(\d\d)(\d\d)", vmrk.read_text(errors="ignore"))
        out["start_s"] = int(m.group(2)) * 3600 + int(m.group(3)) * 60 + int(m.group(4)) if m else None
    checks = connected_contacts(vhdr).get("brainvision")
    if checks is not None and groups:
        pick = checks["time"] == checks["time"].max() #last check before the recording
        on = dict(zip(checks["channel"][pick].tolist(), _plugged(checks["impedance"][pick])))
        if any(dbs.get(c) in on for c in range(1, 17)):
            out["leads"] = {str(lead):[k + 1 for k, c in enumerate(contacts) if on.get(dbs.get(c)) == 1] for lead, contacts in groups.items()}
    return out


def _plugged(kohm):
    #1 where the impedance looks like an electrode (nan = "Out of Range!" etc. = 0)
    kohm = np.asarray(kohm, dtype=np.float64)
    return ((kohm > LEAD_KOHM[0]) & (kohm < LEAD_KOHM[1])).astype(int).tolist()


def _stim_windows(box, store, n_channels):
    #tsq-only read of the stim stores, then short windows of the lead store in the middle of each stim train
    with contextlib.redirect_stdout(io.StringIO()):
        headers = tdt.read_block(str(box), headers=1)
        raw = tdt.read_block(str(box), headers=headers, evtype=["epocs", "scalars"])
    time.sleep(BOX_PAUSE_S)
    #same events read_tdt gives, so the uh3001/uh3006 epocs with one onset too many are fixed the same way (read_tdt._match_epoc)
    result = _events(raw, headers)
    ev = stim_events(result, box)
    trains = []
    for voice in np.unique(ev["voice"]):
        pick = (ev["voice"] == voice) & ~(ev["amplitude_ua"] == 0) #amplitude 0 = logged but no current; nan (unknown) is kept
        t, amp, chn = ev["time"][pick], ev["amplitude_ua"][pick], ev["channel"][pick]
        if len(t) == 0:
            continue
        cut = np.flatnonzero(np.diff(t) > 1.0) + 1 #a pause over 1 s starts a new train
        for s, e in zip(np.r_[0, cut], np.r_[cut, len(t)]):
            trains.append({"voice":str(voice), "t0":float(t[s]), "t1":float(t[e - 1]), "amplitude_ua":float(np.nanmax(amp[s:e])) if np.isfinite(amp[s:e]).any() else None,
                           "channel":int(chn[s])})
    groups = lead_groups(store, n_channels)
    out = []
    if not trains or not groups:
        return [{k:v for k, v in tr.items()} for tr in trains]
    #spread the windows over the trains: every train gets one, long trains get more (contacts can change inside one)
    spots = []
    for i, tr in enumerate(trains):
        dur = tr["t1"] - tr["t0"]
        n = int(np.clip(dur // 30, 1, 2))
        spots += [(i, tr["t0"] + dur * (k + 1) / (n + 1)) for k in range(n)]
    if len(spots) > MAX_WINDOWS:
        spots = [spots[int(j)] for j in np.linspace(0, len(spots) - 1, MAX_WINDOWS)]
    for i, mid in spots:
        tr = trains[i]
        w = min(WINDOW_S, max(0.05, (tr["t1"] - tr["t0"]) * 0.8))
        with contextlib.redirect_stdout(io.StringIO()):
            raw = tdt.read_block(str(box), t1=max(0.0, mid - w / 2), t2=mid + w / 2, evtype=["streams"], store=store)
        time.sleep(BOX_PAUSE_S)
        if store not in raw.streams.keys(): #tdt StructType only answers "in" through keys()
            continue
        x = np.atleast_2d(raw.streams[store].data).astype(np.float64)
        sd = x.std(axis=1) * 1e6 #uV
        out.append({**tr, "time":round(float(mid), 2), "std":{str(l):sd[np.asarray(c) - 1].round(1).tolist() for l, c in groups.items()}})
    return out


def flat_contacts(std):
    """
    Input:
    std: std in uV of one lead's 8 contacts during a stim window (contact 1 first)
    Returns: list of contacts (1-8) that are flat, i.e. unplugged from the amp because they're on the stimulator (rules above)
    """
    s = np.asarray(std, dtype=np.float64)
    if len(s) != 8 or s.max() < MIN_ARTIFACT_UV:
        return []
    flat = {c for c in range(2, 8) if s[c - 1] < FLAT_RATIO * s.max()}
    for level, ratio in LEVEL_RATIO.items():
        for c in level:
            if s[c - 1] < ratio * np.median([s[o - 1] for o in level if o != c]):
                flat.add(c)
    if np.mean(s[1:4]) < GANG_RATIO * np.mean(s[4:7]):
        flat |= {2, 3, 4}
    for level in LEVEL_RATIO:
        if set(level) <= flat: #only next to a whole level that's off; one quiet segment in an increment sits close to a plugged ring 8
            flat |= {c for c in (1, 8) if s[c - 1] < RING_RATIO * max(s[o - 1] for o in level) or s[c - 1] < RING_MAX_RATIO * s.max()}
    return sorted(flat)


def parse_contacts(text, target, case):
    """
    Turns the sheet's stim_contacts into the lead contacts that should be unplugged during stim.
    Input:
    text: stim_contacts as typed ("2-3-4", "DBS1 234; DBS2 234", "Stim1 STN 234; Stim2 GPi 234; Stim3 GPi 8", "GPi 234; lead2 contact 8", ...)
    target: the row's target (picks the lead when the text doesn't say)
    case: the cases.csv row (target1 = lead 1, target2 = lead 2)
    Returns: (contacts, order): contacts = set of (lead, contact); order = list of single contacts when written as a sequence ("2-3-4"), else None
      None instead of the tuple if the notation can't be checked this way (bipolar "4+ 1-", patch, directional letters)
    "2-3-4" = one contact at a time in that order, "234" = 2, 3 and 4 together; GPe is lead 2's dorsal contact, so it counts as the GPi lead
    """
    text = text.strip()
    if not text or re.search(r"\+|patch|gnd|[a-c]-|abc", text, re.I):
        return None
    targets = [_norm_target(case.get("target1", "")), _norm_target(case.get("target2", ""))]
    default = _lead_of(_norm_target(target), targets)
    contacts, order = set(), []
    for part in text.split(";"):
        part = re.sub(r"to stim c\d", "", part, flags=re.I) #"DBS1 234 to stim c1" = which stimulator output, not a contact
        lead = default
        m = re.search(r"DBS\s*([12])|lead\s*([12])", part, re.I)
        if m:
            lead = int(m.group(1) or m.group(2))
        else:
            for word in re.findall(r"STN|GPi|GPe|GP\b", part, re.I):
                lead = _lead_of(_norm_target(word), targets) or lead
        part = re.sub(r"Stim\d|DBS\s*\d|lead\s*\d|contact|STN|GPi|GPe|GP\b|crash|/", " ", part, flags=re.I)
        tokens = re.findall(r"c?(\d+)", part)
        if lead is None or not tokens:
            return None
        sequence = "-" in part and all(len(t) == 1 for t in tokens)
        for t in tokens:
            contacts.update((lead, int(c)) for c in t)
            if sequence:
                order.append((lead, int(t)))
    if any(c < 1 or c > 8 for _, c in contacts):
        return None
    return contacts, (order or None)


def _norm_target(t):
    t = t.strip().upper()
    return "GPI" if t in ("GPE", "GP", "GPI") else t


def _lead_of(target, targets):
    #lead 1 = target1, lead 2 = target2; a target that names both (STN+GPi) or neither gives None
    return targets.index(target) + 1 if target in targets and target else None


def data_contacts(ev):
    """
    Input:
    ev: from block_evidence
    Returns: (contacts, order): set of (lead, contact) flat during stim, and the order in which single flat contacts showed up;
      None if no stim windows could be read
    """
    windows = [w for w in ev.get("trains", []) if "std" in w]
    if not ev.get("stim_read") or not windows or switched(ev):
        return None
    #per window: {lead: (std, flat contacts)}, starting from the rules in flat_contacts; railed or mis-scaled leads are left out
    per = []
    for w in sorted(windows, key=lambda w: w["time"]):
        leads = {int(l):np.asarray(std, dtype=np.float64) for l, std in w["std"].items()}
        per.append({l:(s, set(flat_contacts(s))) for l, s in leads.items() if len(s) == 8 and _readable(s)})
    if not any(per):
        return None
    _lone_rings(per)
    contacts, order = set(), []
    for leads in per:
        now = {(l, c) for l, (_, flat) in leads.items() for c in flat}
        contacts |= now
        if len(now) == 1 and (not order or order[-1] != next(iter(now))):
            order.append(next(iter(now)))
    return contacts, order


def switched(ev):
    """
    True when the stim log names a stimulator channel above 1 (the old uh3_1dbs_directional setup, uh3002/uh3003: channel 2, 3, 4 in turn).
    There a switch routes the stim and the stimulated contact stays on the amp, so a flat channel isn't evidence of anything.
    noSwitch experiments log channel -1 (unknown) or 1.
    """
    return any(w.get("channel", -1) >= 2 for w in ev.get("trains", []))


def _readable(s):
    #a lead whose std reads exactly 0 or over 1 V on any contact is railed or mis-scaled (uh3002 dbs1 reads 1e9 uV), it can't say what's flat
    return bool(np.all(np.isfinite(s)) and s.min() > 0 and s.max() < MAX_UV)


def _ring_ratio(s, ring, level):
    return s[ring - 1] / max(np.median([s[c - 1] for c in level]), 1e-9)


def _lone_rings(per):
    """
    Adds a ring (contact 1 or 8) that was unplugged on its own, which flat_contacts can't see (a plugged ring 8 far from the stim reads
    down to 0.15 of its level, the same as an unplugged one). Only on a lead carrying the stim with nothing else flat, when the ring is under
    RING_ALONE_RATIO of its neighboring level and either
      reads like a contact that is off the amp on the other lead in the same window (unplugged inputs share one floor: uh3005 sedation
        pre_pp_trial2 lead2 c8 4992 uV next to lead1 c2-4 5384-5766), or
      drops under RING_DROP of its own ratio in the block's other windows (unplugged for one train of a sequence: uh3008 increment1_depth1_v2
        c1 0.12 of its level, 0.97 in the next three trains)
    Input: per = list of windows in time order, each {lead: (std, set of flat contacts)}; the sets are updated in place
    """
    flat_before = [{l:set(f) for l, (_, f) in leads.items()} for leads in per]
    for i, leads in enumerate(per):
        for l, (s, flat) in leads.items():
            if flat_before[i][l] or s.max() < MIN_ARTIFACT_UV:
                continue
            floor = [per[i][o][0][c - 1] for o in leads if o != l for c in flat_before[i][o]]
            for ring, level in ((1, (2, 3, 4)), (8, (5, 6, 7))):
                r = _ring_ratio(s, ring, level)
                if r >= RING_ALONE_RATIO:
                    continue
                twin = bool(floor) and TWIN_RANGE[0] <= s[ring - 1] / np.median(floor) <= TWIN_RANGE[1]
                own = [_ring_ratio(p[l][0], ring, level) for j, p in enumerate(per) if j != i and l in p and p[l][0].max() >= MIN_ARTIFACT_UV]
                drop = bool(own) and r < RING_DROP * np.median(own)
                if twin or drop:
                    flat.add(ring)


def _fmt(contacts):
    by = {}
    for lead, c in sorted(contacts):
        by.setdefault(lead, []).append(str(c))
    return "; ".join(f"lead{l} c{''.join(cs)}" for l, cs in by.items()) or "none"


def crosscheck(cases, blocks, data_root, box_root=None, crosswalk=None, cache=None, read_stim=True, progress=False):
    """
    Checks every row of the case tables against the recordings.
    Input:
    cases, blocks: the tables as read with pd.read_csv(..., dtype=str, keep_default_na=False)
    data_root, box_root, crosswalk, cache, read_stim: see block_evidence
    progress: print one line per block (study ID and block only)
    Returns: (cases, blocks) copies with new columns
      blocks: check_block, check_leads_in_place, check_stim_contacts, check_target, status, check_detail
      cases: check_implant, status, check_detail
      each check is match, conflict, sheet_only (the data can't say) or no_data (recording not found);
      status is no_data if the block isn't found, conflict if any check conflicts, match if at least one matches, else sheet_only
    """
    cases, blocks = cases.copy(), blocks.copy()
    evidence = []
    for i, r in blocks.iterrows():
        if progress:
            print(f"{i + 1}/{len(blocks)} {r['study_id']} {r['block']}", flush=True)
        evidence.append(block_evidence(r["study_id"], r["block"], data_root, box_root, crosswalk, cache, read_stim, r.get("stage", "")))
    case_rows = {r["study_id"]:r for _, r in cases.iterrows()}
    #leads counted up through each stage: a lead seen plugged in stays in place for every later block of that stage
    seen = [_leads_seen(ev) for ev in evidence]
    upto, most = [None] * len(evidence), [None] * len(evidence)
    for (sid, stage), idx in pd.Series(range(len(evidence))).groupby([[ev.get("study_id") for ev in evidence], [ev.get("stage", "") for ev in evidence]]).groups.items():
        top = max((seen[i] for i in idx if seen[i] is not None), default=None) #None = no block of this stage has a check
        for i in idx:
            most[i] = top #most leads plugged in at any point this stage
        idx = sorted((i for i in idx if evidence[i].get("start_s") is not None), key=lambda i: evidence[i]["start_s"])
        best = 0
        for i in idx:
            best = max(best, seen[i] or 0)
            upto[i] = best
    out = {c:[] for c in ("check_block", "check_leads_in_place", "check_stim_contacts", "check_target", "status", "check_detail")}
    for i, (_, r) in enumerate(blocks.iterrows()):
        ev, case = evidence[i], case_rows.get(r["study_id"], {})
        detail = []
        if not ev["found"]:
            checks = ["no_data"] * 4
            detail.append("recording not found locally or on Box")
        else:
            checks = ["match", _check_leads(r, ev, seen[i], upto[i], most[i], detail), *_check_stim(r, ev, case, detail)]
        status = "no_data" if checks[0] == "no_data" else "conflict" if "conflict" in checks else "match" if "match" in checks[1:] else "sheet_only"
        for c, v in zip(list(out)[:4], checks):
            out[c].append(v)
        out["status"].append(status)
        out["check_detail"].append("; ".join(detail))
    for c, v in out.items():
        blocks[c] = v
    #implant: dual if any block of the case had both lead groups plugged in or stimulated, single if leads were seen but never two
    implant, status, detail = [], [], []
    for _, r in cases.iterrows():
        mine = [i for i, ev in enumerate(evidence) if ev.get("study_id") == r["study_id"] and ev["found"]]
        both = [i for i in mine if (seen[i] or 0) >= 2 or {l for l, _ in (data_contacts(evidence[i]) or (set(), []))[0]} == {1, 2}]
        any_lead = [i for i in mine if seen[i] is not None or data_contacts(evidence[i]) is not None]
        if not mine:
            implant.append("no_data"); status.append("no_data"); detail.append("no recordings found")
            continue
        if not any_lead:
            implant.append("sheet_only"); status.append("sheet_only"); detail.append("no impedance checks or stim windows to count leads")
            continue
        data = "dual" if both else "single"
        sheet = r.get("implant", "").strip().lower()
        ok = sheet == data
        implant.append("match" if ok else "conflict")
        status.append("match" if ok else "conflict")
        why = f"two leads plugged in or stimulated in {len(both)} block(s), e.g. {evidence[both[0]]['block']}" if both else f"never more than one lead group connected or stimulated ({len(any_lead)} blocks with lead evidence)"
        detail.append(("" if ok else f"sheet {sheet or 'blank'} vs data {data}: ") + why)
    cases["check_implant"], cases["status"], cases["check_detail"] = implant, status, detail
    return cases, blocks


def _leads_seen(ev):
    """
    Leads the block itself shows: lead groups plugged in at its impedance check (tdt PZ5 csv or .vhdr), and any lead with a contact
    that went flat during stim (it has to be in the brain and on the stimulator). None = this block can't say.
    """
    counts = []
    if ev.get("leads"):
        counts.append(sum(len(c) >= MIN_CONNECTED for c in ev["leads"].values()))
    stim = data_contacts(ev)
    if stim is not None and stim[0]:
        counts.append(max(l for l, _ in stim[0])) #lead 2 stimulated means two leads are in
    return max(counts) if counts else None


def _check_leads(r, ev, seen, upto, most, detail):
    """
    leads_in_place against the impedance checks. The data gives a range: at least the leads plugged in so far this stage
    (a lead doesn't come out mid-surgery), at most the leads ever plugged in this stage.
    sheet below the low end or above the high end = conflict; inside it = match when this block has its own check, else sheet_only
    (a lead can sit in place unplugged, e.g. on the rep's clinical device, so a block's own check alone can't lower the count)
    """
    sheet = r["leads_in_place"].strip()
    if not sheet.isdigit():
        return "sheet_only"
    sheet = int(sheet)
    low = max(seen or 0, upto or 0)
    if seen is None and upto is None and not most:
        return "sheet_only"
    if sheet < low:
        detail.append(f"leads_in_place: sheet {sheet}, but {low} lead(s) plugged in by this block" + (f" ({_plugged_text(ev)})" if seen is not None else ""))
        return "conflict"
    if most is not None and sheet > most:
        detail.append(f"leads_in_place: sheet {sheet}, but never more than {most} lead(s) plugged in this stage")
        return "conflict"
    return "match" if seen is not None or upto else "sheet_only"


def _plugged_text(ev):
    return "impedance: " + ", ".join(f"lead{l} {len(c)}/8" for l, c in ev["leads"].items())


def _check_stim(r, ev, case, detail):
    #returns (check_stim_contacts, check_target)
    if ev.get("format") == "bv":
        return "sheet_only", "sheet_only" #brainvision rigs: contacts are switched outside the recording, nothing flat to see
    sheet_text = r["stim_contacts"].strip()
    if switched(ev):
        detail.append("stim_contacts: stim routed by a switch (stim log names stimulator channel 2 or up), flat channels aren't evidence")
        return "sheet_only", "sheet_only"
    data = data_contacts(ev)
    trains = ev.get("trains", [])
    if data is None:
        if not ev.get("stim_read"):
            return "sheet_only", "sheet_only"
        if sheet_text and not trains and parse_contacts(sheet_text, r["target"], case) is None:
            return "sheet_only", "sheet_only" #bipolar / patch notation = an external stimulator, which tdt doesn't log
        if sheet_text and not trains:
            detail.append(f"stim_contacts: sheet {sheet_text}, but the recording has no stim pulses")
            return "conflict", "sheet_only"
        return ("match" if not sheet_text and not trains else "sheet_only"), "sheet_only"
    contacts, order = data
    targets = [case.get("target1", ""), case.get("target2", "")]
    #target: which lead(s) carried the stim
    target_check = "sheet_only"
    if contacts and "mer" not in r["block"].lower(): #in mer blocks the target is where the microelectrode is, stim can come from the other lead
        leads = {l for l, _ in contacts}
        data_targets = {_norm_target(targets[l - 1]) for l in leads if targets[l - 1]}
        sheet_targets = {_norm_target(t) for t in r["target"].split("+") if t.strip()}
        if data_targets and sheet_targets:
            target_check = "match" if data_targets == sheet_targets or (len(leads) == 1 and data_targets <= sheet_targets and len(sheet_targets) == 1) else "conflict"
            if target_check == "conflict":
                detail.append(f"target: sheet {r['target']}, stim was on {' and '.join(f'lead{l} ({targets[l - 1]})' for l in sorted(leads))}")
    if not sheet_text:
        if contacts:
            detail.append(f"stim_contacts: sheet blank, data flat during stim: {_fmt(contacts)}")
            return "conflict", target_check
        return "match", target_check
    parsed = parse_contacts(sheet_text, r["target"], case)
    if parsed is None:
        return "sheet_only", target_check
    if not contacts:
        loud = max((max(max(v) for v in w["std"].values()) for w in trains if "std" in w), default=0)
        detail.append(f"stim_contacts: sheet {sheet_text}, no lead contact went flat during stim (loudest lead channel {loud:.0f} uV)")
        return "conflict", target_check
    want, want_order = parsed
    if want == contacts:
        if want_order and len(want_order) == len(set(want_order)) and order and [c for c in order if c in want_order] != want_order:
            detail.append(f"stim_contacts: same contacts, order in data {'-'.join(str(c) for _, c in order)}")
        return "match", target_check
    detail.append(f"stim_contacts: sheet {sheet_text} ({_fmt(want)}), data {_fmt(contacts)}")
    return "conflict", target_check

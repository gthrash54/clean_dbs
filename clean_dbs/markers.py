import re
from pathlib import Path
import numpy as np
from clean_dbs.bv_channels import cohort_of
from clean_dbs.flags import new_segment_gaps
from clean_dbs import sidecars as sc

#Brainvision S/R markers are the state of the trigger port (S = Stimulus input, R = Response input), written every time it changes.
#Each port line is a TTL: the lab's STG stimulator sync output, the Buzzy vibration box, E-Prime, ... So "S15 then S13 60 us later" is
#one line (value 2) going high then low. Markers are read as edges: a code with bits the previous one didn't have = rising edge,
#a code that only drops bits = falling edge. The time a line stays high (width) tells what it is (survey of 1,940 .vmrk on Box, 2026-10-12):
#  width 0.03-0.5 ms                       = "pulse": the STG sync follows each stim pulse (width = pulse width, 60-240 us on sdcurves)
#  width 0.1-40 ms, 3+ different widths    = "ttl_code": one TTL per stim trial, its width codes the condition (TTL.txt, stim_struct.mat)
#  width 0.5-40 ms, one or two widths      = "fixed_ttl": could be stim (DLPFC itbs: 2 ms every 20 ms) or a task trigger (E-Prime), so not stim
#                                            unless a rule or a stimulator file says so
#  width >= 40 ms                          = "long": not a stim pulse (100 ms = the Buzzy vibration in UG/BB/JOY, task cues elsewhere)
#  no falling edge                         = "single"
MAX_SYNC_MS = 40.0 #longest stim TTL code seen (pp stim_TTL.txt goes to 26 ms)
PULSE_MS = 0.5 #longest pulse-width sync seen (sdcurves go to 0.24 ms)
MIN_GAP_MS = 0.05 #two stim rises closer than this can't be pulses (paired pulses go down to 0.1 ms apart)

#what a code means in a cohort, with the evidence. cohort, code, kind (from the widths) -> type, meaning, evidence, confirmed
#confirmed = 1 when checked against a stimulator file or the OR notes, 0 when it only follows the generic width rule or a pattern;
#an optional last entry is the fewest rising edges the rule needs
RULES = [
    ("UH", "R15", "pulse", "stim", "stim_pulse_start", "uh3021 STN_dbs1: R15 = stimulation protocol.dat pulse onsets within 1 sample, R11 60 us later (2026-10-11); uh3010 contact2: R15 = all 1921 STG sync outputs of the BRAINaim2 protocol .dat, max error 0.02 ms (2026-10-12)", 1),
    ("UH", "S15", "pulse", "stim", "stim_pulse_start", "uh3003/3004/3013/3017: S15 -> S13 after 0.06 ms in 32,432 of 32,898 pairs = the 60 us pulse; same STG protocol as uh3010/uh3021 (2026-10-12)", 1),
    ("UG", "S15", "pulse", "stim", "stim_pulse_start", "UG3010 stage2 erp1: S15 = all 13,146 sync outputs of erp1_events.dat, max error 0.023 ms (2026-10-12)", 1),
    ("UG", "S7", "pulse", "stim", "stim_pulse_start", "ug3001 erpstim1-9 / sdcurve1: S7 -> S5 after the pulse width (0.06 ms; 0.06-0.24 ms on the sdcurve), same STG setup as UG3010 (2026-10-12)", 1),
    ("UG", "S13", "single", "stim", "stim_pulse_start", "ug3003 stage2 erp001: S13 only ever rises (9,844 times) and its spacings are exactly the BRAINaim1 ERP program's (0.3/0.8/2/30 ms pairs, 7.69 ms bursts, 50 and 200 ms gaps); the program file next to it is from another run (dated about two weeks earlier, same design in a different order, best alignment 1,811 of 9,844 at the 2 ms limit), so no per-pulse settings (2026-10-12). Needs 100+ rises: one stray S13 (ug3007 baseline_lfp01) isn't stim", 0, 100),
    ("UG", "S7", "long", "task", "buzz_start", "ug3001 stage1 mer and LFP_PRE: every S7 = a 100 ms sync of 'stim files/buzz_stim.dat' (156/156 and 139/139 matched, max error 0.09 ms): the Buzzy vibration, not DBS (2026-10-12)", 1),
    ("UG", "S15", "long", "task", "buzz_start", "S15 held exactly 100 ms at a jittered ~2 Hz, same as the buzz_stim.dat program next to it (UG3003 LFP_baseline, UG3007 baseline_buzz); the Buzzy vibration, not DBS", 0),
    ("BB", "S15", "long", "task", "buzz_start", "BuzzyBee cohort: buzzy_bee_stim_TTL.txt says pulse_width_ms=100, buzz_stim.dat plays 100 ms buzzes; S15 is held 100 ms", 0),
    ("BB", "S7", "long", "task", "buzz_start", "BuzzyBee cohort: buzzy_bee_stim_TTL.txt says pulse_width_ms=100, buzz_stim.dat plays 100 ms buzzes; S7 is held 100 ms", 0),
    ("MC", "*", "ttl_code", "stim", "stim_trial_start", "MC011 0-1, 1-2, 2-3: S15 = all 1700 sync outputs of results.txt (max error 0.02 ms); the stim pulse starts on the same row, and the S15 -> S7 width = the TTL that TTL.txt maps to the condition (2026-10-12)", 1),
    ("PP", "S13", "pulse", "stim", "stim_pulse_start", "PP017 and PP023: S13 = all 7200 sync outputs of results.txt (max error 0.01 ms), S12 0.06 ms later (2026-10-12)", 1),
]
#generic meaning by width class when no rule matches
GENERIC = {"pulse":("stim", "stim_pulse_start", "generic: line high for a pulse width (0.03-0.5 ms) like the STG sync on UH/UG/PP; not checked for this cohort", 0),
           "ttl_code":("stim", "stim_trial_start", "generic: line high 0.1-40 ms with many widths, the lab's TTL-coded stim trials (checked on MC011); not checked for this cohort", 0),
           "fixed_ttl":("marker", "unknown", "line held 0.5-40 ms with one fixed width: stim or a task trigger, can't tell without a rule or a stimulator file", 0),
           "long":("marker", "unknown", "line held 40 ms or more: not a stim pulse; meaning not documented", 0),
           "single":("marker", "unknown", "no falling edge on this line; meaning not documented", 0)}
#non-trigger markers brainvision writes itself
FIXED = {"NewSegment":("other", "new_segment", "brainvision recording (re)start, see flags/new_segments"),
         "BtnPressed":("marker", "recorder_button_press", "BrainVision 'Comment' marker written when the recorder's push button is pressed (paired with BtnReleased 0.1-0.3 s later)"),
         "BtnReleased":("marker", "recorder_button_release", "BrainVision 'Comment' marker written when the recorder's push button is released")}
#leads for the lab, never used as meanings (unknown stays unknown until someone confirms)
HINTS = {("DLPFC", "R15"):"DLPFC047/091: R15 (or R14) count = 141 = Go_NoGo log frames (138 trials + 3 blocks); candidate E-Prime trial onset, needs the .es3 script",
         ("DLPFC", "R14"):"DLPFC047/091: R15 (or R14) count = 141 = Go_NoGo log frames (138 trials + 3 blocks); candidate E-Prime trial onset, needs the .es3 script",
         ("JOY", "R15"):"JOY: R15 -> R13 held ~100 ms, ~4 s apart; candidate task cue from the reach task, needs the task script"}


def read_vmrk(vhdr):
    """
    Reads the markers of a recording without the signal (for checking codes on headers only).
    Input: a brainvision .vhdr (its .vmrk next to it)
    Returns: epocs like read_bv: {key: {"onset", "offset", "marker"}}, times in s from the first sample
    """
    vhdr = Path(vhdr)
    header = vhdr.read_text(errors="ignore")
    vmrk = vhdr.parent / next((l.split("=", 1)[1].strip() for l in header.splitlines() if l.startswith("MarkerFile=")), vhdr.stem + ".vmrk")
    dt = float(re.search(r"^SamplingInterval=(\S+)", header, re.M).group(1)) / 1e6
    out = {}
    for m in re.finditer(r"^Mk\d+=([^,]*),([^,]*),(\d+),(\d+)", vmrk.read_text(errors="ignore") if vmrk.is_file() else "", re.M):
        name = f"{m.group(1)}/{m.group(2)}"
        key = m.group(2).replace(" ", "") or m.group(1).replace(" ", "")
        e = out.setdefault(key, {"onset":[], "offset":[], "marker":name})
        e["onset"].append((int(m.group(3)) - 1) * dt)
        e["offset"].append((int(m.group(3)) - 1 + int(m.group(4))) * dt)
    return {k:{"onset":np.array(v["onset"]), "offset":np.array(v["offset"]), "marker":v["marker"]} for k, v in out.items()}


def edges(epocs):
    """
    Input: epocs from read_bv / read_vmrk
    Returns: dictionary of arrays, one entry per S/R marker in time order: time, key (epoc name), port (S or R), code, rise (1 = gained a bit),
             width_ms (how long the bit it raised stayed high, until the first of its raised bits fell; nan = never fell / not a rise)
    A code equal to the one before is a rise too: the port went back to 0 in between, which brainvision doesn't mark.
    """
    rows = []
    for key, e in epocs.items():
        m = re.fullmatch(r"([SR])(\d+)", key)
        if m:
            rows += [(t, key, m.group(1), int(m.group(2))) for t in np.asarray(e["onset"], dtype=np.float64)]
    rows.sort()
    n = len(rows)
    out = {"time":np.array([r[0] for r in rows]), "key":np.array([r[1] for r in rows], dtype=object), "port":np.array([r[2] for r in rows], dtype=object),
           "code":np.array([r[3] for r in rows], dtype=np.int64), "rise":np.zeros(n, dtype=np.int8), "width_ms":np.full(n, np.nan)}
    for port in ("S", "R"):
        state, open_bits = 0, {}
        for i in np.flatnonzero(out["port"] == port):
            code = out["code"][i]
            gained = code & ~state if code != state else code
            lost = state & ~code if code != state else 0
            for b in _bits(lost):
                j = open_bits.pop(b, None)
                if j is not None and np.isnan(out["width_ms"][j]):
                    out["width_ms"][j] = (out["time"][i] - out["time"][j]) * 1e3
            if gained:
                out["rise"][i] = 1
                for b in _bits(gained):
                    open_bits[b] = i
            state = code
    return out


def _bits(x):
    return [1 << k for k in range(16) if x & (1 << k)]


def classify_code(widths, times):
    """
    Input: widths (ms) and times (s) of one code's rising edges
    Returns: kind: pulse, ttl_code, long or single (see the top of this file), plus "fast" if rises come closer than MIN_GAP_MS (then never stim)
    """
    w = np.asarray(widths, dtype=np.float64)
    fin = w[np.isfinite(w)]
    if len(w) == 0 or len(fin) < 0.8 * len(w):
        return "single"
    if len(times) > 1 and np.min(np.diff(np.sort(times))) * 1e3 < MIN_GAP_MS:
        return "fast"
    if np.percentile(fin, 95) <= PULSE_MS:
        return "pulse"
    if np.percentile(fin, 95) <= MAX_SYNC_MS:
        values, counts = np.unique(np.round(fin, 1), return_counts=True) #ttl codes step 0.1 ms or more
        return "fixed_ttl" if counts.max() >= 0.9 * len(fin) or len(values) < 3 else "ttl_code"
    return "long"


def label_bv(result, vhdr, cohort=None):
    """
    Works out what every brainvision marker code means for this recording, and builds the stim table from the codes that are stim.
    Input:
    result: the dictionary from read_bv (its epocs get type, meaning, evidence, kind and hint added in place)
    vhdr: the .vhdr (for the cohort, the .vmrk clock and sidecars next to it)
    cohort: default from the path (bv_channels.cohort_of)
    Returns: {"stim": stim table (see stim_events.py, plus marker_lag_ms = how long the line stayed high, and sidecar columns when one matched),
              "sidecars": {name: entry from sidecars.attach_sidecars}}
    Steps: 1) each code's rising edges get a width class and a rule (RULES, else GENERIC); 2) codes that aren't stim yet are lined up with
    the stimulator files next to the recording: a match with a stim program makes them stim (confirmed), a buzz program makes them buzz_start;
    3) the rising edges of stim codes become the stim table, and the sidecars that line up add their per-pulse settings.
    Only rising edges become stim rows; the falling partner (S13 after S15 ...) is dropped, never a second pulse.
    """
    from clean_dbs.stim_events import _table, _stack, bv_clock_time
    vhdr = Path(vhdr)
    cohort = cohort_of(vhdr) if cohort is None else cohort
    epocs = result["epocs"]
    ed = edges(epocs)
    gaps = new_segment_gaps(vhdr) if vhdr.is_file() else None
    def clock(t): #brainvision drops time at real gaps, sidecars don't
        return bv_clock_time(t, gaps) if gaps is not None and len(gaps["time"]) else t
    confirmed = {}
    for key, e in epocs.items():
        fixed = FIXED.get(key, FIXED.get(key.split("_")[-1]))
        if fixed:
            e["type"], e["meaning"], e["evidence"] = fixed
            continue
        if not re.fullmatch(r"[SR]\d+", key):
            e["type"], e["meaning"], e["evidence"] = "marker", "unknown", "not a trigger code"
            continue
        mine = np.flatnonzero(ed["key"] == key)
        rises = mine[ed["rise"][mine] == 1]
        if len(rises) < 0.5 * len(mine): #mostly falling edges: the end of another code's TTL
            same_port = np.flatnonzero(ed["port"][:mine[0]] == ed["port"][mine[0]]) #the S and R ports are separate lines
            before = ed["key"][same_port[-1]] if len(same_port) else ""
            e.update(type="marker", meaning="ttl_end", kind="end", evidence=f"falling edge: the TTL raised by {before} goes low here (kept as a marker, never a pulse)")
            continue
        kind = classify_code(ed["width_ms"][rises], ed["time"][rises])
        #an optional 8th entry is the fewest rises the rule needs
        rule = next((r for r in RULES if r[0] == cohort and r[1] in (key, "*") and r[2] == kind and len(rises) >= (r[7] if len(r) > 7 else 0)), None)
        etype, meaning, evidence, ok = rule[3:7] if rule else GENERIC.get(kind, ("marker", "unknown", f"rises closer than {MIN_GAP_MS} ms apart, can't be pulses", 0))
        e.update(type=etype, meaning=meaning, evidence=evidence, kind=kind)
        confirmed[key] = ok
        if (cohort, key) in HINTS:
            e["hint"] = HINTS[(cohort, key)]
    #2) codes that aren't stim: do their rising edges line up with a stimulator program next to the recording?
    files = {}
    for f in sc.find_sidecars(vhdr) if vhdr.is_file() else []:
        try:
            kind, table = sc.read_sidecar(f)
        except Exception: #a broken sidecar is reported by attach_sidecars below
            continue
        if kind in ("mcs", "protocol_csv"):
            files[f.name] = (kind, table)
    for key, e in epocs.items():
        if e.get("type") not in ("marker", "task") or e.get("kind") not in ("single", "long", "ttl_code", "fixed_ttl"):
            continue
        rises = (ed["key"] == key) & (ed["rise"] == 1)
        t = clock(ed["time"][rises])
        for name, (kind, table) in files.items():
            m = sc.match_times(t, sc.event_times(kind, table))
            if not m["matched"]:
                continue
            e["evidence"] = f"{e['evidence']}; its times match {name} ({m['n_pairs']} of {len(t)}, max error {m['max_err_ms']:.3f} ms)"
            if "buzz" in name.lower():
                e["type"], e["meaning"] = "task", "buzz_start"
            elif kind == "protocol_csv" or np.nanmax(np.abs(table["amplitude"]), initial=0) > 0:
                e["type"], e["meaning"] = "stim", "stim_pulse_start" if e["kind"] == "single" else "stim_trial_start"
                confirmed[key] = 1
            break
    #2b) a code can carry only some of a program's sync outputs (DLPFC047 itbs: 40 S15 among 1,103 S7), too few to line up on their own;
    #it's stim when the program lines up with it and the stim codes together, and nearly all of its own rises pair up
    stim_keys = [k for k, e in epocs.items() if e.get("type") == "stim"]
    for key, e in epocs.items():
        if not stim_keys or e.get("type") not in ("marker", "task") or e.get("kind") not in ("single", "long", "ttl_code", "fixed_ttl"):
            continue
        rises = ed["rise"] == 1
        union = rises & np.isin(ed["key"], stim_keys + [key])
        mine = (ed["key"][union] == key)
        for name, (kind, table) in files.items():
            if "buzz" in name.lower() or not (kind == "protocol_csv" or np.nanmax(np.abs(table["amplitude"]), initial=0) > 0):
                continue
            m = sc.match_times(clock(ed["time"][union]), sc.event_times(kind, table))
            paired = int(mine[m["pairs"][:, 0]].sum()) if m["matched"] else 0
            if paired >= 3 and paired >= 0.9 * mine.sum():
                e["evidence"] = (f"{e['evidence']}; together with {', '.join(stim_keys)} its times match {name} ({paired} of its {int(mine.sum())} rises, "
                                 f"max error {m['max_err_ms']:.3f} ms)")
                e["type"], e["meaning"] = "stim", "stim_pulse_start" if e["kind"] == "single" else "stim_trial_start"
                confirmed[key] = 1
                break
    #3) stim table from the rising edges of stim codes
    rows = []
    for key, e in epocs.items():
        if e.get("type") == "stim":
            rises = (ed["key"] == key) & (ed["rise"] == 1)
            table = _table(ed["time"][rises], key, "", {}, confirmed=confirmed.get(key, 0))
            table["marker_lag_ms"] = ed["width_ms"][rises]
            rows.append(table)
    stim = _stack(rows) if rows else _stack([])
    if "marker_lag_ms" not in stim:
        stim["marker_lag_ms"] = np.full(len(stim["time"]), np.nan)
    found = {}
    if vhdr.is_file():
        work = dict(stim, time=clock(stim["time"]))
        found = sc.attach_sidecars(work, vhdr)
        stim.update({k:v for k, v in work.items() if k != "time"})
        if len(stim["time"]) and "sidecar" in stim:
            stim["confirmed"] = np.where(np.asarray(stim["sidecar"]) != "", 1, stim["confirmed"]).astype(np.int8)
    #what the measured width is: the matched sidecar's pulse width or its ttl code
    lag = stim["marker_lag_ms"]
    for col, label, scale in (("pulse_width_us", "pulse_width", 1e3), ("ttl_ms", "ttl_code", 1)):
        if col in stim and len(lag):
            ref = np.asarray(stim[col], dtype=np.float64) / scale
            ok = np.isfinite(ref) & np.isfinite(lag)
            if ok.sum() >= 3 and np.mean(np.abs(lag[ok] - ref[ok]) <= 0.025) >= 0.9:
                stim["marker_lag_meaning"] = np.full(len(lag), label, dtype=object)
    return {"stim":stim, "sidecars":found}


def code_table(cohort=None):
    """
    Input: cohort (e.g. "UG") or None for all
    Returns: list of (cohort, code, kind, type, meaning, evidence, confirmed) rows from RULES, for docs and checks
    """
    return [r[:7] for r in RULES if cohort is None or r[0] == cohort.upper()]

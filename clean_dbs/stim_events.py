from pathlib import Path
import numpy as np

#Synapse eStim gizmo settings, 6 per voice (A, B, C, D) = the 24 rows of eS1p, in this order
#read from the lab's Synapse experiments (WalkerStim, uh3_stimRec_dual): par1-24 = PerA, CntA, AmpA, DurA, DelA, ChnA, PerB, ...
ESTIM_PARAMS = ["period_ms", "count", "amplitude_ua", "duration_ms", "delay_ms", "channel"]
VOICES = "ABCD"
#the gizmo can also save any setting as its own epoc store, named by the user (AmpA or AmA/, CntA or CnA/, ...), so match the front of the name
EPOC_PREFIX = {"Per":"period_ms", "Pe":"period_ms", "Cnt":"count", "Cn":"count", "Amp":"amplitude_ua", "Am":"amplitude_ua",
               "Dur":"duration_ms", "Du":"duration_ms", "Del":"delay_ms", "De":"delay_ms", "Chn":"channel", "Ch":"channel"}


def stim_events(result, block):
    """
    One row per stim pulse, the same table for tdt and brainvision blocks.
    Input:
    result: the dictionary from read_tdt / read_bv
    block: the raw block folder or .vhdr it came from (only used to tell tdt from brainvision)
    Returns: dictionary of equal-length arrays, one entry per pulse (all empty if the block has no stim):
      time (s on the recording's time axis), source (store / marker name, e.g. PeA_, PerB or R15), voice (A-D for tdt, "" for brainvision),
      channel (stimulator channel, -1 = unknown), amplitude_ua (nan = unknown), period_ms (time to the next pulse in its train, nan = unknown),
      count (pulses in its train, -1 = unknown), confirmed (1 = tdt pulse store, 0 = brainvision marker whose meaning is inferred)
    Brainvision: one row per rising edge of a code that markers.py calls stim for this cohort (see markers.label_bv), confirmed = 1 when a rule
    or a stimulator file checked it. The falling marker of the same TTL (S13 after S15, ...) is never a second pulse; its delay is marker_lag_ms.
    Codes held 40 ms or more (buzz, task cues) are never stim. Sidecar columns (amplitude_ua, pulse_width_us, protocol_*, ...) are added when
    a stimulator file lines up with the markers.
    A tdt store whose onsets and values don't line up is trimmed to the shorter one and noted in result["info"]["stim_warning"].
    Heads up: tdt pulses with amplitude 0 are kept (the gizmo logs them but no current goes out), filter on amplitude_ua if needed.
    """
    if Path(block).suffix.lower() == ".vhdr":
        from clean_dbs.markers import label_bv
        return label_bv(result, block)["stim"]
    warnings = []
    out = _tdt_events(result["epocs"], result.get("scalars", {}), warnings)
    if warnings: #kept with the file (info stim_warning) so nobody has to trust a silent fix
        result.setdefault("info", {})["stim_warning"] = "; ".join(warnings)
    return out


def _tdt_events(epocs, scalars, warnings):
    rows = []
    for voice in VOICES:
        #PeA/ = one event per pulse on voice A, data = its period in ms (checked against MonA on uh3003); uh3001 calls voice B's PerB
        source = next((name for name in (f"Pe{voice}_", f"Per{voice}") if name in epocs), None)
        pulses = epocs.get(source)
        if pulses is None or len(pulses["onset"]) == 0:
            continue
        t, period = _same_length(source, pulses, warnings)
        info = {"period_ms":period}
        #eS1p = all 24 settings once per train, a few us before its first pulse; period (row 0) lags one train behind, so it isn't used
        if "eS1p" in scalars and len(scalars["eS1p"]["ts"]) > 0:
            ts, data = scalars["eS1p"]["ts"], scalars["eS1p"]["data"]
            train = np.searchsorted(ts, t + 1e-6, side="right") - 1 #+1 us so a setting saved at the same instant still counts
            for j in (1, 2, 5): #count, amplitude, channel
                values = data[6 * VOICES.index(voice) + j][np.clip(train, 0, None)]
                info[ESTIM_PARAMS[j]] = np.where(train >= 0, values, np.nan)
        #settings saved as their own epoc stores win over eS1p, they are logged at (or just before) the pulse itself
        for name, e in epocs.items():
            key = EPOC_PREFIX.get(name.rstrip("_")[:-1])
            if name.rstrip("_")[-1:] != voice or key in (None, "period_ms") or len(e["onset"]) == 0:
                continue
            onset, data = _same_length(name, e, warnings)
            last = np.searchsorted(onset, t + 1e-6, side="right") - 1
            info[key] = np.where(last >= 0, data[np.clip(last, 0, None)], np.nan) if len(data) else np.full(len(t), np.nan)
        rows.append(_table(t, source, voice, info, confirmed=1))
    return _stack(rows)


def _same_length(name, e, warnings):
    """
    Input: epoc name, the epoc (onset, data), list to add a warning to
    Returns: (onset, data) as float arrays of the same length. read_tdt already fixes the known case (uh3001/uh3006 CntA, one extra onset);
    if a store still has more onsets than values (or fewer) the extra ones at the end are dropped and it's written down, never silently
    """
    onset, data = np.asarray(e["onset"], dtype=np.float64), np.asarray(e["data"], dtype=np.float64)
    if len(onset) != len(data):
        n = min(len(onset), len(data))
        warnings.append(f"{name}: {len(onset)} onsets but {len(data)} values, used the first {n}")
        onset, data = onset[:n], data[:n]
    return onset, data


def _table(t, source, voice, info, confirmed):
    n = len(t)
    nan = np.full(n, np.nan)
    return {"time":t, "source":np.full(n, source, dtype="<U32"), "voice":np.full(n, voice, dtype="<U1"),
            "channel":np.nan_to_num(info.get("channel", nan), nan=-1).astype(np.int64),
            "amplitude_ua":info.get("amplitude_ua", nan), "period_ms":info.get("period_ms", nan),
            "count":np.nan_to_num(info.get("count", nan), nan=-1).astype(np.int64), "confirmed":np.full(n, confirmed, dtype=np.int8)}


def _stack(rows):
    #join the voices / markers and sort by time; empty arrays (right types) when there is no stim
    if not rows:
        rows = [_table(np.array([]), "", "", {}, 0)]
    out = {key:np.concatenate([r[key] for r in rows]) for key in rows[0]}
    order = np.argsort(out["time"], kind="stable")
    return {key:values[order] for key, values in out.items()}


def bv_clock_time(times, new_segments):
    """
    Brainvision glues the data together after a real recording gap (see flags.new_segment_gaps), so after each gap
    the recording's time axis falls behind the wall clock. This adds the missing time back.
    Input:
    times: times in s on the recording's time axis (e.g. stim_events(...)["time"] or epoc onsets)
    new_segments: the dictionary from new_segment_gaps (time, missing, is_gap)
    Returns: times in s on a clock-consistent axis (same 0 as the recording), same shape as times
    """
    times = np.asarray(times, dtype=np.float64)
    gap = np.asarray(new_segments["is_gap"]) == 1 #only real gaps, the tiny clock jitter at other NewSegments is ignored
    gap_time = np.asarray(new_segments["time"], dtype=np.float64)[gap]
    missing = np.cumsum(np.asarray(new_segments["missing"], dtype=np.float64)[gap])
    after = np.searchsorted(gap_time, times, side="right") #how many gaps came at or before each time
    return times + np.concatenate(([0.0], missing))[after]

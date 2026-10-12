from datetime import datetime
from pathlib import Path
import numpy as np


def find_saturation(stream, min_ms=1.0):
    """
    Finds stretches where a channel sits flat at its own maximum or minimum, the signature of the amplifier hitting its ceiling
    (e.g. some ecos channels pinned at 0.5369 V for ~2 s on uh3021 clinical_dbs1). Flags only, the data isn't changed.
    Input:
    stream: one entry from result["streams"] (needs data, fs, start_time, channels)
    min_ms: shortest flat stretch that counts, in milliseconds
    Returns: dictionary of arrays, one entry per stretch: channel, start, end (seconds from the start of the recording)
    Heads up: streams stored as whole mV (some older blocks) are flat a lot just from rounding, so treat flags there as a hint.
    """
    data, fs = stream["data"], float(stream["fs"])
    min_samples = max(2, int(round(min_ms / 1000 * fs)))
    channel, start, end = [], [], []
    for row, number in enumerate(stream["channels"]):
        x = data[row]
        if x.size == 0 or x.max() == x.min(): #a channel flat the whole time (unused, unplugged) isnt saturation, qc_block flags those
            continue
        for extreme in (x.max(), x.min()):
            hits = np.flatnonzero(x == extreme)
            if hits.size < min_samples:
                continue
            #split the hits into runs of back-to-back samples
            breaks = np.flatnonzero(np.diff(hits) > 1)
            run_starts = np.concatenate(([hits[0]], hits[breaks + 1]))
            run_ends = np.concatenate((hits[breaks], [hits[-1]]))
            keep = (run_ends - run_starts + 1) >= min_samples
            channel += [int(number)] * int(keep.sum())
            start += list(stream["start_time"] + run_starts[keep] / fs)
            end += list(stream["start_time"] + (run_ends[keep] + 1) / fs)
    return {"channel":np.asarray(channel, dtype=np.int64), "start":np.asarray(start, dtype=np.float64), "end":np.asarray(end, dtype=np.float64)}


def make_flags(result, block):
    """
    Input:
    result: the dictionary from read_tdt / read_bv
    block: the raw block folder or .vhdr it came from
    Returns: dictionary for result["flags"]:
      saturation_<stream>: channel, start, end of flat-at-the-ceiling stretches
      new_segments (brainvision only): time, missing, is_gap for every NewSegment marker, see new_segment_gaps
      dnu: 1 if the block name says do not use
    """
    flags = {}
    for name, stream in result["streams"].items():
        flags[f"saturation_{name}"] = find_saturation(stream)
    if Path(block).suffix.lower() == ".vhdr":
        flags["new_segments"] = new_segment_gaps(block)
    flags["dnu"] = np.int8("dnu" in Path(block).name.lower()) #Path so it works with Windows backslashes too
    return flags


def new_segment_gaps(vhdr_path, min_gap=0.001):
    """
    BrainVision writes a NewSegment marker with a wall-clock timestamp. Comparing the clock with the sample count tells
    whether recording actually stopped there. Checked on the inventory (2026-10-11): most NewSegments are NOT gaps
    (thousands of few-ms segments in the ERP/stim recordings, clock and samples agree), but a few per recording are
    real gaps of up to ~0.65 s that MNE silently glues together, so the time axis falls behind the real clock after each.
    Input:
    vhdr_path: the .vhdr (its .vmrk is found from the MarkerFile= line)
    min_gap: missing time (s) needed to count as a real gap; smaller differences are clock jitter
    Heads up: every non-gap NewSegment shows exactly one sample of missing time, either brainvision drops one sample
    per segment or the stamp is one sample late; we couldnt tell which, so only gaps over min_gap are treated as real
    Returns: dictionary of arrays, one entry per NewSegment after the first:
      time (s, on the recording's own time axis), missing (s of clock time with no samples before it), is_gap (1/0)
    """
    vhdr_path = Path(vhdr_path)
    header = vhdr_path.read_text(errors="ignore")
    marker_file = next((line.split("=", 1)[1].strip() for line in header.splitlines() if line.startswith("MarkerFile=")), vhdr_path.stem + ".vmrk")
    interval = next((float(line.split("=", 1)[1]) for line in header.splitlines() if line.startswith("SamplingInterval=")), None)
    empty = {"time":np.array([]), "missing":np.array([]), "is_gap":np.array([], dtype=np.int8)}
    vmrk = vhdr_path.parent / marker_file
    if interval is None or not vmrk.is_file():
        return empty
    dt = interval / 1e6 #SamplingInterval is in microseconds
    positions, clocks = [], []
    for line in vmrk.read_text(errors="ignore").splitlines():
        if line.startswith("Mk") and "New Segment" in line:
            fields = line.split("=", 1)[1].split(",")
            stamp = fields[5].strip() if len(fields) > 5 else ""
            if len(stamp) >= 20: #YYYYMMDDhhmmss + microseconds
                positions.append(int(fields[2]))
                clocks.append(datetime.strptime(stamp[:20], "%Y%m%d%H%M%S%f"))
    if len(positions) < 2:
        return empty
    positions = np.asarray(positions)
    clock_s = np.array([(c - clocks[0]).total_seconds() for c in clocks])
    missing = np.diff(clock_s) - np.diff(positions) * dt
    return {"time":(positions[1:] - 1) * dt, "missing":missing, "is_gap":(missing > min_gap).astype(np.int8)}

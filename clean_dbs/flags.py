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
        if x.size == 0:
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
      new_segments: times of brainvision NewSegment markers (usually a pause/restart, but erp001 has 334 of them
        some only 23 ms apart, so check what they mean on your rig before treating them as gaps)
      dnu: 1 if the block name says do not use
    """
    flags = {}
    for name, stream in result["streams"].items():
        flags[f"saturation_{name}"] = find_saturation(stream)
    new_segments = result["epocs"].get("NewSegment", {}).get("onset", np.array([]))
    flags["new_segments"] = {"time":np.asarray(new_segments, dtype=np.float64)[np.asarray(new_segments) > 0]} #the first NewSegment is just the start
    flags["dnu"] = np.int8("dnu" in Path(block).name.lower()) #Path so it works with Windows backslashes too
    return flags

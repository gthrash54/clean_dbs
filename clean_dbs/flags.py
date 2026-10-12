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
    A channel flat for the whole block counts only if it sits at the stream's rail (stream["rail"], from read_tdt); flat at 0 is just unused.
    """
    runs = Saturation(stream, min_ms)
    runs.add(stream["data"])
    return runs.result()


class Saturation:
    """
    find_saturation one time window at a time, so a recording never has to be in memory all at once (clean_block uses it).
    Input:
    stream: one entry from result["streams"] (needs fs, start_time, channels; data isn't used; rail optional)
    min_ms: shortest flat stretch that counts, in milliseconds
    Use: add(window) for each channels x samples window in time order, then result() gives what find_saturation gives
    Per channel it keeps the max and min so far and the runs sitting at them; a window with a new max (or min) drops the old runs,
    and a run still going at the end of a window is joined to one starting on the next window's first sample.
    """
    def __init__(self, stream, min_ms=1.0):
        self.fs, self.start_time, self.channels = float(stream["fs"]), stream["start_time"], stream["channels"]
        self.rail = float(stream.get("rail", np.nan)) #amplifier ceiling in volts (tdt), NaN = unknown (brainvision)
        self.min_samples = max(2, int(round(min_ms / 1000 * self.fs)))
        n = len(self.channels)
        self.n = 0 #samples seen so far
        self.has_nan = np.zeros(n, dtype=bool) #a NaN anywhere makes max/min NaN, so the whole-array version never flags that channel
        self.extreme = [[None] * n, [None] * n] #[max per channel, min per channel] so far
        self.runs = [[[] for _ in range(n)], [[] for _ in range(n)]] #[start, stop) sample pairs at that extreme

    def add(self, x):
        if x.shape[1] == 0:
            return
        for k, value in enumerate((x.max(axis=1), x.min(axis=1))):
            for row in range(x.shape[0]):
                v = value[row]
                if np.isnan(v):
                    self.has_nan[row] = True
                    continue
                best = self.extreme[k][row]
                if best is None or (v > best if k == 0 else v < best):
                    self.extreme[k][row], self.runs[k][row] = v, []
                elif v != best:
                    continue
                hits = np.flatnonzero(x[row] == v)
                #split the hits into runs of back-to-back samples
                breaks = np.flatnonzero(np.diff(hits) > 1)
                starts = np.concatenate(([hits[0]], hits[breaks + 1])) + self.n
                stops = np.concatenate((hits[breaks], [hits[-1]])) + 1 + self.n
                runs = self.runs[k][row]
                if runs and runs[-1][1] == starts[0]: #the last run reached the end of the previous window, so it carries on here
                    runs[-1][1] = int(stops[0])
                    starts, stops = starts[1:], stops[1:]
                runs += [[int(a), int(b)] for a, b in zip(starts, stops)]
                #drop short runs now (whole-mV streams can sit at their max thousands of times), but not one still going at the window's end
                self.runs[k][row] = [r for r in runs if r[1] - r[0] >= self.min_samples or r[1] == self.n + x.shape[1]]
        self.n += x.shape[1]

    def result(self):
        channel, start, end = [], [], []
        for row, number in enumerate(self.channels):
            if self.n == 0 or self.has_nan[row]:
                continue
            #a channel flat the whole time is saturated if it sits at the rail (uh3002 dbs1_clinicaltesting: 7 of 8 channels the whole block);
            #anywhere else (0, unused, unplugged) it isnt, qc_block flags those
            if self.extreme[0][row] == self.extreme[1][row]:
                if abs(self.extreme[0][row]) >= self.rail * (1 - 1e-3): #False when rail is NaN
                    channel.append(int(number))
                    start.append(self.start_time)
                    end.append(self.start_time + self.n / self.fs)
                continue
            for k in (0, 1): #max then min, same order as before
                runs = np.asarray(self.runs[k][row], dtype=np.int64).reshape(-1, 2)
                keep = (runs[:, 1] - runs[:, 0]) >= self.min_samples
                channel += [int(number)] * int(keep.sum())
                start += list(self.start_time + runs[keep, 0] / self.fs)
                end += list(self.start_time + runs[keep, 1] / self.fs)
        return {"channel":np.asarray(channel, dtype=np.int64), "start":np.asarray(start, dtype=np.float64), "end":np.asarray(end, dtype=np.float64)}


def make_flags(result, block, saturation=None):
    """
    Input:
    result: the dictionary from read_tdt / read_bv
    block: the raw block folder or .vhdr it came from
    saturation: optional {stream name: find_saturation result} already worked out window by window (clean_block);
                left out, it's found from each stream's data
    Returns: dictionary for result["flags"]:
      saturation_<stream>: channel, start, end of flat-at-the-ceiling stretches
      new_segments (brainvision only): time, missing, is_gap for every NewSegment marker, see new_segment_gaps
      dnu: 1 if the block name says do not use
    """
    flags = {}
    for name, stream in result["streams"].items():
        flags[f"saturation_{name}"] = saturation[name] if saturation is not None else find_saturation(stream)
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

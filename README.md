# clean_dbs

This is the Walker lab toolbox for loading and cleaning our intraoperative DBS recordings. Right now it reads both of our rigs: stage 1 surgery (TDT) and stage 2 surgery (BrainVision). It currently reads tdt and bv files but is growing to become a cleaning and preprocessing pipeline for all current and future recordings.

Both readers give back the same thing, so whatever you write on top of them works no matter which rig the data came from. It runs on Windows, Mac and Linux, and you can call it from Python or MATLAB.

If something breaks or doesn't make sense, message me (Garrett Thrash, garrettthrash54@gmail.com).

## Setup (one time)

You need two things: a copy of this `clean_dbs` folder, and [uv](https://docs.astral.sh/uv/getting-started/installation/), which installs Python and the packages for you.

Open a terminal inside the `clean_dbs` folder (the outer one, with this README in it) and run:

```
uv venv --python 3.12
uv pip install tdt numpy mne
```

That makes a private Python for the toolbox in a folder called `.venv`. You only do this once.

Every time you open a new terminal to use the toolbox, switch that Python on first:

- Mac / Linux: `source .venv/bin/activate`
- Windows: `.venv\Scripts\activate`

## Reading a recording in Python

Run Python from the outer `clean_dbs` folder, then:

```python
from clean_dbs.read_tdt import read_tdt
from clean_dbs.read_bv import read_bv

# TDT (stage 1): point at ONE block folder, not the whole "stage1 surgery" folder
r = read_tdt("<path to>/stage1 surgery/clinical_dbs1", t1=0, t2=10)

# BrainVision (stage 2): point at the .vhdr file, the .eeg and .vmrk have to sit next to it
r = read_bv("<path to>/stage 2 surgery/STN_dbs1.vhdr", t1=0, t2=10)
```

`t1` and `t2` are the start and end in seconds. Leave them out (or use 0) to read the whole recording. Full recordings are big (around 2 GB), so use a window when you're just looking around.

## What you get back

A dictionary with three parts. It's the same layout for both rigs.

```
r["streams"]   one entry per signal (TDT: "ecos", "emgg", "mic_", ...; BrainVision: "dbs")
    ["data"]        channels x samples, in VOLTS, always 2-D (even one channel is 1 x samples)
    ["fs"]          exact sampling rate in Hz (e.g. 48828.125, not rounded)
    ["start_time"]  time of the first sample, in seconds from the start of the recording
    ["channels"]    channel numbers

r["epocs"]     one entry per marker type (TDT: "Cam1", "PeA_", ...; BrainVision: "R15", "S13", ...)
    ["onset"]       when each marker starts, in seconds from the start of the recording
    ["offset"]      when it ends
    ["data"]        TDT only for now: the value tied to each marker (e.g. frame number for Cam1)

r["info"]      block name, start date and time, duration in seconds
               (TDT also gives subject and the rig template under "experiment")
```

To get the time of every sample in a stream: `start_time + np.arange(n_samples) / fs`. Don't assume the first sample is exactly at `t1`; on TDT it's a few microseconds after.

## Using it from MATLAB

MATLAB has to know which Python to use. Run this once per MATLAB session, with your own paths:

```matlab
pyenv(Version="<path to clean_dbs>/.venv/bin/python", ExecutionMode="OutOfProcess");
% Windows: Version="<path to clean_dbs>\.venv\Scripts\python.exe"
P = py.sys.path; P.insert(int32(0), "<path to clean_dbs>");
```

Then:

```matlab
m = py.importlib.import_module("clean_dbs.read_tdt");
r = struct(m.read_tdt("<path to block>", t2=10));
streams = struct(r.streams);
ecos = struct(streams.ecos);
x = single(ecos.data);      % 28 x samples, volts
fs = ecos.fs;
```

MATLAB only converts one level at a time, so you need a `struct()` at each level for now. I'm writing a small MATLAB helper that does it all in one call. Use `single()` for the data, not `double()`; `double()` doubles the memory.

## Things that will trip you up

- **Wrong folder for TDT.** Point at a single block (like `clinical_dbs1`), not `stage1 surgery`. It'll tell you if you get this wrong.
- **Incomplete copies.** A lot of the local copies only have the small index files and are missing the actual signal (`.tev` for TDT, `.eeg` for BrainVision). Box has the full recordings. It'll tell you which file is missing.
- **Everything is in volts.** BrainVision files store microvolts, but the reader converts them, so both rigs match. Multiply by 1e6 when you want µV for a plot.
- **Marker names depend on the rig.** On the recording I tested (uh3 cohort, stage 2), R15 looks like the start of each stim pulse and R11 comes 60 µs later. Other patients use S13/S15 and so on. Check what they mean for your rig before you rely on them.
- **NewSegment markers mean gaps.** If a BrainVision recording has `NewSegment` markers, the recording was paused and restarted there. Don't filter across those points without handling the gap.
- **Some channels saturate.** Some TDT blocks hit the amplifier ceiling (about 0.537 V) for a few seconds. 

## What's in here

```
clean_dbs/
    read_tdt.py    reads a TDT block (stage 1)
    read_bv.py     reads a BrainVision recording (stage 2)
scripts/           one-off scripts
examples/          examples (coming)
```

Recordings and anything made from them never go in this folder or in git. The `.gitignore` blocks the usual data files, but double check before you commit.

# clean_dbs

This is the Walker lab toolbox for loading and cleaning our intraoperative DBS recordings. Right now it reads both of our rigs: stage 1 surgery (TDT) and stage 2 surgery (BrainVision). It currently reads tdt and bv files but is growing to become a cleaning and preprocessing pipeline for all current and future recordings.

Both readers give back the same thing, so whatever you write on top of them works no matter which rig the data came from. It runs on Windows, Mac and Linux, and you can call it from Python or MATLAB.

If something breaks or doesn't make sense, message me (Garrett Thrash, garrettthrash54@gmail.com).

## Setup (one time)

You need two things: a copy of this `clean_dbs` folder, and [uv](https://docs.astral.sh/uv/getting-started/installation/), which installs Python and the packages for you.

Open a terminal inside the `clean_dbs` folder (the outer one, with this README in it) and run:

```
uv venv --python 3.12
uv pip install -e .
```

That makes a private Python for the toolbox in a folder called `.venv` and installs the toolbox into it, along with what it needs (tdt, numpy, mne, h5py). The `-e` means edits to the code take effect right away, no reinstall. You only do this once.

Every time you open a new terminal to use the toolbox, switch that Python on first:

- Mac / Linux: `source .venv/bin/activate`
- Windows: `.venv\Scripts\activate`

## Reading a recording in Python

With the `.venv` switched on, this works from any folder:

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
    ["device"]      which device recorded it (e.g. ecogdbs, emg, microphone, brainvision)
    ["units"]       always "V"

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

MATLAB only converts one level at a time, so you need a `struct()` at each level. Use `single()` for the data, not `double()`; `double()` doubles the memory. If you just want cleaned data in MATLAB, skip Python entirely and use `clean_dbs_load` (below).

## Cleaning a whole folder

```
python scripts/clean.py "<raw folder>" "<cleaned folder>"
```

The raw folder can be the whole inventory, one patient, or one surgery folder. Each recording becomes `<cleaned folder>/<patient>/<stage>/<block>.h5` with:

- the signals in volts, tagged by type (lfp, mer, emg, audio, stim_monitor) with channel names where we know them
- the events (Cam1 video frames, stim pulses, markers) and the stimulator settings per pulse (eS1p)
- flags: saturated stretches, real gaps in BrainVision recordings, impedance checks, dnu
- where the raw data and the video are, and which version of this code made the file

Blocks that are already cleaned get skipped, so adding a new patient is the same command. `clean_report.csv` in the cleaned folder lists every block as cleaned, skipped or failed with the reason. To run the whole inventory on Cheaha, see `docs/cheaha.md`.

To find things across the cleaned set:

```python
from clean_dbs import build_index, find_clean, qc_report
idx = build_index("<cleaned folder>")          # one row per block
files = find_clean(idx, type="lfp", target="STN", implant="dual", max_saturated_s=0)
```

Case info (dual or single, hemisphere, targets, depths, stim contacts) comes from `<cleaned folder>/_cases/cases.csv` and `blocks.csv`, transcribed from the OR data sheets. Only rows someone has checked against the sheet (`verified` filled in) get used.

## Loading a cleaned file in MATLAB (no Python needed)

```matlab
addpath("<path to clean_dbs>/matlab")
r = clean_dbs_load("<cleaned folder>/uh3021/stage1/clinical_dbs1.h5", "types", {"lfp"}, "t1", 10, "t2", 20);
x = r.streams.ecos.data;   % channels x samples, single, volts
```

Only the window you ask for is read from disk. Flags (like `r.flags.saturation_ecos`) cover the whole recording.

## Running the tests

```
uv pip install -e ".[dev]"
python -m pytest tests
```

They use small made-up files, so no recordings are needed.

## What's in the cleaned folder

- `<patient>/<stage>/<block>.h5`: one cleaned recording
- `_cases/cases.csv`, `_cases/blocks.csv`: patient and recording info from the OR sheets (only rows with `verified` filled in get used)
- `clean_report.csv`: what the last run did with every recording
- `_about/`: the user guide PDF

## Things that will trip you up

- **Wrong folder for TDT.** Point at a single block (like `clinical_dbs1`), not `stage1 surgery`. It'll tell you if you get this wrong.
- **Incomplete copies.** A lot of the local copies only have the small index files and are missing the actual signal (`.tev` for TDT, `.eeg` for BrainVision). Box has the full recordings. It'll tell you which file is missing.
- **Everything is in volts.** BrainVision files store microvolts, and some older TDT blocks (uh3002, uh3003) store whole numbers in µV or mV. The readers convert all of it, so every rig matches. Multiply by 1e6 when you want µV for a plot.
- **Some older microelectrode channels are basically unusable.** On uh3002 the spike and micro/macro LFP stores were saved as whole millivolts, so they only take a few flat values. That's how they were recorded, the code can't get the detail back.
- **Marker names depend on the rig.** On the recording I tested (uh3 cohort, stage 2), R15 looks like the start of each stim pulse and R11 comes 60 µs later. Other patients use S13/S15 and so on. Check what they mean for your rig before you rely on them.
- **Check NewSegment markers.** BrainVision writes a `NewSegment` marker when a recording is paused and restarted, but one recording I checked has 334 of them, some only 23 ms apart. They're saved under `/flags/new_segments` in the cleaned file. Figure out what they mean for your rig before you filter across them.
- **Some channels saturate.** Some TDT blocks hit the amplifier ceiling (about 0.537 V) for a few seconds. Cleaning flags these under `/flags/saturation_<stream>` (channel, start, end) but doesn't change the data. 

## What's in here

```
clean_dbs/
    read_tdt.py        reads a TDT block (stage 1)
    read_bv.py         reads a BrainVision recording (stage 2)
    find_blocks.py     finds every recording under any folder
    check_blocks.py    table of what's there and what's missing, before cleaning
    clean_block.py     reads, tags and saves one recording
    save_h5.py         writes the cleaned .h5
    flags.py           saturation, BrainVision gaps, dnu
    channel_map.py     channel names and impedance checks
    datatypes.py       which streams are lfp, mer, emg, ...
    block_info.py      patient, stage and output path for a recording
    index.py           table of the cleaned set, find_clean, quick QC
scripts/
    clean.py           cleans a folder (or one Cheaha array task)
    make_manifest.py   fixed list of recordings for Cheaha
    clean_slurm.sh     the Cheaha job
matlab/
    clean_dbs_load.m   loads a cleaned file into a struct
tests/                 small tests, no recordings needed
docs/cheaha.md         running the whole inventory on Cheaha
examples/              examples (coming)
```

Recordings and anything made from them never go in this folder or in git. The `.gitignore` blocks the usual data files, but double check before you commit.

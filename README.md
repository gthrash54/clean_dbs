# clean_dbs

This is the Walker lab toolbox for loading and cleaning our intraoperative DBS recordings. Right now it reads both of our rigs: stage 1 surgery (TDT) and stage 2 surgery (BrainVision). It currently reads tdt and bv files but is growing to become a cleaning and preprocessing pipeline for all current and future recordings.

Both readers give back the same thing, so whatever you write on top of them works no matter which rig the data came from. It runs on Windows, Mac and Linux, and you can call it from Python or MATLAB.

If something breaks or doesn't make sense, message me (Garrett Thrash, garrettthrash54@gmail.com).

## Quick start

You can try it without any recordings. In a fresh virtual environment:

```
pip install git+https://github.com/gthrash54/clean_dbs
```

Then run this from an empty folder. It makes a small fake BrainVision recording in `demo_raw/`, cleans it into `demo_clean/`, and averages the signal around each stim pulse:

```python
from pathlib import Path
import h5py
import numpy as np
from clean_dbs import clean_block, load_type

# 1. A fake BrainVision recording: 8 lead contacts (DBS1-8), 10 s at 10 kHz, in microvolts.
#    A stim pulse every 50 ms, marked like the lab's stimulator sync (S 15 up, S 13 down one sample later),
#    each followed by a small 300 Hz response in 5 uV of noise.
fs, n = 10_000, 100_000
pulses = np.arange(5_000, n - 1_000, 500)  # sample of each pulse
t = np.arange(150) / fs
response = 20 * np.sin(2 * np.pi * 300 * t) * np.exp(-t / 0.004)
data = np.random.default_rng(0).normal(0, 5, (8, n))
for p in pulses:
    data[:, p + 10:p + 160] += response  # starts 1 ms after the pulse

folder = Path("demo_raw/demo/stage2 surgery")  # <patient>/<stage> folders, like the lab's raw data
folder.mkdir(parents=True, exist_ok=True)
(folder / "rec.vhdr").write_text(
    "Brain Vision Data Exchange Header File Version 1.0\n[Common Infos]\nCodepage=UTF-8\n"
    "DataFile=rec.eeg\nMarkerFile=rec.vmrk\nDataFormat=BINARY\nDataOrientation=MULTIPLEXED\n"
    f"NumberOfChannels=8\nSamplingInterval={1e6 / fs:g}\n[Binary Infos]\nBinaryFormat=IEEE_FLOAT_32\n"
    "[Channel Infos]\n" + "".join(f"Ch{k}=DBS{k},,1,µV\n" for k in range(1, 9)), encoding="utf-8")
data.T.astype("<f4").tofile(folder / "rec.eeg")
markers = "".join(f"Mk{2 * i + 2}=Stimulus,S 15,{p + 1},1,0\nMk{2 * i + 3}=Stimulus,S 13,{p + 2},1,0\n"
                  for i, p in enumerate(pulses))
(folder / "rec.vmrk").write_text(
    "Brain Vision Data Exchange Marker File, Version 1.0\n[Common Infos]\nCodepage=UTF-8\nDataFile=rec.eeg\n"
    "[Marker Infos]\nMk1=New Segment,,1,1,0,20240101120000000000\n" + markers, encoding="utf-8")

# 2. Clean it: one .h5 with the channels labeled, the data in volts, and a stim table made from the markers
h5 = clean_block(folder / "rec.vhdr", "demo_clean")
print(h5)

# 3. Load it back and average contact 1 over the 20 ms after each stim pulse
lead = load_type(h5, "lfp")["lead"]
print(lead["channel_names"], lead["channel_group"][0], lead["units"])
with h5py.File(h5, "r") as f:
    stim = f["stim/time"][:]  # seconds
onsets = np.round((stim - lead["start_time"]) * lead["fs"]).astype(int)
avg = np.mean([lead["data"][0, k:k + 200] for k in onsets], axis=0) * 1e6  # uV
peak = np.abs(avg).argmax()
print(f"{len(stim)} pulses; the average peaks at {avg[peak]:.1f} uV, {1e3 * peak / lead['fs']:.1f} ms after the pulse")
```

It prints:

```
demo_clean/demo/stage2/rec.h5
['dbs1', 'dbs2', 'dbs3', 'dbs4', 'dbs5', 'dbs6', 'dbs7', 'dbs8'] lead1 V
188 pulses; the average peaks at 17.2 uV, 1.8 ms after the pulse
```

`clean_block` read DBS1-8 as lead 1, contacts 1-8, converted microvolts to volts, and called `S 15` a stim pulse because its line stayed high only for a pulse width (see `markers.py`; `confirmed` is 0 in `/stim`, since `demo` isn't a cohort with a checked marker rule). The fake response peaks at 16.3 uV, 1.8 ms after each pulse, so the average finds it at the right time, with a little noise left over.

`clean_block` cleans one recording. To clean a whole folder of real recordings, use `scripts/clean.py` (below), which needs the case tables. `pip install` gives you the package only; for `scripts/`, `matlab/` and the tests, get a copy of this repository and follow the setup below.

## Setup (one time)

You need two things: a copy of this `clean_dbs` folder, and [uv](https://docs.astral.sh/uv/getting-started/installation/), which installs Python and the packages for you.

Open a terminal inside the `clean_dbs` folder (the outer one, with this README in it) and run:

```
uv venv --python 3.12
uv pip install -e .
```

That makes a private Python for the toolbox in a folder called `.venv` and installs the toolbox into it, along with what it needs (tdt, numpy, mne, h5py). The `-e` means edits to the code take effect right away, no reinstall.

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
r["streams"]   one entry per signal (TDT: "ecos", "emgg", "mic_", ...; BrainVision: one per channel type, "lead", "ecog", "eeg", "emg", "audio", "mer", ...)
    ["data"]        channels x samples, in VOLTS, always 2-D (even one channel is 1 x samples)
    ["fs"]          exact sampling rate in Hz (e.g. 48828.125, not rounded)
    ["start_time"]  time of the first sample, in seconds from the start of the recording
    ["channels"]    channel numbers
    ["units"]       "V" (a BrainVision channel recorded in some other unit, like BIP, gets its own stream that keeps that unit)
    ["rail"]        the amplifier ceiling in volts (TDT 0.537 V or 32000 counts; BrainVision actiCHamp 0.4096 V), NaN when unknown
    BrainVision streams also say what each channel is: type, subtype, channel_names, channel_group (lead1, ecog_anterior, emg, ...),
                    channel_contact, channel_muscle, channel_side, channel_units, channel_source. read_bv(..., split=False) gives the
                    old single "dbs" stream with every channel instead

r["epocs"]     one entry per marker type (TDT: "Cam1", "PeA_", ...; BrainVision: "R15", "S13", "BtnPressed", ...)
    ["onset"]       when each marker starts, in seconds from the start of the recording
    ["offset"]      when it ends
    ["data"]        TDT only for now: the value tied to each marker (e.g. frame number for Cam1)

r["info"]      block name, start date and time (the rig's Birmingham clock time, whatever computer you run it on; BrainVision
               takes it from the .vmrk as written, both as "2025-01-31T09:15:02.123456"), duration in seconds (from the longest
               stream); TDT also gives subject, the rig template under "experiment", duration_header (stop minus start from the
               block header, NaN when the header has no stop time) and epoc_warning when tdt gave a store more onsets than values;
               BrainVision also gives cohort
```

To get the time of every sample in a stream: `start_time + np.arange(n_samples) / fs`. Don't assume the first sample is exactly at `t1`; on TDT it's a few microseconds after.

## Channel labels and contact locations

Every lead and ECoG channel in a cleaned file says what it is, right next to the data:

- `channel_names`, `channel_group` (ecog_anterior, ecog_posterior, lead1, lead2, emg, ...), `channel_target` (STN, GPi, from the case table), `channel_contact`
- contact numbering is the lab's 1-8 for every lead: 1 = bottom (distal) ring, 2/3/4 = lower segments, 5/6/7 = upper segments (5 sits above 2, 6 above 3, 7 above 4), 8 = top ring. On a Medtronic SenSight lead that's levels 0, 1a-c, 2a-c, 3. Which way each segment faces isn't recorded.
- ECoG: anterior strip = 1-6, posterior = 7-12 (from the OR sheets). TDT lead 1 = 13-20, lead 2 = 21-28, worked out from impedance and stim in the data (`channel_source` says so) because the rig wiring isn't written down anywhere. Lead 1 is the lead implanted and stimulated first; stage 2 BrainVision has the STN lead on DBS1-8 and GPi on DBS9-16 (lab collection guide).
- **Contact coordinates** (`channel_x`, `channel_y`, `channel_z` in mm, `channel_position`, and the stream's `coord_space`): where each lead contact was during that recording, from the lab's Brainlab reconstructions (`_cases/contacts.csv`). The space is Brainlab ACPC in the patient's own hemisphere (+x right, +y anterior, +z superior). Only the final lead position is imaged; inching blocks (`increment1_depth3`, ...) get that position slid along the lead to the depth they were recorded at. A lead that wasn't in the brain yet gets NaN. Only streams with lead channels get these fields, so an ecog or emg stream has none. So far only 18 uh patients have reconstructions. `coord_status` in the info says which case a file is: `recon` (the patient is in `contacts.csv`), `no_recon` (the table was there but the patient isn't in it) or `no_contacts_csv` (there was no `_cases/contacts.csv` when it was cleaned, so nothing was looked up).
- **Do not use** (`channel_dnu`, 1 = don't use it, and the stream's `dnu_list`): channels the lab marked do not use for that recording (bad or noisy contacts, contacts still in the cannula, bad ECoG or EMG), from `_cases/dnu_channels.csv`. `dnu_list` says where the list came from, or `none` when the recording isn't on it, so then nothing is known and every `channel_dnu` is 0. The data itself isn't changed. So far the list covers uh stage 1, from the lab's sheet.
- TDT microelectrode (MER) streams carry `analog_gain` = 4000: they come from the AlphaOmega analog output, so electrode volts = data / 4000. Which AO output fed which TDT input isn't documented, so check before comparing MER amplitudes to LFP.

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

The cleaned folder has to have the case tables in it first (`<cleaned folder>/_cases/cases.csv` and `blocks.csv`, see below): every file gets its lead targets and a `case_status` from them, so `clean.py` won't start without them. For a cohort with no case tables, add `--no-cases` (targets stay blank, `case_status` is `missing`).

For uh, put `_cases/contacts.csv` in too. That's where the contact coordinates come from. `clean.py` still runs without it, but every uh file comes out with no coordinates, so `clean_block` prints a warning once per run and the files say `coord_status` = `no_contacts_csv`.

Put `_cases/dnu_channels.csv` in as well (the lab's do-not-use list, see below). Once `cases.csv` is there, `clean.py` won't start without it, so a run can't quietly go out with every file saying "not listed". Add `--no-dnu` to clean without it on purpose. `--no-cases` skips it too, since a cohort with no case tables has no list.

The raw folder can be the whole inventory, one patient, or one surgery folder. Each uh/UG recording becomes `<cleaned folder>/<patient>/<stage>/<block>.h5`; recordings from the other BrainVision cohorts (DLPFC, BI, PP, MC, ...) become `<cleaned folder>/<study ID>/<session>/<name>.h5`, where the session is the subfolder path below the study ID (lowercase, spaces as `_`, nothing when the file sits right in the ID folder). Each file has:

- the signals in volts, tagged by type (lfp, eeg, mer, emg, audio, accel, joystick, stim_monitor, unknown) with channel names where we know them; BrainVision recordings are split into one stream per channel type (`lead`, `ecog`, `emg`, ...)
- per-channel labels: TDT ecog/lead channels get `ecog_a1 ... lead2_c8` with the case's targets; a lead that was never implanted (single-lead cases, uh3023's GPi) keeps its channels as `lead2_not_implanted`, and TDT blocks from a stage 2 folder (uh3005) mark the strip channels `ecog_*_not_placed`
- the events (Cam1 video frames, stim pulses, markers, each with a `type`: video, stim, stim_params, task, marker, other) and the stimulator settings per pulse (eS1p)
- `/stim`: one row per stim pulse (TDT pulse stores, or the BrainVision codes `markers.py` calls stim) with voice, channel, amplitude and period where known, and `/sidecars`: the stimulator and task files found next to a BrainVision recording
- flags: saturated stretches, real gaps in BrainVision recordings, impedance checks, dnu
- where the raw data and the video are, and which version of this code made the file (`clean_version`, the git commit, package versions)
- `case_status` in the info: the cross-check status of the patient's `cases.csv` row when it was cleaned (`match`, `conflict`, `sheet_only`, `unchecked` if never checked, `missing` if there's no row), so you know how far to trust the targets

Cleaning reads, checks and writes the recording one window at a time (`window_mb` in `clean_block`, 128 MB by default), so it needs about 0.5 GB of memory no matter how big the file is. The cleaned file is the same whatever the window size.

Blocks that are already cleaned get skipped, so adding a new patient is the same command. A cleaned file is only skipped if it opens and was made by the current `CLEAN_VERSION` (in `clean_dbs/version.py`, saved as `clean_version` at the top of every `.h5`); a file that won't open or comes from an older version is cleaned again. `--force` cleans everything again. A recording copied under more than one number folder on Box is cleaned once; the other copies show up in the report as `duplicate` (same size) or `conflict` (same name and place but a different size, the bigger one is cleaned; check these by hand).

**Editing a case table doesn't make a file stale.** A file only counts as out of date when `CLEAN_VERSION` goes up, so a fix to `cases.csv`, `blocks.csv` or `contacts.csv` only reaches the files cleaned after it. To get the fix into files that are already there, clean them again with `--force` (pointing at just that patient's raw folder is enough), or wait for the next version bump.

**The CLEAN_VERSION rule:** whenever you change the code in a way that changes what goes into a cleaned file (a new or renamed field, different values, a different flag), add 1 to `CLEAN_VERSION` in `clean_dbs/version.py` in the same commit, and add a line to the list below. Otherwise the next run skips every old file and the cleaned set ends up a mix of old and new. Changes that leave the `.h5` the same (speed, messages, docs) don't need a bump.

- `1`: first numbered version (truncated `.tev` check, Chicago start time, stream-based duration, `case_status`, saturation at the rail)
- `2`: BrainVision split into typed streams with per-channel labels, marker meanings and stim from `markers.py`, `/sidecars`; a `rail` on every stream (BrainVision actiCHamp too); BrainVision start time straight from the `.vmrk`; `cohort` in every file's info; `CnA_` typed as a setting, `PerB` read as voice B pulses, uneven TDT epocs fixed (`epoc_warning`); `lead<N>_not_implanted` and stage 2 strip labels
- `3`: lead contact coordinates per channel (`channel_x/y/z`, `channel_position`, `coord_space`) from `_cases/contacts.csv`; `analog_gain` on TDT MER streams
- `4`: fixes from a trial run on the other cohorts: better `condition`/`state` from names (a folder counts when the file name says nothing, a study ID in a file name isn't a condition, `impedance_check`, `off_meds`/`on_meds`), the same ID nested inside itself keeps its session; stimulator programs found by their header whatever they're called, stricter sidecar matching and matching per marker code, no condition when two TTL tables disagree, binary MC_Stimulus files listed under `/sidecars`; headers that give two channels the same name can be read
- `5`: inching blocks named without the underscore (uh3002's `increment1depth3`) get their own depth's contact coordinates instead of the final ones
- `6`: an implanted lead that was never plugged into the amplifier (uh3016's GPi lead) is labeled `lead<N>_not_connected` with no target, from a cases.csv note "<target> lead never connected"
- `7`: `channel_dnu` and `dnu_list` on every stream (the lab's do-not-use list, from `_cases/dnu_channels.csv`). Only streams with lead channels get contact coordinates (ecog, emg and the rest used to get all-NaN `channel_x/y/z` and a `coord_space`), `coord_status` in the info (`recon`, `no_recon`, `no_contacts_csv`), and "post op", "postop" and "post_op" folders give the state `postop`. The two-word states (pre op, post dbs, off meds, on meds) also match when written with a space or a dash

A TDT `.tev` that was cut short by a stopped copy fails with "cut short" instead of being read (tdt would fill the missing part with zeros without saying anything). `python scripts/check_copy.py "<raw folder>"` lists those, and recordings missing their signal file, before you clean.

Each cohort's layout and what its recordings seem to be (rest, AX-CPT, iTBS, paired pulses, awake/asleep, follow-up visit, ...) is written down in `clean_dbs/cohorts.py`. The `condition` and `state` that `describe_block` gives (and that go in the manifest and the cleaned file's info) come from folder and file names only, so treat them as a first guess. When a study ID folder sits inside another one (like `PP012/CHR003/`), the recording belongs to the inner ID and the outer one is kept as `parent_id`. `clean_report.csv` in the cleaned folder lists every block with its status and the reason: `cleaned`, `skipped` (already done by this version), `recleaned` (there but broken, older, or `--force`), `duplicate`, `conflict`, `failed` (paths in the reason start at the study ID), `excluded` (never cleaned: no study ID above it or in its name, or a folder or file name that isn't cleared, see below), or `not run` (array jobs only). It's saved after every block, and `clean.py` exits with an error code if any block failed. To run the whole inventory on Cheaha, see `docs/cheaha.md`.

Some folder and file names on Box have a person's first name in them, and names go into the cleaned folders, file names, the info in every file and the reports. So every name that would go into the cleaned set has to be made of words that describe recordings (`WORDS` in `clean_dbs/privacy.py`), or be renamed in `<cleaned folder>/_cases/renames.csv`, or have its words in `_cases/allowed_words.txt`. `make_manifest.py` stops otherwise and lists what it couldn't clear in `_cases/names_to_check.csv`, and `clean.py` checks every recording again before cleaning it. How to fill those in is in `docs/cheaha.md`, "Cleaning the other cohorts".

To find things across the cleaned set:

```python
from clean_dbs import build_index, find_clean, qc_report
idx = build_index("<cleaned folder>")          # one row per block
files = find_clean(idx, type="lfp", target="STN", implant="dual", max_saturated_s=0)
```

Case info (dual or single, hemisphere, targets, depths, stim contacts) comes from `<cleaned folder>/_cases/cases.csv` and `blocks.csv`, transcribed from the OR data sheets by AI, so any value can be a misreading. `scripts/crosscheck_cases.py` checks every row against the recordings and adds a `status` to each one:

- `match`: the data agrees (leads in place from the impedance checks, stim contacts from which lead contact goes quiet while it's on the stimulator, target from which lead was stimulated)
- `conflict`: it doesn't; `check_detail` says why, so check that row against the paper sheet
- `sheet_only`: nothing in the data to check it with
- `no_data`: the recording wasn't found

`build_index` carries the status through, so `find_clean(idx, target="STN", status="match")` only uses rows the data confirms. `verified` is still there for rows someone has checked by hand.

## Loading a cleaned file in MATLAB (no Python needed)

```matlab
addpath("<path to clean_dbs>/matlab")
r = clean_dbs_load("<cleaned folder>/uh3021/stage1/clinical_dbs1.h5", "types", {"lfp"}, "t1", 10, "t2", 20);
x = r.streams.ecos.data;   % channels x samples, single, volts
```

Only the window you ask for is read from disk. Flags (like `r.flags.saturation_ecos`) cover the whole recording. `r.stim` is the stim table (one entry per pulse inside the window, text columns as cell arrays), `r.info` has every info attribute (including `duration_header`, `case_status` and `coord_status`) and `r.provenance` has `clean_version`.

In Python, `load_type("<cleaned file>.h5", "lfp")` (from `clean_dbs`) gives every lfp stream in a file with its data, `fs`, channel labels (`channel_names`, `channel_group`, ...) and contact coordinates, the same per-channel fields `clean_dbs_load` gives.

## Running the tests

```
uv pip install -e ".[dev]"
python -m pytest tests
```

They use small made-up files, so no recordings are needed.

## What's in the cleaned folder

- `<patient>/<stage>/<block>.h5` (uh/UG) or `<study ID>/<session>/<name>.h5` (other cohorts): one cleaned recording
- `_cases/cases.csv`, `_cases/blocks.csv`: patient and recording info from the OR sheets, with the cross-check `status` of every row (needed before cleaning). blocks.csv has a `stage` column, because a few UG blocks have the same name in stage 1 and stage 2
- `_cases/contacts.csv`: where each lead contact is, from the Brainlab reconstructions (uh only so far). Without it the uh files get no contact coordinates
- `_cases/dnu_channels.csv`: the channels the lab marked do not use, one row per channel per recording (a blank channel means the recording was checked and nothing is bad). `scripts/make_dnu_table.py` builds it from the lab's `UH3 or notes - data to not use.xlsx` (I keep a copy in `_cases/lab_sheets/`): `python scripts/make_dnu_table.py "<sheet>.xlsx" <cleaned folder>/_cases`. Reading the .xlsx needs openpyxl (`pip install -e ".[sheets]"`), cleaning doesn't. It prints any sheet row that matches no recording
- `_cases/dnu_block_fixes.csv`: hand-checked fixes for sheet rows whose recording name doesn't match a real recording (typos, a note typed in the trial column, renamed recordings). Columns `study_id`, `stage`, `sheet_block`, `block` (the real recording, blank = drop the row, `*` = every recording of that patient and stage) and `decided` (who decided and why). `make_dnu_table.py` applies it, so the lab's sheet itself stays untouched
- `clean_report.csv`: what the last run did with every recording (cleaned, skipped, recleaned, duplicate, conflict, failed, excluded, not run); with a manifest named `manifest_<cohort>.csv` it's `clean_report_<cohort>.csv`
- `_cases/renames.csv`, `_cases/allowed_words.txt`, `_cases/names_to_check.csv`: the names check (see above). `renames.csv` and `names_to_check.csv` have the original names in them, so leave them (and the `manifest*.csv` files, which have the raw paths) out when you share the cleaned set
- `_about/`: the user guide PDF

## Things that will trip you up

- **Wrong folder for TDT.** Point at a single block (like `clinical_dbs1`), not `stage1 surgery`. It'll tell you if you get this wrong.
- **Everything is in volts.** BrainVision files store microvolts, and some older TDT blocks (uh3002, uh3003) store whole numbers in µV or mV. The readers convert all of it, so every rig matches. Multiply by 1e6 when you want µV for a plot.
- **Some older microelectrode channels are basically unusable.** On uh3002 the spike and micro/macro LFP stores were saved as whole millivolts, so they only take a few flat values. That's how they were recorded, the code can't get the detail back.
- **Marker codes mean different things in different cohorts.** For BrainVision, `markers.py` reads each S/R code as a TTL line going up and down and decides per recording: a line up for a pulse width (60-240 µs) is a stim pulse, a line up 0.1-40 ms with many widths is a TTL-coded stim trial, and a line held 100 ms (the Buzzy vibration) or a fixed-width trigger is never called stim. Every marker in the cleaned file says its `type`, `meaning` and the `evidence`. Anything we couldn't check says `unknown`.
- **Generic BrainVision channel names stay unknown.** A lot of headers say Ch12 or just 17 (the amplifier input), or Leg, Mouth, SCS3. Nothing says what those are, so they go in the `unknown` stream until someone tells me the montage (add it to `COHORT_OVERRIDES` in `bv_channels.py`).
- **Stimulator files get lined up with the markers.** The MC_Stimulus programs, TTL tables, stim_struct.mat and stimuli_ordered.csv next to a BrainVision recording are matched to the stim markers by their timing. When they line up, the amplitude, pulse width and condition of every pulse go into `/stim`. Every file found is also kept as it is under `/sidecars/<name>`, matched or not.
- **Check NewSegment markers.** BrainVision writes a `NewSegment` marker when a recording is paused and restarted, but one recording I checked has 334 of them, some only 23 ms apart. They're saved under `/flags/new_segments` in the cleaned file. Figure out what they mean for your rig before you filter across them.
- **Some channels saturate.** Some TDT blocks hit the amplifier ceiling (about 0.537 V, or 32000 counts on whole-number stores) for a few seconds, and some sit there the whole block (uh3002 dbs1_clinicaltesting). Cleaning flags these under `/flags/saturation_<stream>` (channel, start, end) but doesn't change the data. A channel flat the whole block is only flagged if it sits at or past the stream's `rail`; flat somewhere else (like 0, unused) isn't. For BrainVision the rail is known for the actiCHamp (0.4096 V per electrode; a referenced channel can reach about twice that, -819 mV on ug3003 erp001) and for INT_16 files (32767 x the channel's resolution); any other amplifier writing floats has rail NaN, so only stretches inside the block are flagged there.

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
    bv_channels.py     what each BrainVision channel is, from its name (lead, ecog, eeg, emg, mic, MER, ...)
    markers.py         what each BrainVision marker code means per cohort, and which ones are stim
    sidecars.py        reads the stimulator and task files next to a BrainVision recording and lines them up with the markers
    stim_events.py     the /stim table (one row per pulse) for TDT, and the BrainVision clock fix
    block_info.py      patient, stage/session and output path for a recording, duplicate copies
    cohorts.py         one profile per cohort (layout, project, conditions), study ID patterns
    privacy.py         keeps names out of the cleaned set: known words, renames.csv, allowed_words.txt
    version.py         CLEAN_VERSION, bump it when the cleaned output changes
    index.py           table of the cleaned set, find_clean, quick QC
scripts/
    clean.py           cleans a folder (or one Cheaha array task)
    make_manifest.py   fixed list of recordings for Cheaha
    clean_slurm.sh     the Cheaha job
    copy_uhug.sh       copies only the uh/UG study folders from Box to Cheaha
    copy_cohorts.sh    copies the study folders of any other cohorts from Box to Cheaha's scratch
    check_copy.py      checks a copied raw folder for missing or cut-short signal files
matlab/
    clean_dbs_load.m   loads a cleaned file into a struct
tests/                 small tests, no recordings needed
docs/cheaha.md         running the whole inventory on Cheaha
examples/              examples (coming)
```

Recordings and anything made from them never go in this folder or in git. The `.gitignore` blocks the usual data files, but double check before you commit.

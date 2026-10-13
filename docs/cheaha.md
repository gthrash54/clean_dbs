###### WRITTEN FOR A NEW LAB MEMBER

# Cleaning the whole inventory on Cheaha

On my computer, `scripts/clean.py` cleans one block at a time. A full TDT block takes about 30 s and writes about 1.3 GB. Cleaning reads one window of the recording at a time, so it uses about 0.5 GB of RAM whatever the file size (an 8.4 GB, 69-channel brainvision file cleaned in 85 s at 0.56 GB, under a 4 GB cap), so on one computer the whole uh/UG inventory would take several hours. On Cheaha the first full run (574 uh/UG recordings in 58 tasks of 10) took about 15 minutes. Cheaha is UAB's computing cluster. It runs the cleaning as a SLURM **job array**: many copies of the same script run at once, and each copy (a "task") cleans its own few rows of a fixed list of blocks.

Don't put real patient folder names or paths in this repo, in issues, or in messages. The manifest and reports have patient IDs in them, so keep them in the out folder and never commit them.

Everything below runs in a terminal on Cheaha. Commands in `<angle brackets>` are placeholders for your own names.

## 1. Log in

You need a Cheaha account (request one at https://rc.uab.edu). Log in with your BlazerID, either:

- in the browser at https://rc.uab.edu (Open OnDemand), then **Clusters > Cheaha Shell**, or
- `ssh <BlazerID>@cheaha.rc.uab.edu` from a terminal.

You land on the **login node**. It's shared by everyone, so only do small things there (copying code, installing, submitting jobs). The actual cleaning always goes through `sbatch`.

If you fail to log in several times in a row, Cheaha can block your network for a while (ssh says "Connection refused", OnDemand says "Proxy Error"). Don't keep retrying, that makes it longer. Wait about an hour or use another network, like a phone hotspot. The login banner has a Disk Quota Report, which is quicker than running `du` on Cheaha's file system.

In a Linux terminal, copy is Ctrl+Shift+C and paste is Ctrl+Shift+V. Long pasted commands can break into two lines in a narrow window, so widen it first or put long paths in short variables (`M=$USER_SCRATCH/clean_out/manifest.csv`). A lone `sbatch` with no script after it just waits for input forever; press Ctrl+C and check `squeue -u $USER` before submitting again.

## 2. Get the code and its Python

Cheaha has a few storage areas (from the Research Computing storage docs, checked Oct 2026):

| where | path | good for |
|---|---|---|
| user data | `$USER_DATA` (`/data/user/<BlazerID>`) | the code, the venv, raw data if there's room |
| project | `/data/project/<lab project>` | shared lab data, has to be requested by the PI |
| scratch | `$USER_SCRATCH` (`/scratch/<BlazerID>`) | cleaned output while you work on it. **Files not touched for 90 days get deleted**, so copy results somewhere permanent |

Put the code in `$USER_DATA`:

```
cd $USER_DATA
git clone https://github.com/gthrash54/clean_dbs.git clean_dbs
cd clean_dbs
```

The repo is private, so `git clone` and later `git pull` ask for a password. GitHub won't take your normal one. Make a fine-grained token on GitHub (Settings > Developer settings) with read-only access to this repo's contents and paste that as the password. If git complains it "cannot open display", run `unset SSH_ASKPASS` first.

Then make the venv the same way as in the README. Cheaha doesn't (as far as I know) have a uv module, but uv installs into your own folder without admin rights:

```
curl -LsSf https://astral.sh/uv/install.sh | sh
source ~/.bashrc
uv venv --python 3.12
uv pip install --only-binary :all: -e .
```

`--only-binary :all:` makes it use prebuilt packages. Cheaha's system is old, so building numpy or pandas from source there fails. With it you get numpy 2.2.6 and pandas 2.3.2, and the tests pass on those. Everything else on Cheaha is old too (bash 4.2, git 1.8.3). The scripts are written to work with those, so you don't need to load a newer git or bash. If something fails there that works on your own computer, an old tool is the first thing I'd suspect.

If uv doesn't work there, use Cheaha's Anaconda module to get a Python and make a plain venv instead:

```
module load Anaconda3
python -m venv .venv
source .venv/bin/activate
pip install --only-binary :all: -e .
```

Either way you end up with `clean_dbs/.venv`, which is what the job script switches on. Check it works:

```
source .venv/bin/activate
python -c "import clean_dbs; print('ok')"
```

## 3. Copy the recordings from Box

Research Computing's way to get data from Box is rclone. Set up a Box "remote" once:

```
module load rclone
rclone config
```

Make a new remote called `box`, pick **Box** as the type, leave the rest blank, and answer **n** to "Use auto config?". It then waits for a token. On **your own computer** (it needs a browser), install rclone and run `rclone authorize "box"`, sign in with Single Sign On using your UAB email and BlazerID, and paste the token it prints (between `--->` and `<---`) back into Cheaha. Check it with `rclone lsd box:`.

Then copy **only the uh/UG study folders** (uh30xx and UG30xx, in any case). The Box inventory also holds other studies and some legacy folders named after patients; those must not be copied, so never copy the whole inventory folder. `Imaging` folders are left out too, because DICOM headers can carry patient names. rclone's include filters didn't pick folders reliably, so `scripts/copy_uhug.sh` lists the study folders first and copies them one by one:

```
rclone lsf --dirs-only --max-depth 2 "box:UAB DBS Data/UAB DBS Raw Data" | grep -i -E '^[0-9.]+/(uh|ug)30[0-9]{2}/$' > $USER_DATA/dbs_raw/uhug_folders.txt
while read d; do rclone copy --progress --ignore-case --exclude "imaging/**" --exclude "**/imaging/**" "box:UAB DBS Data/UAB DBS Raw Data/$d" "$USER_DATA/dbs_raw/$d"; done < $USER_DATA/dbs_raw/uhug_folders.txt
```

The 31 uh/UG folders come to about 1 TB without imaging, so the cleaned set can be bigger than you'd guess too; scratch has plenty of room for it. The list should have 31 lines (the script prints the count; don't print the list itself, it has the Box number folders in it). It's kept in `$USER_DATA/dbs_raw`, never in the repo. The copy is big and slow, so run it as a job; that way it keeps going after you log out:

```
cd $USER_DATA/clean_dbs
mkdir -p logs
sbatch --partition=short --time=12:00:00 --mem=2G --job-name=uhug_copy --output=logs/uhug_copy_%j.out --wrap "bash scripts/copy_uhug.sh"
```

Rerunning it is safe: `rclone copy` skips files that are already there, so a stopped copy just picks up where it left off. Use `copy`, never `sync`. `sync` deletes files on the destination to make it match, and you can't undo that. If rclone later says the token expired, run `rclone config reconnect box:`.

### Check the copy

A copy that stopped partway can leave a `.tev` cut short. tdt reads it without an error and fills the missing part with zeros, so check the whole copy before cleaning:

```
source .venv/bin/activate
python scripts/check_copy.py $USER_DATA/dbs_raw
```

It lists every recording missing its signal file (`.tev` / `.eeg`) or with a `.tev` shorter than its `.tsq` index says, then a count. Rerun the copy job for anything it lists, then check again. Cleaning refuses a cut-short `.tev` too (it fails with the same message), so nothing zero-filled gets into the cleaned set.

## 4. Make the manifest

The manifest is a csv with one row per recording (path, patient, stage, block, format, which signal files are missing). Task 0 cleans the first rows, task 1 the next ones, and so on, so the list can't change while the job runs.

```
cd $USER_DATA/clean_dbs
source .venv/bin/activate
mkdir -p $USER_SCRATCH/clean_out
python scripts/make_manifest.py $USER_DATA/dbs_raw $USER_SCRATCH/clean_out/manifest.csv --chunk 10
```

It prints how many recordings it found, how many are missing signal files, how much signal there is, the biggest peak memory estimate (`est_peak_gb`, also a column, next to `signal_gb`), and the `--array` range to use. It stops with an error if two recordings would write the same cleaned file (names that differ only in capitals count as the same, since Mac and Windows treat them that way); fix those folder names first. `--chunk 10` means 10 blocks per task, so about 550 blocks becomes 55 tasks of roughly 5 minutes each. Recordings inside an `Imaging` folder are skipped, so scans left over from an older copy don't matter. Add `--skip-missing` to leave out blocks that are missing their `.tev` / `.eeg` (it prints how many it skipped). Without it, those blocks just fail quickly and show up in the report as failed, so you can see what still needs copying.

### Put the case tables in the out folder

Every cleaned file gets the targets and a `case_status` from `_cases/cases.csv`, so `clean.py` won't start without it. Copy `_cases` (with `cases.csv` and `blocks.csv`) into the out folder before submitting:

```
rclone copy "box:<cleaned folder on Box>/_cases" $USER_SCRATCH/clean_out/_cases
ls $USER_SCRATCH/clean_out/_cases     # cases.csv  blocks.csv
```

For a cohort with no case tables, add `--no-cases` after the chunk size in the `sbatch` line (the targets then stay blank and `case_status` is `missing`).

## 5. Submit

The code has to be committed: every cleaned file records the git commit that made it, so `clean.py --manifest` refuses to run if the commit is unknown or the code has uncommitted changes. `git status` in the `clean_dbs` folder should say nothing to commit (`--allow-dirty` overrides it, for test runs only).

From the `clean_dbs` folder (the job looks for `.venv` in the folder you submit from):

```
mkdir -p logs
sbatch --array=0-0 scripts/clean_slurm.sh $USER_SCRATCH/clean_out/manifest.csv $USER_SCRATCH/clean_out 10
```

Try just task 0 first like that. When it works, submit the rest with the range `make_manifest.py` printed (it won't redo task 0's blocks, see step 8):

```
sbatch --array=0-54 scripts/clean_slurm.sh $USER_SCRATCH/clean_out/manifest.csv $USER_SCRATCH/clean_out 10
```

The three things after the script are the manifest, the out folder, and the chunk size. Anything after the chunk size goes to `clean.py` (for example `--force`, see step 8). The chunk **has to be the same number** you gave `make_manifest.py`, or the tasks won't line up with the `--array` range. To cap how many tasks run at the same time, add `%`, for example `--array=0-54%20`.

The defaults are set at the top of `scripts/clean_slurm.sh`: `express` partition (2 hour max), 1 hour per task, 4 GB of memory, 1 CPU. 4 GB is enough for every recording: cleaning reads one window at a time, so a block peaks at about 0.5 GB whatever its size (`est_peak_gb` in the manifest). In the first full run each task used 230 to 505 MB and took 1 to 9 minutes, so the defaults have plenty of room. You can change any of them when you submit, for example `sbatch --time=02:00:00 ...`.

## 6. Watch it

```
squeue -u $USER                                             # what's waiting (PD) and running (R)
sacct -j <job id> --format=JobID,State,Elapsed,MaxRSS       # how each task went, after it finishes
tail logs/clean_dbs_<job id>_<task>.out                     # what one task printed
```

`sbatch` prints the job id when you submit. In `sacct`, `COMPLETED` means every block in that task cleaned (or was already done), `FAILED` means at least one block in it failed (the reason is in the report and the `.out` log), `OUT_OF_MEMORY` means it needs more `--mem`, and `TIMEOUT` means it needs more `--time` or a smaller chunk. The other recordings in a FAILED task still clean.

sacct also lists extra step lines for each task (the ones with a `.` in the JobID). To see one line per task:

```
sacct -j <job id> --format=JobID,State,Elapsed,MaxRSS | grep -v '\.'
```

`NODE_FAIL` means the computer running that task died. You don't have to do anything: Slurm puts the task back in the queue and it runs again on another node. Blocks it already finished are skipped and a half-written one is redone, because a `.h5` only appears once it's complete. To see both tries of a task:

```
sacct -j <job id>_<task> --duplicates --format=JobID,State,Elapsed,NodeList
```

Each task writes its own report to `<out>/_reports/clean_report_<task>_chunk<chunk>_<manifest hash>.csv` after every block, so tasks never write over each other, and a task that gets killed still reports the blocks it finished. The chunk size and a short hash of the manifest are in the name, so reports from another manifest never get mixed in. A task number past the end of the manifest prints a note and exits without error.

## 7. Merge the reports

When everything is done:

```
python scripts/clean.py --merge-reports $USER_SCRATCH/clean_out --manifest $USER_SCRATCH/clean_out/manifest.csv
```

`--manifest` is required. That combines the task reports made from this manifest into `<out>/clean_report.csv`, one row per manifest row, in manifest order, and prints how many of each status. When a row was run more than once (a rerun, or another chunk size), only the newest try is kept; rows that aren't in the manifest are left out with a warning. A block with no report at all (usually a task that was killed before it finished its first block) gets a "not run" row.

| status | means |
|---|---|
| `cleaned` | cleaned for the first time |
| `skipped` | already cleaned by this version of the code, left as is |
| `recleaned` | there already, but cleaned again; `reason` says why (wouldn't open, older `clean_version`, or `--force`) |
| `duplicate` | another copy of the same recording (same size) under a different number folder, cleaned once from the copy named in `reason` |
| `conflict` | same name and place as another copy but a different size; the bigger one is cleaned, check these by hand |
| `failed` | `reason` says why (paths in it start at the study ID) |
| `excluded` | never cleaned: no study ID above it or in its name, or a folder or file name that isn't cleared (see "Cleaning the other cohorts") |
| `not run` | no report for that row |

## 8. Rerun what failed

Rerunning is safe. `clean_block` skips a block whose `.h5` is already there **and** was made by the current version of the cleaning (`clean_version`, see the README); a `.h5` that won't open or comes from an older version is cleaned again and reported as `recleaned`. It only writes the `.h5` once the whole file is done (it writes `<block>.h5.part` and renames it at the end), so a killed task never leaves a half-finished `.h5` that looks done.

- To rerun everything that didn't finish: submit the same `sbatch` command again. Finished blocks are skipped in seconds and show up as "skipped" in the new report.
- To rerun only some tasks: give just those task numbers, for example `sbatch --array=3,17,40 ...`. Find them with `sacct` or in the "row" column of the report (task = row // chunk).
- To clean everything again on purpose (say a fix that didn't bump `clean_version`), add `--force` after the chunk size: `sbatch --array=0-54 scripts/clean_slurm.sh <manifest> <out> 10 --force`.
- If tasks ran out of time, add `--time=02:00:00` to the `sbatch` line.

Then merge again (step 7). The merged report keeps the newest try of every row, so it's fine to rerun with a different chunk size.

## 9. Keep the results

Scratch deletes files that haven't been touched in 90 days. When the merged report looks right, move the cleaned set somewhere permanent: your `$USER_DATA` (5 TB) or, better, the lab's project folder (`/data/project/<lab>`, which the PI requests from Research Computing), so everyone in the lab can read it. Keep a second copy somewhere else too, for example on Box:

```
rclone copy --progress $USER_SCRATCH/clean_out "box:<cleaned folder on Box>"
```

RC's storage docs don't say whether user or project folders are backed up, so don't count on it until RC says so.

You don't have to move the data to look at it. Analysis can run on Cheaha too (JupyterLab or MATLAB through OnDemand, or `sbatch`).

## Cleaning the other cohorts

Everything above was written for uh/UG. Box has 24 more cohorts (DLPFC, BI, PP, MC, BV, ORBV, BPP, BB, CHR, PT, SEP, PNS, TM, TMSCS, LFP, JOY, ECOG, burst, DYS, sensory, PR, MOTIF, SDC and the BrainVision part of OR; `clean_dbs/cohorts.py` has what each one seems to be). They use the same scripts, with three differences: the raw copy goes to scratch, one cohort at a time; there are no case tables for them, so they run with `--no-cases`; and their folder and file names have to pass a names check first, because some of them have a person's first name in them.

All the cleaned files go into the same out folder as uh/UG (`$USER_SCRATCH/clean_out`), one manifest per cohort next to each other (`manifest.csv` stays uh/UG, then `manifest_dlpfc.csv`, `manifest_bi.csv`, ...). Every patient folder is a study ID, so two cohorts can't write the same file unless a study ID sits inside another cohort's folder (LFP and SEP IDs inside PP or PT folders) and also has a folder of its own. `make_manifest.py` checks the other manifests in the out folder for that: the same recording copied twice is marked `duplicate` and cleaned once, and a different recording that would land on the same file stops it.

### How big it is

The sizes are what `copy_cohorts.sh` would copy, and the recordings and tasks are from a dry run of `make_manifest.py`, both worked out from a listing of Box in October 2026.

| cohort | study IDs | copy (GB) | recordings | tasks at `--chunk 10` |
|---|---:|---:|---:|---:|
| bi | 32 | 950 | 851 | 86 |
| dlpfc | 154 | 325 | 528 | 53 |
| pp (with the IDs inside it) | 42 | 357 | 187 | 19 |
| or (BrainVision part only) | 65 | 215 | 29 | 3 |
| mc | 14 | 151 | 98 | 10 |
| bv | 4 | 119 | 51 | 6 |
| bpp | 20 | 59 | 22 | 3 |
| bb | 31 | 45 | 73 | 8 |
| orbv | 2 | 39 | 29 | 3 |
| pt, chr, ecog, lfp, sep, pns, tmscs, joy, burst, tm, sensory, sdc, motif, pr, dys | 2 to 12 each | 1 to 35 each | 1 to 32 each | 1 to 4 each |
| **all 24** | **about 420** | **about 2,500** | **about 2,090** | **about 220** |

That's about 3.5 times the uh/UG run (574 recordings, 637 GB of signal, 58 tasks of 10). The recordings are about the same size on average, so I'd expect each task to take 1 to 15 minutes and use 230 to 505 MB like the uh/UG ones did, and a whole cohort to finish in well under an hour once its tasks start. On a sample of 43 BrainVision recordings from these cohorts the cleaned files came to half the size of the signal, so the cleaned set should be 1 to 2 TB. Scratch has room for all of it, but the raw copies add up, which is why I do one cohort at a time and delete each raw copy when its cohort is done.

### Order of work, one cohort at a time

I'd start with a small cohort (PNS, 20 recordings) to go around the whole loop once, then go from small to big and do BI last. For each cohort (here DLPFC):

1. **Copy it to scratch.** `scripts/copy_cohorts.sh` works like `copy_uhug.sh`, for any cohort. Try `--list-only` first; it prints how many folders and files it would copy and copies nothing:

    ```
    cd $USER_DATA/clean_dbs
    bash scripts/copy_cohorts.sh --list-only dlpfc
    mkdir -p $USER_DATA/copy_logs
    sbatch --partition=short --time=12:00:00 --mem=2G --job-name=copy_dlpfc --output=$USER_DATA/copy_logs/%x_%j.out --wrap "bash scripts/copy_cohorts.sh dlpfc"
    ```

    It copies `<number folder>/<study ID>/` for that cohort's IDs (in any case, `DLPFC045`, `dlpfc045`, `SEP_03`), with any IDs nested inside them, into `$USER_SCRATCH/dbs_raw/dlpfc/`. It only looks inside number folders, never in the folders at the top of Box that are named after patients. Number folders with a short letter suffix after the number count too (49 study IDs are only there); `--plain-only` leaves those out. It also takes recordings sitting straight in a number folder with the cohort's ID in their name (`BPP012.vhdr`), folders named after two study IDs (with the cohort of the first one), and, for MOTIF, PR and burst, the folders with no study ID that belong to them. It leaves out Imaging, Images, DICOM and MRI folders, scan files, documents, photos and videos other than TDT's `.avi`. The lists it makes have the number folders in them, so they stay next to the copy and are never printed, and the log only shows totals. On my laptop the remote is called `Box`, so there it's `REMOTE=Box DEST=<some folder> bash scripts/copy_cohorts.sh --list-only dlpfc`. If a copy job runs out of time, submit it again; it picks up where it stopped. Copy only, never sync.

2. **Check the copy.** `python scripts/check_copy.py $USER_SCRATCH/dbs_raw/dlpfc`. Names it can't clear (below) are printed as `<withheld>`.

3. **Make the manifest.** Same as step 4, with the cohort's raw folder and its own manifest name:

    ```
    python scripts/make_manifest.py $USER_SCRATCH/dbs_raw/dlpfc $USER_SCRATCH/clean_out/manifest_dlpfc.csv --chunk 10
    ```

    It can take more than one raw folder (`make_manifest.py <raw 1> <raw 2> <manifest>`), and so can `check_copy.py`. If it stops at the names check, do the privacy step below and run it again. When it works it prints the whole `sbatch` line to use, with `--no-cases` when none of the cohort's IDs are in `_cases/cases.csv`.

4. **Test task 0**, and check its `Elapsed` and `MaxRSS`:

    ```
    sbatch --array=0-0 scripts/clean_slurm.sh $USER_SCRATCH/clean_out/manifest_dlpfc.csv $USER_SCRATCH/clean_out 10 --no-cases
    ```

5. **Submit the rest** with the range it printed, e.g. `--array=0-52`.

6. **Merge.** `python scripts/clean.py --merge-reports $USER_SCRATCH/clean_out --manifest $USER_SCRATCH/clean_out/manifest_dlpfc.csv` writes `clean_report_dlpfc.csv` (the report is named after the manifest, so cohorts don't write over each other's).

7. **Move the cleaned set off scratch** (step 9). Copying the whole out folder again only adds what's new, for example `rsync -a $USER_SCRATCH/clean_out/ /data/project/<lab>/clean_out/`.

8. **Delete the raw copy** once the merged report looks right: `rm -rf $USER_SCRATCH/dbs_raw/dlpfc`. Box is the master, so it can always be copied again.

### Recordings that can't be cleaned yet

A recording has to be filed under a study ID (a folder above it, or the ID in its file name) for clean_dbs to know whose it is, and it never guesses. Recordings in the MOTIF, Motif protocol, PR, Reversal_Polarity and BurstPairPulses folders, and in folders named after two study IDs, are in the manifest with `excluded` and the reason, and show up in the report as `excluded`; nothing in them is cleaned until they're filed under one study ID on Box. A few things aren't copied at all: recordings sitting loose in a number folder with no ID in their name, a LOC folder, a number folder inside a number folder, and the DBS and AO folders (their files aren't BrainVision or TDT). The OR folders get copied, but only one of them has BrainVision files; clean_dbs doesn't see the `T###.eeg` files with no `.vhdr` in the rest, so you could leave OR until someone knows what recorded those.

### The privacy step: names in folder and file names

Some folder and file names on Box have a person's first name in them (a recording saved under the patient's name, a subfolder named after someone). Whatever is in a name goes into the cleaned folder and file names, the `source` and `block` in every file's info, the stimulator files stored in `/sidecars`, and the reports. So `make_manifest.py` checks every name that would go into the cleaned set: the folders from the study ID down, the recording's name, the signal and marker files its header names, the stimulator and task files that get stored with it, and TDT videos. Every word in a name has to be one that describes recordings (the list is `WORDS` in `clean_dbs/privacy.py`); numbers, single letters and study IDs are fine. If a name has a word it doesn't know, it stops before writing anything, prints how many names per study ID (never the names), and lists them in `_cases/names_to_check.csv` in the out folder. It errs on the side of stopping: an abbreviation it hasn't seen stops it just like a name does.

To fill it in:

1. Open `$USER_SCRATCH/clean_out/_cases/names_to_check.csv` on Cheaha (not on a shared screen). Each row has the study ID, what kind of name it is (folder, recording, data file, sidecar, video), the name as it is on Box (`original`), and a path to find it (`example`).
2. If it has a person's name in it, write a new name in `renamed` that says what the recording is, from words in the list, and keep its numbers so recordings stay apart: `rest_2`, `task_axcpt_3`. If you can't tell what it is, use the kind and a number: `record1`, `folder1`, `stimfile1`.
3. If it isn't a name, just a word the list doesn't have (an abbreviation, a typo, a task name), leave the row out and add the word to `_cases/allowed_words.txt` instead, one word per line. Never add a word that's mostly used as a first name.
4. Copy the rows you filled in into `_cases/renames.csv` (columns `study_id`, `original`, `renamed`; a blank `study_id` applies to every study). A renamed value has to pass the check itself, and can't have a `/` in it.
5. Run `make_manifest.py` again. It prints how many renames and allowed words it used.

Renames only change the cleaned names, never Box or the raw copy. `clean.py` reads the same files from the out folder you give it, checks every recording again before cleaning it, and refuses a recording whose output path isn't the one in the manifest, so if `renames.csv` changes after you made the manifest, make it again. To clean the rest of a cohort while someone sorts out the names, add `--exclude-unchecked` to `make_manifest.py`: those recordings are then in the manifest and the report as `excluded`, with their name withheld, and nothing of them is cleaned.

`renames.csv` and `names_to_check.csv` have the original names in them, and the manifests have the raw paths with the Box number folders, so leave those out whenever the cleaned set goes somewhere other people can open (the case tables in `_cases` stay with it, `build_index` needs them):

```
rclone copy --progress --exclude "_cases/renames.csv" --exclude "_cases/names_to_check.csv" --exclude "manifest*.csv" $USER_SCRATCH/clean_out "box:<cleaned folder on Box>"
```

The check covers uh/UG as well, and some uh/UG names don't pass it yet, so the next time you make the uh/UG manifest it will stop until those are cleared the same way.

## Looking at the data in JupyterLab

Once, from the `clean_dbs` folder:

```
source .venv/bin/activate
uv pip install ipykernel --only-binary :all:
python -m ipykernel install --user --name clean_dbs --display-name "clean_dbs"
```

That only adds packages, it doesn't change the ones the cleaning uses, and `.venv` isn't tracked by git, so `git status` stays clean. Then at https://rc.uab.edu go to **Interactive Apps > JupyterLab** and pick the **clean_dbs** kernel. I'd use JupyterLab over Notebook, since it has a file browser, tabs, a terminal and a CSV viewer.

## Running it on your own computer instead

The array mode works anywhere, it's just slower. Run the tasks one after another by giving `--task` yourself:

```
python scripts/clean.py --manifest manifest.csv --task 0 --chunk 10 <out folder>
```

The original `python scripts/clean.py <raw folder> <out folder>` still cleans everything in one go and writes `clean_report.csv` like before.

## Where this came from

UAB Research Computing docs, checked October 2026. Things marked "not checked" I haven't confirmed on Cheaha myself.

- Partitions and time limits (express 2 h, short 12 h, medium 50 h, long 150 h): https://docs.rc.uab.edu/cheaha/hardware/
- Job arrays, `%A` / `%a` log names, `%` throttling: https://docs.rc.uab.edu/cheaha/slurm/submitting_jobs/
- Storage paths, scratch 90 day purge, Cheaha accepting PHI (Long-Term Storage does not): https://docs.rc.uab.edu/data_management/storage/
- rclone with Box: https://docs.rc.uab.edu/data_management/transfer/rclone/
- `module load Anaconda3`: https://docs.rc.uab.edu/workflow_solutions/using_anaconda/
- Checked on Cheaha in October 2026: uv installs work (with `--only-binary :all:`), the ssh hostname, and the run times and memory above. Still run task 0 alone first and check its `Elapsed` and `MaxRSS` in `sacct` before submitting everything.

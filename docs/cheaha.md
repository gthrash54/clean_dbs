###### WRITTEN FOR A NEW LAB MEMBER

# Cleaning the whole inventory on Cheaha

On my computer, `scripts/clean.py` cleans one block at a time. A full TDT block takes about 30 s, uses up to 1.8 GB of RAM and writes about 1.3 GB, so the whole inventory (about 550 recordings) is several hours and around a terabyte of disk. Cheaha is UAB's computing cluster. It runs the cleaning as a SLURM **job array**: many copies of the same script run at once, and each copy (a "task") cleans its own few rows of a fixed list of blocks.

Don't put real patient folder names or paths in this repo, in issues, or in messages. The manifest and reports have patient IDs in them, so keep them in the out folder and never commit them.

Everything below runs in a terminal on Cheaha. Commands in `<angle brackets>` are placeholders for your own names.

## 1. Log in

You need a Cheaha account (request one at https://rc.uab.edu). Log in with your BlazerID, either:

- in the browser at https://rc.uab.edu (Open OnDemand), then **Clusters > Cheaha Shell**, or
- `ssh <BlazerID>@cheaha.rc.uab.edu` from a terminal.

You land on the **login node**. It's shared by everyone, so only do small things there (copying code, installing, submitting jobs). The actual cleaning always goes through `sbatch`.

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

Then make the venv the same way as in the README. Cheaha doesn't (as far as I know) have a uv module, but uv installs into your own folder without admin rights:

```
curl -LsSf https://astral.sh/uv/install.sh | sh
source ~/.bashrc
uv venv --python 3.12
uv pip install -e .
```

If uv doesn't work there, use Cheaha's Anaconda module to get a Python and make a plain venv instead:

```
module load Anaconda3
python -m venv .venv
source .venv/bin/activate
pip install -e .
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

Then copy the inventory. Do a dry run first, it shows what would copy without copying:

```
mkdir -p $USER_DATA/dbs_raw
rclone copy --dry-run "box:<inventory folder on Box>" $USER_DATA/dbs_raw
rclone copy --progress "box:<inventory folder on Box>" $USER_DATA/dbs_raw
```

Use `copy`, never `sync`. `sync` deletes files on the destination to make it match, and you can't undo that. The copy is close to a terabyte, so it takes a while. If the login node stops it, ask Research Computing what they want us to use for big transfers (they also support Globus). If rclone later says the token expired, run `rclone config reconnect box:`.

## 4. Make the manifest

The manifest is a csv with one row per recording (path, patient, stage, block, format, which signal files are missing). Task 0 cleans the first rows, task 1 the next ones, and so on, so the list can't change while the job runs.

```
cd $USER_DATA/clean_dbs
source .venv/bin/activate
mkdir -p $USER_SCRATCH/clean_out
python scripts/make_manifest.py $USER_DATA/dbs_raw $USER_SCRATCH/clean_out/manifest.csv --chunk 10
```

It prints how many recordings it found, how many are missing signal files, and the `--array` range to use. `--chunk 10` means 10 blocks per task, so about 550 blocks becomes 55 tasks of roughly 5 minutes each. Add `--skip-missing` to leave out blocks that are missing their `.tev` / `.eeg` (it prints how many it skipped). Without it, those blocks just fail quickly and show up in the report as failed, so you can see what still needs copying.

## 5. Submit

From the `clean_dbs` folder (the job looks for `.venv` in the folder you submit from):

```
mkdir -p logs
sbatch --array=0-0 scripts/clean_slurm.sh $USER_SCRATCH/clean_out/manifest.csv $USER_SCRATCH/clean_out 10
```

Try just task 0 first like that. When it works, submit the rest with the range `make_manifest.py` printed (it won't redo task 0's blocks, see step 8):

```
sbatch --array=0-54 scripts/clean_slurm.sh $USER_SCRATCH/clean_out/manifest.csv $USER_SCRATCH/clean_out 10
```

The three things after the script are the manifest, the out folder, and the chunk size. The chunk **has to be the same number** you gave `make_manifest.py`, or the tasks won't line up with the `--array` range. To cap how many tasks run at the same time, add `%`, for example `--array=0-54%20`.

The defaults are set at the top of `scripts/clean_slurm.sh`: `express` partition (2 hour max), 1 hour per task, 4 GB of memory, 1 CPU. You can change any of them when you submit, for example `sbatch --mem=8G --time=02:00:00 ...`.

## 6. Watch it

```
squeue -u $USER                                             # what's waiting (PD) and running (R)
sacct -j <job id> --format=JobID,State,Elapsed,MaxRSS       # how each task went, after it finishes
tail logs/clean_dbs_<job id>_<task>.out                     # what one task printed
```

`sbatch` prints the job id when you submit. In `sacct`, `COMPLETED` means every block in that task cleaned (or was already done), `FAILED` means at least one block in it failed (the reason is in the report and the `.out` log), `OUT_OF_MEMORY` means it needs more `--mem`, and `TIMEOUT` means it needs more `--time` or a smaller chunk.

Each task writes its own report to `<out>/_reports/clean_report_<task>.csv` after every block, so tasks never write over each other, and a task that gets killed still reports the blocks it finished.

## 7. Merge the reports

When everything is done:

```
python scripts/clean.py --merge-reports $USER_SCRATCH/clean_out --manifest $USER_SCRATCH/clean_out/manifest.csv
```

That combines all the task reports into `<out>/clean_report.csv`, one row per manifest row, in manifest order, and prints how many cleaned, skipped, and failed. With `--manifest` it also adds a "not run" row for any block that has no report (usually a task that was killed before it finished its first block).

## 8. Rerun what failed

Rerunning is safe. `clean_block` skips any block whose `.h5` is already there, and it only writes the `.h5` once the whole file is done (it writes `<block>.h5.part` and renames it at the end). So a killed task never leaves a half-finished `.h5` that looks done, and finished blocks are never redone.

- To rerun everything that didn't finish: submit the same `sbatch` command again. Finished blocks are skipped in seconds and show up as "skipped" in the new report.
- To rerun only some tasks: give just those task numbers, for example `sbatch --array=3,17,40 ...`. Find them with `sacct` or in the "row" column of the report (task = row // chunk).
- If tasks ran out of memory or time, add `--mem=8G` or `--time=02:00:00` to the `sbatch` line.

Then merge again (step 7). A rerun task replaces its own report file, so the merged report always shows the latest try. If you change the chunk size, delete `<out>/_reports` first so old reports don't mix with new ones.

## 9. Keep the results

Scratch deletes files after 90 days, so when the cleaned set looks right, copy it somewhere permanent, for example back to Box:

```
rclone copy --progress $USER_SCRATCH/clean_out "box:<cleaned folder on Box>"
```

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
- Not checked: that the uv installer works on Cheaha, the exact ssh hostname, whether big rclone copies are allowed on the login node, and the time per block on Cheaha's file system (30 s is on my computer). Run task 0 alone first and check its `Elapsed` and `MaxRSS` in `sacct` before submitting everything.

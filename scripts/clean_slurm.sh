#!/bin/bash
# runs scripts/clean.py as a slurm job array on cheaha, each array task cleans CHUNK rows of the manifest
#
# submit from the clean_dbs folder (the one with pyproject.toml), like:
#   mkdir -p logs
#   sbatch --array=0-57 scripts/clean_slurm.sh <manifest.csv> <out folder> 10
# anything after the chunk size goes to clean.py, e.g. --force to clean everything again
# make_manifest.py prints the right --array range for your chunk size (0-57 was the uh/UG run: 574 recordings, 10 per task)
# the #SBATCH lines below are defaults, anything given to sbatch on the command line wins

#SBATCH --job-name=clean_dbs
# express is cheaha's shortest queue (2 hour max) so jobs start fast; 10 blocks at ~30 s each fits easily
#SBATCH --partition=express
# max time per task (hours:minutes:seconds), slurm kills the task after this; raise it if you use a bigger chunk
#SBATCH --time=01:00:00
# cleaning reads one window at a time, so a block peaks at ~0.5 GB of RAM whatever its size (13.5 GB brainvision included); 4G leaves room, and blocks run one at a time so this doesnt grow with chunk
#SBATCH --mem=4G
# one task = one program on one cpu, the cleaning code doesnt use more than one cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# where printed output goes, %x = job name, %A = array job id, %a = task number; the logs folder has to exist before you submit
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# stop at the first command that fails, and treat unset variables as errors
set -euo pipefail

# the three things given after the script name in the sbatch command
MANIFEST="$1"
OUT="$2"
CHUNK="${3:-1}" # rows per task, 1 if left out
EXTRA=("${@:4}") # anything after the chunk goes to clean.py as is, e.g. --force or --no-cases

# slurm starts the job in the folder sbatch was run from, which should be the clean_dbs folder
cd "$SLURM_SUBMIT_DIR"

# switch on the toolbox's python (the .venv made with uv during setup)
source .venv/bin/activate

# read_tdt already gives Chicago time whatever the machine's zone, this is just so dates in the log match the rig's clock
export TZ=America/Chicago

# the case tables go into every file (targets, case_status); stop here instead of in every task's python (clean.py checks too)
# ${EXTRA[@]+...} looks odd but is needed: cheaha's bash (4.2) calls an empty array "unbound" under set -u, newer bash doesn't
if [ ! -f "$OUT/_cases/cases.csv" ] && [[ " ${EXTRA[*]+${EXTRA[*]}} " != *" --no-cases "* ]]; then
    echo "no $OUT/_cases/cases.csv, copy _cases (cases.csv, blocks.csv, and contacts.csv for the uh contact coordinates) into the out folder first" >&2
    exit 1
fi
# the lab's do-not-use list goes with them (channel_dnu); --no-cases or --no-dnu clean without it on purpose
if [ -f "$OUT/_cases/cases.csv" ] && [ ! -f "$OUT/_cases/dnu_channels.csv" ] && [[ " ${EXTRA[*]+${EXTRA[*]}} " != *" --no-cases "* ]] \
        && [[ " ${EXTRA[*]+${EXTRA[*]}} " != *" --no-dnu "* ]]; then
    echo "no $OUT/_cases/dnu_channels.csv, copy it in with the other case tables (or add --no-dnu to clean without the lab's do-not-use list)" >&2
    exit 1
fi

# show which node and task this is at the top of the log, handy when something goes wrong
echo "task $SLURM_ARRAY_TASK_ID of job $SLURM_ARRAY_JOB_ID on $(hostname), started $(date)"

# clean.py reads SLURM_ARRAY_TASK_ID itself to pick its rows of the manifest
python scripts/clean.py --manifest "$MANIFEST" --chunk "$CHUNK" ${EXTRA[@]+"${EXTRA[@]}"} "$OUT" #same empty-array guard as above

echo "finished $(date)"

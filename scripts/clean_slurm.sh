#!/bin/bash
# runs scripts/clean.py as a slurm job array on cheaha, each array task cleans CHUNK rows of the manifest
#
# submit from the clean_dbs folder (the one with pyproject.toml), like:
#   mkdir -p logs
#   sbatch --array=0-54 scripts/clean_slurm.sh <manifest.csv> <out folder> 10
# make_manifest.py prints the right --array range for your chunk size
# the #SBATCH lines below are defaults, anything given to sbatch on the command line wins

#SBATCH --job-name=clean_dbs
# express is cheaha's shortest queue (2 hour max) so jobs start fast; 10 blocks at ~30 s each fits easily
#SBATCH --partition=express
# max time per task (hours:minutes:seconds), slurm kills the task after this; raise it if you use a bigger chunk
#SBATCH --time=01:00:00
# one block peaks at ~1.8 GB of RAM, 4G leaves room; blocks run one at a time so this doesnt grow with chunk
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

# slurm starts the job in the folder sbatch was run from, which should be the clean_dbs folder
cd "$SLURM_SUBMIT_DIR"

# switch on the toolbox's python (the .venv made with uv during setup)
source .venv/bin/activate

# show which node and task this is at the top of the log, handy when something goes wrong
echo "task $SLURM_ARRAY_TASK_ID of job $SLURM_ARRAY_JOB_ID on $(hostname), started $(date)"

# clean.py reads SLURM_ARRAY_TASK_ID itself to pick its rows of the manifest
python scripts/clean.py --manifest "$MANIFEST" --chunk "$CHUNK" "$OUT"

echo "finished $(date)"

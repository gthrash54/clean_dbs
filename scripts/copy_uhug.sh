#!/bin/bash
# copies only the uh30xx / UG30xx study folders from Box to $USER_DATA/dbs_raw, skipping any Imaging folders
# other studies and the patient-named legacy folders on Box must not be copied, and DICOM headers in Imaging can carry patient names
#
# run it from the clean_dbs folder as a slurm job, so it keeps going after you log out:
#   mkdir -p logs
#   sbatch --partition=short --time=12:00:00 --mem=2G --job-name=uhug_copy --output=logs/uhug_copy_%j.out --wrap "bash scripts/copy_uhug.sh"
# rerunning is safe, rclone copy skips files that are already there; then check the copy with scripts/check_copy.py

module load rclone 2>/dev/null || true # on cheaha rclone comes from a module (before set -u, module scripts trip over unset variables)
set -euo pipefail

SRC="box:UAB DBS Data/UAB DBS Raw Data"
DEST="$USER_DATA/dbs_raw"
mkdir -p "$DEST"

# the study folders sit one level down (<number folder>/<study ID>/); rclone's include filters didn't pick folders reliably, so list them
# and keep only uh30xx / ug30xx in any case (should be 31). The list has number folder names, so it stays out of the repo and isn't printed
LIST="$DEST/uhug_folders.txt" # in dbs_raw, never in the repo
rclone lsf --dirs-only --max-depth 2 "$SRC" | grep -i -E '^[^/]+/(uh|ug)30[0-9]{2}/$' > "$LIST"
echo "$(wc -l < "$LIST") study folders to copy"

while read -r d; do
    rclone copy --progress --exclude "**/[Ii]maging/**" "$SRC/$d" "$DEST/$d"
done < "$LIST"

echo "copy finished $(date)"

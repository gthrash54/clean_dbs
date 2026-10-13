#!/bin/bash
# copies the study folders of the cohorts you name from Box to $USER_SCRATCH/dbs_raw/<cohort>/, one cohort folder each
# (uh/UG are already in $USER_DATA/dbs_raw, see copy_uhug.sh). Box stays the master: this only ever reads from it.
#
# run it from the clean_dbs folder as a slurm job, so it keeps going after you log out (logs outside the repo):
#   mkdir -p $USER_DATA/copy_logs
#   sbatch --partition=short --time=12:00:00 --mem=2G --job-name=copy_dlpfc --output=$USER_DATA/copy_logs/%x_%j.out --wrap "bash scripts/copy_cohorts.sh dlpfc"
# several cohorts at once: bash scripts/copy_cohorts.sh pns tm tmscs
# options, before the cohort names:
#   --list-only   make the lists and print how many folders and files each cohort has, copy nothing
#   --plain-only  only number folders that are digits and dots (by default also ones with a short letter suffix, like <digits>_<letters>,
#                 which hold 49 study IDs that aren't anywhere else)
# settings, put them in front of the command (REMOTE=Box DEST=~/dbs_raw bash scripts/copy_cohorts.sh --list-only pns):
#   REMOTE  the rclone remote, "box" on cheaha (default), "Box" on my laptop
#   DEST    where the copies go, default $USER_SCRATCH/dbs_raw (scratch deletes files untouched for 90 days, the copy is only for cleaning)
#
# What it copies, per cohort, only from inside the de-identified number folders (never the patient-named folders at the top of Box):
#   - <number folder>/<study ID>/ for this cohort's IDs (DLPFC045, dlpfc045, SEP_03), whole, with any IDs nested inside them
#     (an LFP ID inside a PP ID comes with PP; clean_dbs files it under the inner ID)
#   - <number folder>/<ID>_<other ID>/ (two study IDs in one folder name) with the cohort of the first ID; the manifest lists those as excluded
#     until someone says whose they are
#   - files sitting straight in a number folder with this cohort's ID in their name (BPP012.vhdr)
#   - the folders with no study ID that belong to this cohort, by exact name (MOTIF, Motif protocol, PR, Reversal_Polarity, BurstPairPulses);
#     the manifest lists their recordings as excluded with the reason, nothing in them is cleaned
# Left out everywhere: Imaging/Images/DICOM/MRI folders and scan files (DICOM headers can carry patient names), documents, photos and
# videos other than tdt's .avi (consent forms, key files, pictures). Only rclone copy, never sync (sync deletes on the destination).
# Rerunning is safe, rclone copy skips files that are already there; then check the copy with scripts/check_copy.py.
# The log only gets totals once a minute (no file names: the loose files are listed with their number folder).
# The lists have Box number folders in them, so they stay next to the copy (never in the repo) and are never printed.

module load rclone 2>/dev/null || true # on cheaha rclone comes from a module (before set -u, module scripts trip over unset variables)
set -euo pipefail

REMOTE="${REMOTE:-box}"
SRC="$REMOTE:UAB DBS Data/UAB DBS Raw Data"
DEST="${DEST:-${USER_SCRATCH:?set USER_SCRATCH or DEST}/dbs_raw}"

# every cohort in clean_dbs/cohorts.py (tests/test_allcohorts.py checks this list against it)
COHORTS="uh ug dlpfc bi pp mc bv orbv bpp bb chr pt sep pns tm tmscs lfp joy ecog burst dys sensory pr motif sdc or"
# folders with no study ID, by exact name (any case), and the cohort they go with
NOID_motif="MOTIF|Motif protocol"
NOID_pr="PR|Reversal_Polarity"
NOID_burst="BurstPairPulses"

ANY="$(echo $COHORTS | tr ' ' '|')" # any cohort prefix, for the second ID in a folder named after two
LIST_ONLY=0
NUMBER='[0-9]+(\.[0-9]+)+(_[a-z]+)?' # a number folder: digits and dots, maybe a short letter suffix
while [ $# -gt 0 ]; do
    case "$1" in
        --list-only) LIST_ONLY=1; shift ;;
        --plain-only) NUMBER='[0-9]+(\.[0-9]+)+'; shift ;;
        -*) echo "unknown option $1" >&2; exit 2 ;;
        *) break ;;
    esac
done
[ $# -gt 0 ] || { echo "usage: bash scripts/copy_cohorts.sh [--list-only] [--plain-only] <cohort> [<cohort> ...]" >&2; exit 2; }
for c in "$@"; do
    c="$(echo "$c" | tr 'A-Z' 'a-z')"
    case " $COHORTS " in
        *" $c "*) ;;
        *) echo "unknown cohort $c, pick from: $COHORTS" >&2; exit 2 ;;
    esac
    if [ "$c" = uh ] || [ "$c" = ug ]; then
        echo "uh/UG are copied by scripts/copy_uhug.sh into \$USER_DATA/dbs_raw, not here" >&2; exit 2
    fi
done

# what not to copy, in any case: both imaging patterns are needed, "**/imaging/**" alone misses an Imaging folder right inside the study folder
FILTERS=(--ignore-case)
for d in imaging images dicom mri; do
    FILTERS+=(--exclude "$d/**" --exclude "**/$d/**")
done
for e in dcm nii nii.gz mgz pdf doc docx xls xlsx ppt pptx jpg jpeg png heic mp4 mov; do
    FILTERS+=(--exclude "*.$e")
done

mkdir -p "$DEST"
umask 077 # the lists are only for whoever runs the copy
# one listing of the top two levels for all cohorts, kept only until the lists below are made (it has what's in the patient-named folders too)
TOP="$(mktemp "$DEST/.box_top.XXXXXX")"
trap 'rm -f "$TOP"' EXIT
rclone lsf --max-depth 2 "$SRC" > "$TOP"

for c in "$@"; do
    c="$(echo "$c" | tr 'A-Z' 'a-z')"
    ID="$c[ _-]?[0-9]{2,5}" # same as cohorts.ID_RE for this prefix
    OUT="$DEST/$c"
    mkdir -p "$OUT"
    FOLDERS="$OUT/_box_folders.txt" # <number folder>/<folder>/ per line, never printed
    FILES="$OUT/_box_files.txt"     # <number folder>/<file> per line, never printed
    # study folders, and folders named after two IDs (first one this cohort's)
    grep -i -E "^$NUMBER/$ID([ _-]+($ANY)[ _-]?[0-9]{2,5})?/\$" "$TOP" > "$FOLDERS" || true
    # folders with no study ID that go with this cohort
    noid_var="NOID_$c"
    if [ -n "${!noid_var:-}" ]; then
        grep -i -E "^$NUMBER/(${!noid_var})/\$" "$TOP" >> "$FOLDERS" || true
    fi
    # files straight in a number folder with this cohort's ID in the name (not preceded by a letter: BPP012 isn't PP012), without the excluded types
    grep -i -E "^$NUMBER/([^/]*[^a-z/])?$ID([^0-9/][^/]*)?\$" "$TOP" \
        | grep -i -v -E '\.(dcm|nii|nii\.gz|mgz|pdf|docx?|xlsx?|pptx?|jpe?g|png|heic|mp4|mov)$' > "$FILES" || true
    echo "$c: $(wc -l < "$FOLDERS") folders and $(wc -l < "$FILES") loose files to copy"
done
rm -f "$TOP" # the lists are made, the full listing isn't needed any more
[ "$LIST_ONLY" = 1 ] && exit 0

for c in "$@"; do
    c="$(echo "$c" | tr 'A-Z' 'a-z')"
    OUT="$DEST/$c"
    while read -r d; do
        rclone copy --stats 1m --stats-one-line "${FILTERS[@]}" "$SRC/$d" "$OUT/$d"
    done < "$OUT/_box_folders.txt"
    if [ -s "$OUT/_box_files.txt" ]; then
        # --files-from-raw can't be mixed with other filters, the grep above already left out the excluded types
        rclone copy --stats 1m --stats-one-line --files-from-raw "$OUT/_box_files.txt" "$SRC" "$OUT"
    fi
    echo "$c: copy finished $(date)"
done

#!/bin/bash
# copy_preproc_data.sh  (rewritten 2026-10-17: PreprocData was lost on the grid)
#
# Copies each subject's preprocessed encoding data from the network share
# to BEBO2/PreprocData/<subj>_Enc_allcond.mat (flat, 2026-10-17) -- the path every
# subject-level step reads (pdat.rawEnc in config_and_run_RSA_subject.m,
# config_and_run_alpha_pre2post_subject.m, config_and_run_ebosc_subject.m,
# config_and_run_RSA_singleTrial_merged.m).
#
# Subjects: subj_list_all93.txt (same list as submit_BEBO_pipeline.sh).
# Files that are already complete at the destination are skipped (rsync),
# so the script can simply be run again after an interruption.
#
# Run from inside RSA_pipeline/:
#   bash copy_preproc_data.sh check   only CHECKS: prints one line with how
#                                     many source files exist and how many
#                                     are already at the destination
#   bash copy_preproc_data.sh         copies (one summary line at the end;
#                                     missing subjects -> ~/logs/copy_preproc_missing.txt)

cd "$(dirname "$0")" || exit 1

SRC_ROOT="$HOME/NetworkFolders/private/EEG/EEG-Data/Analyses"
DEST_ROOT="$HOME/BEBO2/PreprocData"
SUBJLIST="subj_list_all93.txt"
MODE=${1:-copy}

if [ ! -f "$SUBJLIST" ]; then echo "Subject list $SUBJLIST not found in $(pwd)"; exit 1; fi
if [ ! -d "$SRC_ROOT" ]; then echo "Source folder not found: $SRC_ROOT (is the network share mounted?)"; exit 1; fi

mkdir -p "$DEST_ROOT" "$HOME/logs"
missing_file="$HOME/logs/copy_preproc_missing.txt"
: > "$missing_file"

n_total=0; n_src=0; n_dest=0; n_ok=0; n_fail=0
while read -r subj; do
    [ -z "$subj" ] && continue
    n_total=$((n_total+1))
    src="$SRC_ROOT/$subj/eeg/${subj}_Enc_allcond.mat"
    dest_dir="$DEST_ROOT"
    dest="$dest_dir/${subj}_Enc_allcond.mat"

    [ -f "$dest" ] && n_dest=$((n_dest+1))
    if [ ! -f "$src" ]; then
        echo "$subj MISSING SOURCE $src" >> "$missing_file"
        continue
    fi
    n_src=$((n_src+1))
    [ "$MODE" = check ] && continue

    mkdir -p "$dest_dir"
    if rsync -a "$src" "$dest"; then
        n_ok=$((n_ok+1))
    else
        echo "$subj COPY FAILED" >> "$missing_file"
        n_fail=$((n_fail+1))
    fi
done < "$SUBJLIST"

if [ "$MODE" = check ]; then
    echo "CHECK: $n_total subjects, $n_src source files found, $n_dest already at destination"
else
    echo "DONE: $n_total subjects, $n_ok copied/up to date, $((n_total-n_src)) missing source, $n_fail failed"
fi

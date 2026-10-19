#!/bin/bash
# submit_BEBO_rerun_from_perAge.sh
#
# Reruns Ispec from the per-age group step onwards, using the
# existing subject-level similarity data (Results/Similarity_Data/):
#
#   Ispec:   per-age YA/OA -> fullSample -> alpha ROI -> single-trial -> aggregation
#            -> plots
#
# Does NOT touch alpha_pre2post / PF or the subject-level jobs.
#
# Thin wrapper: all job definitions live in submit_BEBO_pipeline.sh (so
# the two can never drift apart); this just calls it with start point
# 'perAge'. To rerun the whole pipeline instead: bash submit_BEBO_pipeline.sh
#
# Run with: bash submit_BEBO_rerun_from_perAge.sh   (from inside RSA_pipeline/)

cd "$(dirname "$0")" || exit 1
bash submit_BEBO_pipeline.sh perAge

#!/bin/bash
# submit_eBOSC_rerun.sh
#
# Runs ONLY the eBOSC branch: eBOSC subject array -> results + stats ->
# plots (2026-10-05: all episodes saved; alpha 8-14 and beta 14-25 Hz
# rhythms, ebosc_bands.m). Assumes the Ispec fullSample and Ispec ROI jobs
# have finished (eBOSC results use the Ispec cluster channels and Ispec).
# Thin wrapper around submit_BEBO_pipeline.sh (start point 'ebosc').
# Run with: bash submit_eBOSC_rerun.sh   (from inside RSA_pipeline/)

cd "$(dirname "$0")" || exit 1
bash submit_BEBO_pipeline.sh ebosc

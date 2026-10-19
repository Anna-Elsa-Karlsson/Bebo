#!/bin/bash
# submit_eBOSC_results_plots.sh
#
# Reruns ONLY the eBOSC results (stats, per-band measures, PF topographies)
# and the eBOSC figure -- not the eBOSC subject array. Use after changes to
# ebosc_results.m / ebosc_plots.m. Assumes the eBOSC subject array (current
# version: all episodes saved) and the Ispec ROI step have finished.
# Thin wrapper around submit_BEBO_pipeline.sh (start point 'eboscResults').
# Run with: bash submit_eBOSC_results_plots.sh   (from inside RSA_pipeline/)

cd "$(dirname "$0")" || exit 1
bash submit_BEBO_pipeline.sh eboscResults

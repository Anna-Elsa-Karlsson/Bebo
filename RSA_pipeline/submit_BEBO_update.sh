#!/bin/bash
# submit_BEBO_update.sh
#
# Reruns everything affected by the changes of 2026-10-02..05, WITHOUT
# recomputing any cluster test (Ispec fullSample and alpha group stay):
#   between-category control pairs (subject array) + Ispec ROI (whole
#   cluster) -> single-trial (alpha 8-20 Hz; alpha of 1st/2nd/average) ->
#   aggregation -> Ispec figures + stats; alpha x pair type (results + figure)
#   alpha figure + stats (desync in the desync cluster)
#   eBOSC subject array (all episodes) -> eBOSC results (alpha/beta) -> figure
# Thin wrapper around submit_BEBO_pipeline.sh (start point 'update').
# Run with: bash submit_BEBO_update.sh   (from inside RSA_pipeline/)

cd "$(dirname "$0")" || exit 1
bash submit_BEBO_pipeline.sh update

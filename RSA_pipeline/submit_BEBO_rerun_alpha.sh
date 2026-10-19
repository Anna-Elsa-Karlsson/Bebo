#!/bin/bash
# submit_BEBO_rerun_alpha.sh
#
# Reruns everything that uses alpha (2026-10-02: alpha cluster test with
# minnbchan 2; alpha desync stats and eBOSC results now in the Ispec ROI):
#   alpha subject array -> alpha group cluster test
#   Ispec ROI (core channels, from the existing Ispec cluster stats) ->
#     alpha figure + stats (desync in the Ispec ROI)
#     single-trial extraction -> aggregation -> Ispec plots + stats
#     eBOSC results (Ispec ROI channels) -> eBOSC plots
# NOT recomputed: the Ispec RSA, its per-age and fullSample stats, and the
# eBOSC subject array. Submit only AFTER the eBOSC subject array has
# finished (eBOSC results read its per-subject files).
# Thin wrapper around submit_BEBO_pipeline.sh (start point 'alpha').
# Run with: bash submit_BEBO_rerun_alpha.sh   (from inside RSA_pipeline/)

cd "$(dirname "$0")" || exit 1
bash submit_BEBO_pipeline.sh alpha

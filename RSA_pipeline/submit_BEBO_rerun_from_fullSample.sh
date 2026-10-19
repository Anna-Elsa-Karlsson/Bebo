#!/bin/bash
# submit_BEBO_rerun_from_fullSample.sh
#
# Reruns everything from the Ispec 2nd-level (fullSample) cluster stats on,
# WITHOUT recomputing the subject-level similarity matrices, the per-age
# jobs, or the alpha / eBOSC subject steps:
#   Ispec fullSample cluster test -> Ispec ROI (whole cluster) ->
#     single-trial extraction -> aggregation -> Ispec plots + stats
#     alpha figure + stats (desync in the Ispec cluster)
#     eBOSC results + plots (Ispec cluster channels)
# Assumes the per-age jobs, the alpha group job and the eBOSC subject
# array have finished.
# Thin wrapper around submit_BEBO_pipeline.sh (start point 'fullSample').
# Run with: bash submit_BEBO_rerun_from_fullSample.sh   (from inside RSA_pipeline/)

cd "$(dirname "$0")" || exit 1
bash submit_BEBO_pipeline.sh fullSample

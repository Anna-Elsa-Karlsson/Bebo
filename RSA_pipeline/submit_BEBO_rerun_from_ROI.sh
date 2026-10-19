#!/bin/bash
# submit_BEBO_rerun_from_ROI.sh
#
# Reruns everything AFTER the Ispec 2nd-level cluster stats (which are not
# recomputed):
#   Ispec ROI (whole cluster) -> single-trial extraction -> aggregation ->
#     Ispec plots + stats
#   alpha figure + stats (needs the Ispec ROI)
#   eBOSC results (Ispec ROI channels) -> eBOSC plots
# NOT included: the eBOSC subject array (run it with submit_eBOSC_rerun.sh).
# Assumes Ispec fullSample, the alpha group job and the eBOSC subject
# array have finished.
# Thin wrapper around submit_BEBO_pipeline.sh (start point 'roi').
# Run with: bash submit_BEBO_rerun_from_ROI.sh   (from inside RSA_pipeline/)

cd "$(dirname "$0")" || exit 1
bash submit_BEBO_pipeline.sh roi

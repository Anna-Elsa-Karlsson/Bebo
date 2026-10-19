#!/bin/bash
# submit_PF_rerun.sh
#
# Reruns ONLY the PF chain (PF_derive_roi -> PF subject array -> PF group)
# using the existing alpha_pre2post group results -- e.g. after changing
# the PF search band (alpha_band.m) or ROI rule.
# Thin wrapper around submit_BEBO_pipeline.sh (start point 'pf').
# Run with: bash submit_PF_rerun.sh   (from inside RSA_pipeline/)

cd "$(dirname "$0")" || exit 1
bash submit_BEBO_pipeline.sh pf

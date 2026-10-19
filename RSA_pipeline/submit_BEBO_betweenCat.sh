#!/bin/bash
# submit_BEBO_betweenCat.sh
#
# Between-category CONTROL pairs for the alpha x pair-type analysis
# (2026-10-05): lag-matched exemplar 1 of category A vs exemplar 2 of a
# different category B, same TFR files / background removal / baseline
# method / Ispec cluster as Isim and Csim. Runs:
#   subject array (Bsim pairs) -> single-trial extraction -> aggregation
#   -> Ispec plots + alpha x pair type (results + figure)
# Assumes the Ispec fullSample stats and the Ispec ROI step have finished.
# Thin wrapper around submit_BEBO_pipeline.sh (start point 'bsim').
# Run with: bash submit_BEBO_betweenCat.sh   (from inside RSA_pipeline/)

cd "$(dirname "$0")" || exit 1
bash submit_BEBO_pipeline.sh bsim

#!/bin/bash
# submit_BEBO_2_afterSubjectRSA.sh
#
# Step 2 of 2: everything after the subject-level RSA, assuming step 1
# (submit_BEBO_1_subjectRSA.sh) has finished:
#
#   alpha/PF: alpha_pre2post subject array -> group (cluster test) ->
#             PF_derive_roi -> PF subject array -> PF group; alpha figure
#   behaviour: figure + IM/LDI stats (Plots_Behavior/)
#   Ispec:    per-age YA/OA -> fullSample -> alpha ROI -> single-trial -> aggregation
#             -> plots
#
# Thin wrapper around submit_BEBO_pipeline.sh (start point 'afterRsaSubj').
# Run with: bash submit_BEBO_2_afterSubjectRSA.sh   (from inside RSA_pipeline/)

cd "$(dirname "$0")" || exit 1
bash submit_BEBO_pipeline.sh afterRsaSubj

#!/bin/bash
# submit_BEBO_1_subjectRSA.sh
#
# Step 1 of 2: ONLY the Ispec subject-level RSA array (one task per
# subject in subj_list_all93.txt): loads the original TFR files, BOSC,
# pairing (Isim + Csim), similarity + baseline correction -> 
# Results/Similarity_Data/Ispec_30to100Hz/<ID>_within-item/-cat_sim_ST.mat
#
# When it has finished: bash submit_BEBO_2_afterSubjectRSA.sh
#
# Thin wrapper around submit_BEBO_pipeline.sh (start point 'rsaSubj').
# Run with: bash submit_BEBO_1_subjectRSA.sh   (from inside RSA_pipeline/)

cd "$(dirname "$0")" || exit 1
bash submit_BEBO_pipeline.sh rsaSubj

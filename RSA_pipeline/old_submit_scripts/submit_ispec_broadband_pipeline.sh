#!/bin/bash
# submit_ispec_broadband_pipeline.sh
#
# THE submission script for this pipeline -- runs the complete,
# validated chain for Ispec broadband (30-100Hz), the only variant this
# project uses going forward. Replaces all earlier fragmented/variant-
# generic submission scripts (submit_rsa_all_variants.sh,
# submit_rsa_downstream.sh, submit_subject_level.sh,
# submit_alpha_singleTrial.sh, submit_singleTrial_merged.sh,
# submit_plots_ispec_broadband.sh -- all removed).
#
# Chain (fully dependency-chained, submit once and let SLURM handle the
# ordering):
#   1. subject-level array job (93 tasks) -- per-subject similarity
#   2. per-age jobs (YA, OA) -- grand average + 1st-level within-person
#      test (independent_samples_ttest: all raw item-sim pairs vs all raw
#      cat-sim pairs, unpaired -- the TRUE method the manuscript is based
#      on, restored 2026-09-XX after a category-level paired alternative
#      was tried and came out empirically weaker -- see
#      rsa_group_level_perAge.m headers for the full story)
#   3. fullSample job -- 2nd-level cluster stat (500 permutations) +
#      3rd-level extraction
#
# Runs the PRIMARY (baseline-corrected, simField='powspctrm') analysis
# AND, automatically and in parallel, the UNCORRECTED comparison
# (simField='powspctrmOrig', to see what the baseline correction actually
# does) -- steps 2b/3b below. Both share the SAME subject-level output
# (step 1 already saves both fields, so it only runs once) and are fully
# independent of each other (2b/3b depend only on step 1, not on
# 2/3) -- SLURM runs them concurrently. Uncorrected outputs are tagged
# with an '_uncorrected' suffix throughout, so nothing overwrites the
# primary results. Steps 4/5/6 (plots, single-trial merge, aggregation)
# are downstream of the PRIMARY chain only -- the uncorrected comparison
# is a standalone diagnostic, not part of that path.
#
# The single-trial glmer comparison (Sim vs Sim_noBL) needs no separate
# submission at all -- both columns already exist in the merged
# single-trial tables from step 5/6; just refit the R models on Sim_noBL.
#   4. plots job -- Fig3/Fig4-style figures (PRIMARY/corrected chain only)
#   5. single-trial merged array job (93 tasks) -- per-subject alpha +
#      similarity extraction (rsa_singleTrial_merged.m -- the validated,
#      working version; faithful merge of the original singleTrl_alpha.m
#      + self_sim_singlTrl_fullSample.m logic, mask NaN-propagation bug
#      fixed)
#   6. aggregation job -- concatenates the 93 single-trial merged tables
#      into the final combined within-item/within-cat tables
#
# Time budgets reflect what was actually needed during development (see
# session notes) -- per-age/fullSample/plots all previously timed out at
# shorter guesses, so are set generously here.
#
# Run with: bash submit_ispec_broadband_pipeline.sh   (from inside RSA_pipeline/)

FREQLO=30
FREQHI=100
CTYPE=Ispec
BGMODE=legacy
SUBJLIST="subj_list_all93.txt" # TEMPORARILY reverted (2026-09-XX) to test replication on the original 93 before re-adding the 3 new IDs
NSUBJ=$(wc -l < "$SUBJLIST")
mkdir -p "$HOME/logs"

echo "=== Submitting Ispec_broadband (${FREQLO}-${FREQHI}Hz) full pipeline ==="

# --- 1. Subject-level array job (93 tasks) ---
cat > tmp_subj.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%A_%a.out
#SBATCH --cpus-per-task 2
#SBATCH --time 6:0:0
#SBATCH --mem 60gb
#SBATCH --array=1-${NSUBJ}
#SBATCH --job-name RSA_IspecBB_subj
# Time/mem bumped 2026-09-04: this job now runs the raw DPSS gamma TFR
# decomposition inline (from Enc_allcond, cached to disk) instead of
# just loading a pre-made file -- generous default, not yet timed.
module load matlab/R2017b
SUBJ=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $SUBJLIST)
matlab -r "config_and_run_RSA_subject('\$SUBJ','$CTYPE','$FREQLO','$FREQHI','$BGMODE'); exit"
EOF
SUBJ_JOBID=$(sbatch --parsable tmp_subj.batch)
echo "  subject-level array job ($NSUBJ tasks): $SUBJ_JOBID"

# --- 2. Per-age jobs (YA, OA) -- depend on the full array finishing ---
PERAGE_JOBIDS=()
for GROUP in YA OA; do
  cat > tmp_group.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 2
#SBATCH --time 22:0:0
#SBATCH --mem 60gb
#SBATCH --job-name RSA_IspecBB_${GROUP}
#SBATCH --dependency=afterok:${SUBJ_JOBID}
module load matlab/R2017b
matlab -r "config_and_run_RSA_group_perAge('$GROUP','$CTYPE','$FREQLO','$FREQHI'); exit"
EOF
  JID=$(sbatch --parsable tmp_group.batch)
  PERAGE_JOBIDS+=("$JID")
  echo "  per-age ($GROUP) job: $JID"
done

# --- 3. fullSample job -- depends on BOTH per-age jobs ---
DEP=$(IFS=:; echo "${PERAGE_JOBIDS[*]}")
cat > tmp_full.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 22:0:0
#SBATCH --mem 85gb
#SBATCH --job-name RSA_IspecBB_full
#SBATCH --dependency=afterok:${DEP}
module load matlab/R2017b
matlab -r "config_and_run_RSA_group_fullSample('$CTYPE','$FREQLO','$FREQHI'); exit"
EOF
FULL_JOBID=$(sbatch --parsable tmp_full.batch)
echo "  fullSample job (500 permutations): $FULL_JOBID"

# --- 2b/3b. UNCORRECTED comparison (simField='powspctrmOrig') -- runs
#     independently and in parallel with 2/3, depending only on the
#     subject-level array (step 1), which already saved both fields ---
PERAGE_UNCORR_JOBIDS=()
for GROUP in YA OA; do
  cat > tmp_group_uncorr.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 2
#SBATCH --time 22:0:0
#SBATCH --mem 60gb
#SBATCH --job-name RSA_IspecBB_${GROUP}_uncorrected
#SBATCH --dependency=afterok:${SUBJ_JOBID}
module load matlab/R2017b
matlab -r "config_and_run_RSA_group_perAge('$GROUP','$CTYPE','$FREQLO','$FREQHI','powspctrmOrig'); exit"
EOF
  JID=$(sbatch --parsable tmp_group_uncorr.batch)
  PERAGE_UNCORR_JOBIDS+=("$JID")
  echo "  per-age ($GROUP, uncorrected) job: $JID"
done

DEP_UNCORR=$(IFS=:; echo "${PERAGE_JOBIDS[*]} ${PERAGE_UNCORR_JOBIDS[*]}" | tr ' ' ':')
# Depends on BOTH the corrected AND uncorrected per-age jobs: the GA files
# (including sim_GA_orig.mat) are only ever written by the CORRECTED
# per-age job (see rsa_group_level_perAge.m's simField guard -- GA
# computation is simField-independent, so it only runs once, on the
# corrected call, to avoid two parallel jobs racing on the same files),
# while the uncorrected 1st-level stat file comes from the uncorrected
# per-age job. This fullSample (uncorrected) run needs both.
cat > tmp_full_uncorr.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 22:0:0
#SBATCH --mem 85gb
#SBATCH --job-name RSA_IspecBB_full_uncorrected
#SBATCH --dependency=afterok:${DEP_UNCORR}
module load matlab/R2017b
matlab -r "config_and_run_RSA_group_fullSample('$CTYPE','$FREQLO','$FREQHI','powspctrmOrig'); exit"
EOF
FULL_UNCORR_JOBID=$(sbatch --parsable tmp_full_uncorr.batch)
echo "  fullSample job (uncorrected, 500 permutations): $FULL_UNCORR_JOBID"

# --- 4. plots job -- depends on fullSample (PRIMARY/corrected chain only) ---
cat > tmp_plots.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 22:0:0
#SBATCH --mem 20gb
#SBATCH --job-name RSA_IspecBB_plots
#SBATCH --dependency=afterok:${FULL_JOBID}
module load matlab/R2017b
matlab -r "config_and_run_RSA_plots('$CTYPE','$FREQLO','$FREQHI'); exit"
EOF
PLOTS_JOBID=$(sbatch --parsable tmp_plots.batch)
echo "  plots job: $PLOTS_JOBID"

# --- 4b. Ispec_derive_alphaROI -- derives the non-circular alpha
#         extraction ROI (channels + time window) from the fullSample
#         2nd-level cluster stat. Depends on FULL_JOBID only. ---
cat > tmp_alpharoi.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 0:15:0
#SBATCH --mem 10gb
#SBATCH --job-name RSA_IspecBB_alphaROI
#SBATCH --dependency=afterok:${FULL_JOBID}
module load matlab/R2017b
matlab -r "config_and_run_Ispec_derive_alphaROI('$FREQLO','$FREQHI'); exit"
EOF
ALPHAROI_JOBID=$(sbatch --parsable tmp_alpharoi.batch)
echo "  Ispec_derive_alphaROI job: $ALPHAROI_JOBID"

# --- 5. Single-trial merged array job (93 tasks) -- depends on the
#        subject-level array, the fullSample job (needs the cluster
#        mask), AND the derived alpha ROI ---
cat > tmp_stmerged.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%A_%a.out
#SBATCH --cpus-per-task 2
#SBATCH --time 2:0:0
#SBATCH --mem 40gb
#SBATCH --array=1-${NSUBJ}
#SBATCH --job-name RSA_IspecBB_singleTrialMerged
#SBATCH --dependency=afterok:${SUBJ_JOBID}:${FULL_JOBID}:${ALPHAROI_JOBID}
module load matlab/R2017b
SUBJ=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $SUBJLIST)
matlab -r "config_and_run_RSA_singleTrial_merged('\$SUBJ','$FREQLO','$FREQHI'); exit"
EOF
ST_JOBID=$(sbatch --parsable tmp_stmerged.batch)
echo "  single-trial merged array job ($NSUBJ tasks): $ST_JOBID"

# --- 6. Aggregation job -- depends on the single-trial merged array ---
cat > tmp_stagg.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 0:30:0
#SBATCH --mem 20gb
#SBATCH --job-name RSA_IspecBB_singleTrialAggregate
#SBATCH --dependency=afterok:${ST_JOBID}
module load matlab/R2017b
matlab -r "aggregate_singleTrial_merged($FREQLO,$FREQHI); exit"
EOF
AGG_JOBID=$(sbatch --parsable tmp_stagg.batch)
echo "  single-trial aggregation job: $AGG_JOBID"

rm -f tmp_subj.batch tmp_group.batch tmp_full.batch tmp_group_uncorr.batch tmp_full_uncorr.batch tmp_plots.batch tmp_alpharoi.batch tmp_stmerged.batch tmp_stagg.batch
echo ""
echo "Full Ispec_broadband pipeline submitted (corrected + uncorrected comparison). Check progress with: squeue -u \$USER"

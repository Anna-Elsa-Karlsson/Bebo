#!/bin/bash
# submit_ispec_lowhigh_gamma_comparison.sh
#
# Exploratory band-specificity check (2026-09-05): does the Ispec effect
# hold up in low gamma (30-50Hz) vs high gamma (60-100Hz) separately?
# Motivated by two convergent (if theoretically distinct) literatures --
# Colgin et al. 2009's slow/fast gamma retrieval/encoding dissociation,
# and the ~40Hz/STDP-plasticity encoding literature (Fell & Axmacher
# 2011) -- both implicate LOW/slow gamma as the more likely locus of a
# genuine memory-related effect, regardless of which mechanism is
# invoked, making low-vs-high a real, motivated hypothesis test.
#
# LEAN chain, reusing the EXISTING cached 30-100Hz raw gamma TFR for
# free (bgFitMode='legacy', sub-band fully inside the cached range -- no
# new decomposition needed, see rsa_subject_level.m's
# load_bosc_correct_restrict_rawGamma header):
#   1. subject-level array job (93 tasks) -- per band
#   2. per-age jobs (YA, OA) -- per band
#   3. fullSample job -- 2nd-level cluster stat + 3rd-level extraction, per band
# Deliberately does NOT include plots, single-trial merge, aggregation,
# or the corrected/uncorrected comparison -- this is a quick
# band-specificity check on the primary (corrected) analysis only, not a
# second full reportable pipeline.
#
# Run with: bash submit_ispec_lowhigh_gamma_comparison.sh   (from inside RSA_pipeline/)

CTYPE=Ispec
BGMODE=legacy
SUBJLIST="subj_list_all93.txt" # TEMPORARILY reverted (2026-09-XX) to test replication on the original 93 before re-adding the 3 new IDs
NSUBJ=$(wc -l < "$SUBJLIST")
mkdir -p "$HOME/logs"

for BAND in "30 50" "60 100"; do
  read -r FREQLO FREQHI <<< "$BAND"
  LABEL="${FREQLO}to${FREQHI}Hz"
  echo "=== Submitting Ispec ${LABEL} (band-specificity check) ==="

  # --- 1. Subject-level array job (93 tasks) -- reuses the cached 30-100Hz
  #        raw TFR, so this should be much faster than the broadband run's
  #        first pass (no decomposition, just BOSC + correlation on a
  #        frequency subset) ---
  cat > "tmp_subj_${LABEL}.batch" << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%A_%a.out
#SBATCH --cpus-per-task 2
#SBATCH --time 2:0:0
#SBATCH --mem 40gb
#SBATCH --array=1-${NSUBJ}
#SBATCH --job-name RSA_Ispec_${LABEL}_subj
module load matlab/R2017b
SUBJ=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $SUBJLIST)
matlab -r "config_and_run_RSA_subject('\$SUBJ','$CTYPE','$FREQLO','$FREQHI','$BGMODE'); exit"
EOF
  SUBJ_JOBID=$(sbatch --parsable "tmp_subj_${LABEL}.batch")
  echo "  subject-level array job ($NSUBJ tasks): $SUBJ_JOBID"

  # --- 2. Per-age jobs (YA, OA) -- depend on the full array finishing ---
  PERAGE_JOBIDS=()
  for GROUP in YA OA; do
    cat > "tmp_group_${LABEL}.batch" << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 2
#SBATCH --time 22:0:0
#SBATCH --mem 60gb
#SBATCH --job-name RSA_Ispec_${LABEL}_${GROUP}
#SBATCH --dependency=afterok:${SUBJ_JOBID}
module load matlab/R2017b
matlab -r "config_and_run_RSA_group_perAge('$GROUP','$CTYPE','$FREQLO','$FREQHI'); exit"
EOF
    JID=$(sbatch --parsable "tmp_group_${LABEL}.batch")
    PERAGE_JOBIDS+=("$JID")
    echo "  per-age ($GROUP) job: $JID"
  done

  # --- 3. fullSample job -- depends on BOTH per-age jobs ---
  DEP=$(IFS=:; echo "${PERAGE_JOBIDS[*]}")
  cat > "tmp_full_${LABEL}.batch" << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 22:0:0
#SBATCH --mem 85gb
#SBATCH --job-name RSA_Ispec_${LABEL}_full
#SBATCH --dependency=afterok:${DEP}
module load matlab/R2017b
matlab -r "config_and_run_RSA_group_fullSample('$CTYPE','$FREQLO','$FREQHI'); exit"
EOF
  FULL_JOBID=$(sbatch --parsable "tmp_full_${LABEL}.batch")
  echo "  fullSample job (500 permutations): $FULL_JOBID"

  rm -f "tmp_subj_${LABEL}.batch" "tmp_group_${LABEL}.batch" "tmp_full_${LABEL}.batch"
  echo ""
done

echo "Low- and high-gamma band-specificity chains submitted. Check progress with: squeue -u \$USER"
echo "Results will land in Results/Stats/Ispec_30to50Hz/ and Results/Stats/Ispec_60to100Hz/"

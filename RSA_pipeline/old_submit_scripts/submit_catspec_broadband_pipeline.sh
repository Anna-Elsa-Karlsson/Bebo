#!/bin/bash
# submit_catspec_broadband_pipeline.sh
#
# Catspec (within-cat vs between-cat) broadband gamma (30-100Hz) pipeline
# -- anticipated reviewer question alongside the Ispec (within-item vs
# within-cat) result. Reuses the same rsa_subject_level.m/
# config_and_run_RSA_group_perAge.m/config_and_run_RSA_group_fullSample.m
# machinery as Ispec, just with contrastType='Catspec'.
#
# The within-cat side automatically inherits the current, final pairing
# rule (exemplar duos (1,2) and (3,4), OR-validity each -- see
# pairs_within_cat_exemplar1v2), since Catspec's within-cat and Ispec's
# cat-sim share that exact function. The between-cat side is a separate,
# untouched block in pairs_within_vs_between_cat (exemplar 1 vs exemplar
# 2 of a DIFFERENT category, both sides require their own valid response
# -- not part of any of the recent pairing discussion).
#
# COMPUTE NOTE: the raw gamma TFR cache (pdat.tfrOut, keyed by subject +
# frequency band only, not by contrastType) was already populated by the
# Ispec run for these same subjects/band -- this pipeline's subject-level
# step should hit that cache and skip the expensive DPSS decomposition
# entirely, just rerunning BOSC correction + pairing/correlation. Time
# budget below is set lower than Ispec's original run accordingly; bump
# it if the cache turns out to be stale/missing for any subject.
#
# Chain (subject-level -> per-age -> fullSample only -- no plots, no
# alpha ROI derivation, no single-trial merge; those are specific to the
# Isim/Csim joint-model work, not needed for this standalone contrast):
#   1. subject-level array job
#   2. per-age jobs (YA, OA)
#   3. fullSample job (2nd-level cluster stat + 3rd-level extraction)
#
# Currently set to subj_list_all93.txt to match the ACTIVE state of
# ids_ALL.m (93-subject replication test, 3111/3213/4210 temporarily
# removed) -- update to subj_list_all96.txt alongside Ispec's own move
# back to 96, don't let the two drift apart.
#
# Run with: bash submit_catspec_broadband_pipeline.sh   (from inside RSA_pipeline/)

FREQLO=30
FREQHI=100
CTYPE=Catspec
BGMODE=legacy
SUBJLIST="subj_list_all93.txt"
NSUBJ=$(wc -l < "$SUBJLIST")
mkdir -p "$HOME/logs"

echo "=== Submitting Catspec_broadband (${FREQLO}-${FREQHI}Hz) pipeline ($NSUBJ subjects) ==="

# --- 1. Subject-level array job ---
cat > tmp_catspec_subj.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%A_%a.out
#SBATCH --cpus-per-task 2
#SBATCH --time 2:0:0
#SBATCH --mem 40gb
#SBATCH --array=1-${NSUBJ}
#SBATCH --job-name RSA_CatspecBB_subj
module load matlab/R2017b
SUBJ=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $SUBJLIST)
matlab -r "config_and_run_RSA_subject('\$SUBJ','$CTYPE','$FREQLO','$FREQHI','$BGMODE'); exit"
EOF
SUBJ_JOBID=$(sbatch --parsable tmp_catspec_subj.batch)
echo "  subject-level array job ($NSUBJ tasks): $SUBJ_JOBID"

# --- 2. Per-age jobs (YA, OA) -- depend on the full array finishing ---
PERAGE_JOBIDS=()
for GROUP in YA OA; do
  cat > tmp_catspec_group.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 2
#SBATCH --time 22:0:0
#SBATCH --mem 60gb
#SBATCH --job-name RSA_CatspecBB_${GROUP}
#SBATCH --dependency=afterok:${SUBJ_JOBID}
module load matlab/R2017b
matlab -r "config_and_run_RSA_group_perAge('$GROUP','$CTYPE','$FREQLO','$FREQHI'); exit"
EOF
  JID=$(sbatch --parsable tmp_catspec_group.batch)
  PERAGE_JOBIDS+=("$JID")
  echo "  per-age ($GROUP) job: $JID"
done

# --- 3. fullSample job -- depends on BOTH per-age jobs ---
DEP=$(IFS=:; echo "${PERAGE_JOBIDS[*]}")
cat > tmp_catspec_full.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 22:0:0
#SBATCH --mem 85gb
#SBATCH --job-name RSA_CatspecBB_full
#SBATCH --dependency=afterok:${DEP}
module load matlab/R2017b
matlab -r "config_and_run_RSA_group_fullSample('$CTYPE','$FREQLO','$FREQHI'); exit"
EOF
FULL_JOBID=$(sbatch --parsable tmp_catspec_full.batch)
echo "  fullSample job (500 permutations): $FULL_JOBID"

rm -f tmp_catspec_subj.batch tmp_catspec_group.batch tmp_catspec_full.batch
echo ""
echo "Catspec_broadband pipeline submitted. Check progress with: squeue -u \$USER"

#!/bin/bash
# submit_BEBO_pipeline.sh
#
# THE submission script for the whole BEBO EEG pipeline (2026-09-26) --
# replaces submit_ispec_broadband_pipeline.sh,
# submit_catspec_broadband_pipeline.sh, submit_alpha_pre2post_pipeline.sh,
# submit_ispec_from_fullsample.sh, submit_ispec_after_fullsample.sh and
# submit_ispec_lowhigh_gamma_comparison.sh (dropped -- never showed
# anything). Old scripts moved to old_submit_scripts/. Catspec (within-
# vs between-category) removed from the chain 2026-09-29 -- never showed
# an effect; its code (contrastType 'Catspec') is still in the MATLAB
# scripts if it's ever needed again.
#
# Everything is submitted at once; SLURM dependencies start each job as
# soon as the jobs it needs have finished (afterok = only if they
# succeeded). Two independent branches run in parallel:
#
#   A. Alpha_pre2post (alpha desync existence test + cluster)
#      A1 subject array -> A2 group (depsamplesT vs 0)
#      A6 alpha figure + stats, desync in the desync cluster   <- A2
#
#   E. eBOSC alpha rhythm (replaces the FFT-based PF branch, 2026-09-30)
#      E1 eBOSC subject array (all channels, no dependencies)
#      E2 eBOSC results: Ispec ROI channels, PF/abundance/
#         amplitude, stats incl. PF x Ispec                    <- E1, B4
#      E3 eBOSC plots (plotting only)                          <- E2
#
#   D. Behaviour: D1 figure + IM/LDI stats (Plots_Behavior/), no
#      dependencies (runs with 'all' and 'afterRsaSubj')
#
#   B. Ispec broadband gamma (30-100Hz), within-item vs within-cat
#      B1 subject array (similarity; raw gamma TFR cached per subject)
#      BS between-category control pairs, subject array (2026-10-05)
#      B2 per-age YA/OA (grand average + 1st-level)          <- B1
#      B3 fullSample (2nd-level cluster + 3rd-level means)   <- B2
#      B4 Ispec_derive_alphaROI                              <- B3
#      B5 single-trial array (T_item/T_cat/T_joint + alpha)  <- B1,B3,B4
#      B6 aggregation                                        <- B5
#      B7 plots (needs aggregated per-subject alpha, and the  <- B6, A1
#         alpha_pre2post per-subject files for the BL-alpha comparison)
#      B8 alpha x pair type, results (pair table from B6)    <- B6
#      B9 alpha x pair type, figure                          <- B8
#      (The uncorrected comparison, simField='powspctrmOrig', was removed
#      from the chain 2026-09-30 -- the MATLAB code still supports it.)
#
# Resources (2026-09-29): gamma TFR back to the original 8 ms grid ->
# similarity matrices 213x213 per channel (was 86x86), ~6x compute and
# ~10 GB per subject/condition sim file. Subject arrays 22h/100gb,
# per-age and fullSample 22h (max used), single-trial 6h/80gb.
#
# Trial rule used by every analysis: exemplars 1-2 only (Isim: exemplar
# 1-2 items; Csim: exemplar 1 vs 2), presentations 1-2, no 'new' responses.
#
# SUBJLIST must match ids_ALL.m (group-level steps read ids_ALL('all')).
#
# Run with (from inside RSA_pipeline/):
#   bash submit_BEBO_pipeline.sh              everything (A, B, C)
#   bash submit_BEBO_pipeline.sh rsaSubj      ONLY the Ispec subject-level
#                                             RSA array (B1), nothing after
#   bash submit_BEBO_pipeline.sh afterRsaSubj everything EXCEPT B1: the whole
#                                             alpha/PF branch (A1-A6) + Ispec
#                                             from per-age on (B2-B7) --
#                                             assumes B1 has finished
#   bash submit_BEBO_pipeline.sh ebosc        ONLY the eBOSC branch (E1-E3)
#                                             -- assumes alpha_pre2post group
#                                             and Ispec fullSample are done
#   bash submit_BEBO_pipeline.sh perAge       start at the per-age jobs
#                                             (B2) -- assumes the Ispec
#                                             subject array has finished;
#                                             skips A
#   bash submit_BEBO_pipeline.sh fullSample   start at the Ispec 2nd-level
#                                             job (B3) and rerun everything
#                                             after it (B4-B7, eBOSC E2-E3)
#                                             + the alpha figure (A6, e.g.
#                                             for new ROI weights) -- assumes
#                                             per-age, alpha group and the
#                                             subject arrays have finished
#   bash submit_BEBO_pipeline.sh roi          start AFTER the Ispec 2nd-level
#                                             stats: alpha ROI (B4) -> single-
#                                             trial (B5) -> aggregation (B6) ->
#                                             Ispec plots (B7) + alpha figure
#                                             (A6) + eBOSC results + plots
#                                             (E2-E3, Ispec ROI channels; NOT
#                                             the eBOSC subject array E1) --
#                                             assumes fullSample, alpha group
#                                             and E1 have finished
#   bash submit_BEBO_pipeline.sh alpha        rerun everything that uses alpha
#                                             (2026-10-02: alpha cluster with
#                                             minnbchan 2, desync + eBOSC in the
#                                             Ispec ROI): alpha subject array
#                                             (A1) -> group (A2); Ispec ROI (B4)
#                                             -> alpha figure (A6, <- A2+B4),
#                                             single-trial (B5) -> aggregation
#                                             (B6) -> Ispec plots (B7); eBOSC
#                                             results + plots (E2-E3). NOT the
#                                             eBOSC subject array (E1) -- assumes
#                                             E1 and Ispec fullSample (B3)
#                                             have finished
#   bash submit_BEBO_pipeline.sh bsim         between-category control pairs
#                                             (BS) -> single-trial (B5) ->
#                                             aggregation (B6) -> Ispec plots
#                                             (B7) + alpha x pair type (B8,
#                                             B9) -- assumes fullSample + ROI
#                                             have finished
#   bash submit_BEBO_pipeline.sh update       everything changed 2026-10-02..05
#                                             without recomputing any cluster
#                                             test: BS + Ispec ROI (B4) ->
#                                             single-trial (B5) -> B6 -> B7, B8,
#                                             B9; alpha figure (A6); eBOSC
#                                             E1 -> E2 -> E3 -- assumes Ispec
#                                             fullSample + alpha group done
#   bash submit_BEBO_pipeline.sh eboscResults only eBOSC results + figure (E2,
#                                             E3) -- assumes the eBOSC subject
#                                             array and the Ispec ROI are done
#   bash submit_BEBO_pipeline.sh stats        start at ALL cluster tests
#                                             (alpha group A2, B3)
#                                             and rerun everything after
#                                             them -- assumes the subject
#                                             arrays (A1, B1) and per-age
#                                             jobs have finished
# Everything after the start point follows automatically.

cd "$(dirname "$0")" || exit 1 # run from RSA_pipeline/ no matter where it is called from (relative paths below)

FREQLO=30
FREQHI=100
BGMODE=legacy # Ispec broadband: BOSC fit on the analysis band itself
SUBJLIST="subj_list_all96.txt" # all 96 with usable EEG (2026-10-19); must match ids_ALL.m
NSUBJ=$(wc -l < "$SUBJLIST")
mkdir -p "$HOME/logs"

START=${1:-all}
case "$START" in
  all|rsaSubj|afterRsaSubj|ebosc|perAge|fullSample|roi|alpha|stats|bsim|update|eboscResults) ;;
  *) echo "Unknown start point '$START' -- use: all, rsaSubj, afterRsaSubj, ebosc, perAge, fullSample, roi, alpha, stats, bsim, update or eboscResults"; exit 1 ;;
esac

# Dependency line for a batch file: dep <type> <jobid> <jobid> ...
# Empty job IDs (steps skipped because of the start point) are left out;
# with no job IDs left, no dependency line is written at all.
dep() {
  local type=$1; shift
  local ids
  ids=$(echo "$@" | xargs | tr ' ' ':')
  [ -n "$ids" ] && echo "#SBATCH --dependency=${type}:${ids}"
}

# Submit a batch file and stop everything if sbatch refuses it (e.g. a
# bad dependency) -- otherwise later jobs would get an empty dependency.
submit() {
  local jid
  jid=$(sbatch --parsable "$1") || { echo "sbatch FAILED for $1 -- aborting (already-submitted jobs keep running; scancel them if needed)" >&2; exit 1; }
  echo "$jid"
}

HDR="#!/bin/bash
#SBATCH --mail-type NONE"

echo "=== Submitting BEBO pipeline ($NSUBJ subjects, ${FREQLO}-${FREQHI}Hz, start: $START) ==="

A1=""; A2=""; B1=""; BS=""; B2=(); B3=""; E1=""


# =====================================================================
# A. Alpha_pre2post + PF
# =====================================================================
if [ "$START" = all ] || [ "$START" = stats ] || [ "$START" = afterRsaSubj ] || [ "$START" = fullSample ] || [ "$START" = roi ] || [ "$START" = alpha ]; then

if [ "$START" = all ] || [ "$START" = afterRsaSubj ] || [ "$START" = alpha ]; then
cat > tmp_A1.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%A_%a.out
#SBATCH --cpus-per-task 2
#SBATCH --time 4:0:0
#SBATCH --mem 50gb
#SBATCH --array=1-${NSUBJ}
#SBATCH --job-name AlphaPre2Post_subj
module load matlab/R2017b
SUBJ=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $SUBJLIST)
matlab -r "config_and_run_alpha_pre2post_subject('\$SUBJ'); exit"
EOF
A1=$(submit tmp_A1.batch) || exit 1
echo "  A1 alpha_pre2post subject array: $A1"
fi

if [ "$START" != fullSample ] && [ "$START" != roi ]; then
cat > tmp_A2.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 2
#SBATCH --time 12:0:0
#SBATCH --mem 40gb
#SBATCH --job-name AlphaPre2Post_group
$(dep afterok $A1)
module load matlab/R2017b
matlab -r "config_and_run_alpha_pre2post_group(); exit"
EOF
A2=$(submit tmp_A2.batch) || exit 1
echo "  A2 alpha_pre2post group: $A2"
fi


fi


if [ "$START" != ebosc ] && [ "$START" != eboscResults ]; then

# =====================================================================
# D. Behaviour
# =====================================================================
if [ "$START" = all ] || [ "$START" = afterRsaSubj ]; then
cat > tmp_D1.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 0:30:0
#SBATCH --mem 8gb
#SBATCH --job-name Behavior
module load matlab/R2017b
matlab -r "config_and_run_behavior(); exit"
EOF
D1=$(submit tmp_D1.batch) || exit 1
echo "  D1 behaviour figure + stats: $D1"
fi

# =====================================================================
# B. Ispec
# =====================================================================
CTYPE=Ispec

if [ "$START" = all ] || [ "$START" = rsaSubj ]; then
cat > tmp_B1.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%A_%a.out
#SBATCH --cpus-per-task 2
#SBATCH --time 22:0:0
#SBATCH --mem 100gb
#SBATCH --array=1-${NSUBJ}
#SBATCH --job-name RSA_IspecBB_subj
module load matlab/R2017b
SUBJ=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $SUBJLIST)
matlab -r "config_and_run_RSA_subject('\$SUBJ','$CTYPE','$FREQLO','$FREQHI','$BGMODE'); exit"
EOF
B1=$(submit tmp_B1.batch) || exit 1
echo "  B1 Ispec subject array: $B1"
fi

# BS (2026-10-05): lag-matched BETWEEN-category control pairs (subject level,
# same TFR files / background removal / baseline method as B1) for the
# alpha x pair-type analysis; saved next to the Ispec files.
if [ "$START" = all ] || [ "$START" = bsim ] || [ "$START" = update ]; then
cat > tmp_BS.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%A_%a.out
#SBATCH --cpus-per-task 2
#SBATCH --time 22:0:0
#SBATCH --mem 100gb
#SBATCH --array=1-${NSUBJ}
#SBATCH --job-name RSA_IspecBB_betweenCat_subj
module load matlab/R2017b
SUBJ=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $SUBJLIST)
matlab -r "config_and_run_RSA_subject('\$SUBJ','Bsim','$FREQLO','$FREQHI','$BGMODE'); exit"
EOF
BS=$(submit tmp_BS.batch) || exit 1
echo "  BS between-category subject array: $BS"
fi
if [ "$START" = rsaSubj ]; then
  rm -f tmp_B1.batch
  echo ""
  echo "Subject-level RSA submitted (nothing after it). Check progress with: squeue -u \$USER"
  exit 0
fi

if [ "$START" = all ] || [ "$START" = perAge ] || [ "$START" = afterRsaSubj ]; then
for GROUP in YA OA; do
  cat > tmp_B2.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 2
#SBATCH --time 22:0:0
#SBATCH --mem 60gb
#SBATCH --job-name RSA_IspecBB_${GROUP}
$(dep afterok $B1)
module load matlab/R2017b
matlab -r "config_and_run_RSA_group_perAge('$GROUP','$CTYPE','$FREQLO','$FREQHI'); exit"
EOF
  JID=$(submit tmp_B2.batch) || exit 1
  B2+=("$JID"); echo "  B2 Ispec per-age ($GROUP): $JID"

done
fi

if [ "$START" != roi ] && [ "$START" != alpha ] && [ "$START" != bsim ] && [ "$START" != update ]; then
cat > tmp_B3.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 22:0:0
#SBATCH --mem 85gb
#SBATCH --job-name RSA_IspecBB_full
$(dep afterok ${B2[@]})
module load matlab/R2017b
matlab -r "config_and_run_RSA_group_fullSample('$CTYPE','$FREQLO','$FREQHI'); exit"
EOF
B3=$(submit tmp_B3.batch) || exit 1
echo "  B3 Ispec fullSample: $B3"
fi

if [ "$START" != bsim ]; then
cat > tmp_B4.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 0:15:0
#SBATCH --mem 10gb
#SBATCH --job-name RSA_IspecBB_alphaROI
$(dep afterok $B3)
module load matlab/R2017b
matlab -r "config_and_run_Ispec_derive_alphaROI('$FREQLO','$FREQHI'); exit"
EOF
B4=$(submit tmp_B4.batch) || exit 1
echo "  B4 Ispec_derive_alphaROI: $B4"
fi

# A6 alpha figure + stats: needs the alpha group test (A2) only (2026-10-05:
# the desync values are taken in the desync cluster itself, no Ispec ROI)
if [ "$START" = all ] || [ "$START" = stats ] || [ "$START" = afterRsaSubj ] || [ "$START" = fullSample ] || [ "$START" = roi ] || [ "$START" = alpha ] || [ "$START" = update ]; then
cat > tmp_A6.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 1:0:0
#SBATCH --mem 20gb
#SBATCH --job-name AlphaPre2Post_plots
$(dep afterok $A2)
module load matlab/R2017b
matlab -r "config_and_run_alpha_pre2post_plots(); exit"
EOF
A6=$(submit tmp_A6.batch) || exit 1
echo "  A6 alpha figure + stats: $A6"
fi

cat > tmp_B5.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%A_%a.out
#SBATCH --cpus-per-task 2
#SBATCH --time 6:0:0
#SBATCH --mem 80gb
#SBATCH --array=1-${NSUBJ}
#SBATCH --job-name RSA_IspecBB_singleTrialMerged
$(dep afterok $B1 $BS $B3 $B4)
module load matlab/R2017b
SUBJ=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $SUBJLIST)
matlab -r "config_and_run_RSA_singleTrial_merged('\$SUBJ','$FREQLO','$FREQHI'); exit"
EOF
B5=$(submit tmp_B5.batch) || exit 1
echo "  B5 single-trial array: $B5"

cat > tmp_B6.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 0:30:0
#SBATCH --mem 20gb
#SBATCH --job-name RSA_IspecBB_singleTrialAggregate
$(dep afterok $B5)
module load matlab/R2017b
matlab -r "aggregate_singleTrial_merged($FREQLO,$FREQHI); exit"
EOF
B6=$(submit tmp_B6.batch) || exit 1
echo "  B6 single-trial aggregation: $B6"

cat > tmp_B7.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 22:0:0
#SBATCH --mem 20gb
#SBATCH --job-name RSA_IspecBB_plots
$(dep afterok $B6 $A1)
module load matlab/R2017b
matlab -r "config_and_run_RSA_plots('$CTYPE','$FREQLO','$FREQHI'); exit"
EOF
B7=$(submit tmp_B7.batch) || exit 1
echo "  B7 Ispec plots: $B7"

# B8/B9 (2026-10-04): alpha x pair type -- does alpha desync go with more
# item-specific representations? Results (computation + stats) <- B6,
# plots <- B8 only.
cat > tmp_B8.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 0:30:0
#SBATCH --mem 10gb
#SBATCH --job-name RSA_IspecBB_alphaPairResults
$(dep afterok $B6)
module load matlab/R2017b
matlab -r "config_and_run_alpha_pair_results(); exit"
EOF
B8=$(submit tmp_B8.batch) || exit 1
echo "  B8 alpha x pair type, results: $B8"

cat > tmp_B9.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 0:30:0
#SBATCH --mem 8gb
#SBATCH --job-name RSA_IspecBB_alphaPairPlots
$(dep afterok $B8)
module load matlab/R2017b
matlab -r "config_and_run_alpha_pair_plots(); exit"
EOF
B9=$(submit tmp_B9.batch) || exit 1
echo "  B9 alpha x pair type, figure: $B9"

fi # end of D + B (skipped for start point 'ebosc')

# =====================================================================
# E. eBOSC alpha rhythm
# =====================================================================
if [ "$START" = all ] || [ "$START" = afterRsaSubj ] || [ "$START" = ebosc ] || [ "$START" = update ]; then
cat > tmp_E1.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%A_%a.out
#SBATCH --cpus-per-task 1
#SBATCH --time 12:0:0
#SBATCH --mem 20gb
#SBATCH --array=1-${NSUBJ}
#SBATCH --job-name eBOSC_subj
module load matlab/R2017b
SUBJ=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $SUBJLIST)
matlab -r "config_and_run_ebosc_subject('\$SUBJ'); exit"
EOF
E1=$(submit tmp_E1.batch) || exit 1
echo "  E1 eBOSC subject array: $E1"
fi

if [ "$START" != rsaSubj ] && [ "$START" != bsim ]; then
cat > tmp_E2.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 2:0:0
#SBATCH --mem 20gb
#SBATCH --job-name eBOSC_results
$(dep afterok $E1 $B4)
module load matlab/R2017b
matlab -r "config_and_run_ebosc_results(); exit"
EOF
E2=$(submit tmp_E2.batch) || exit 1
echo "  E2 eBOSC results + stats: $E2"

cat > tmp_E3.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 1:0:0
#SBATCH --mem 8gb
#SBATCH --job-name eBOSC_plots
$(dep afterok $E2)
module load matlab/R2017b
matlab -r "config_and_run_ebosc_plots(); exit"
EOF
E3=$(submit tmp_E3.batch) || exit 1
echo "  E3 eBOSC plots: $E3"
fi

rm -f tmp_A?.batch tmp_B?.batch tmp_D?.batch tmp_E?.batch
echo ""
echo "All submitted. Check progress with: squeue -u \$USER"

#!/bin/bash
# submit_ispec_from_fullsample.sh
#
# Resume script for the Ispec broadband pipeline: submits from the
# fullSample step ONWARDS only. Use this when the subject-level array job
# AND both per-age (YA, OA) jobs have already completed successfully
# (e.g. after a subject-list fix mid-pipeline) -- resubmitting the whole
# chain from submit_ispec_broadband_pipeline.sh would needlessly redo
# already-completed, expensive steps.
#
# Chain submitted here (PRIMARY/corrected path only -- does NOT include
# the uncorrected/powspctrmOrig comparison branch, which is independent
# and can be resumed separately if still needed):
#   3.  fullSample job -- 2nd-level cluster stat (500 permutations) +
#       3rd-level extraction
#   4.  plots job -- depends on fullSample
#   4b. Ispec_derive_alphaROI -- derives the non-circular alpha
#       extraction ROI from fullSample's 2nd-level cluster stat
#   5.  single-trial merged array job (96 tasks) -- depends on
#       fullSample AND the derived alpha ROI
#   6.  aggregation job -- depends on the single-trial merged array
#
# Run with: bash submit_ispec_from_fullsample.sh   (from inside RSA_pipeline/)

FREQLO=30
FREQHI=100
CTYPE=Ispec
SUBJLIST="subj_list_all96.txt"
NSUBJ=$(wc -l < "$SUBJLIST")
mkdir -p "$HOME/logs"

echo "=== Resuming Ispec_broadband (${FREQLO}-${FREQHI}Hz) from fullSample onwards ==="

# --- 3. fullSample job -- no dependency, per-age jobs already completed ---
cat > tmp_full.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 22:0:0
#SBATCH --mem 85gb
#SBATCH --job-name RSA_IspecBB_full
module load matlab/R2017b
matlab -r "config_and_run_RSA_group_fullSample('$CTYPE','$FREQLO','$FREQHI'); exit"
EOF
FULL_JOBID=$(sbatch --parsable tmp_full.batch)
echo "  fullSample job (500 permutations): $FULL_JOBID"

# --- 4. plots job -- depends on fullSample ---
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

# --- 4b. Ispec_derive_alphaROI -- depends on fullSample ---
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

# --- 5. Single-trial merged array job (96 tasks) -- depends on
#        fullSample AND the derived alpha ROI (subject-level array
#        itself already completed, no need to depend on it here) ---
cat > tmp_stmerged.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%A_%a.out
#SBATCH --cpus-per-task 2
#SBATCH --time 2:0:0
#SBATCH --mem 40gb
#SBATCH --array=1-${NSUBJ}
#SBATCH --job-name RSA_IspecBB_singleTrialMerged
#SBATCH --dependency=afterok:${FULL_JOBID}:${ALPHAROI_JOBID}
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

rm -f tmp_full.batch tmp_plots.batch tmp_alpharoi.batch tmp_stmerged.batch tmp_stagg.batch
echo ""
echo "Resumed from fullSample onwards. Check progress with: squeue -u \$USER"

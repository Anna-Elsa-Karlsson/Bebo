#!/bin/bash
# submit_ispec_after_fullsample.sh
#
# Resume script for the Ispec broadband pipeline: submits everything
# DOWNSTREAM of the fullSample job only. Use this when the fullSample
# job (2nd-level cluster stat + 3rd-level extraction) has already
# completed successfully -- unlike submit_ispec_from_fullsample.sh, it
# does NOT rerun fullSample.
#
# Chain submitted here (PRIMARY/corrected path only):
#   4.  plots job -- no dependency
#   4b. Ispec_derive_alphaROI -- no dependency (rerun so the ROI is
#       guaranteed to match the current fullSample stat)
#   5.  single-trial merged array job -- depends on the alpha ROI
#   6.  aggregation job -- depends on the single-trial merged array
#
# SUBJLIST must match ids_ALL.m (aggregate_singleTrial_merged.m loops
# over ids_ALL('all')) -- currently 96.
#
# Run with: bash submit_ispec_after_fullsample.sh   (from inside RSA_pipeline/)

FREQLO=30
FREQHI=100
CTYPE=Ispec
RUN_PLOTS=0 # rsa_plots.m still loads the old 93-subject Apow alpha data -- leave off until fixed
SUBJLIST="subj_list_all96.txt"
NSUBJ=$(wc -l < "$SUBJLIST")
PREPROC_ROOT="../PreprocData" # = pdat.rawEnc in config_and_run_RSA_singleTrial_merged.m
mkdir -p "$HOME/logs"

# --- Pre-check: step 5 needs each subject's raw encoding data ---
n_missing=0
while read -r SUBJ; do
    [ -z "$SUBJ" ] && continue
    f="$PREPROC_ROOT/$SUBJ/eeg/${SUBJ}_Enc_allcond.mat"
    if [ ! -f "$f" ]; then
        echo "MISSING: $f"
        n_missing=$((n_missing+1))
    fi
done < "$SUBJLIST"
if [ $n_missing -gt 0 ]; then
    echo "$n_missing subject(s) missing preprocessed encoding data -- fix (copy_preproc_data.sh) before submitting. Aborting."
    exit 1
fi

echo "=== Resuming Ispec_broadband (${FREQLO}-${FREQHI}Hz) after fullSample ($NSUBJ subjects) ==="

# --- 4. plots job (optional) ---
if [ "$RUN_PLOTS" -eq 1 ]; then
cat > tmp_plots.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 22:0:0
#SBATCH --mem 20gb
#SBATCH --job-name RSA_IspecBB_plots
module load matlab/R2017b
matlab -r "config_and_run_RSA_plots('$CTYPE','$FREQLO','$FREQHI'); exit"
EOF
PLOTS_JOBID=$(sbatch --parsable tmp_plots.batch)
echo "  plots job: $PLOTS_JOBID"
fi

# --- 4b. Ispec_derive_alphaROI ---
cat > tmp_alpharoi.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 0:15:0
#SBATCH --mem 10gb
#SBATCH --job-name RSA_IspecBB_alphaROI
module load matlab/R2017b
matlab -r "config_and_run_Ispec_derive_alphaROI('$FREQLO','$FREQHI'); exit"
EOF
ALPHAROI_JOBID=$(sbatch --parsable tmp_alpharoi.batch)
echo "  Ispec_derive_alphaROI job: $ALPHAROI_JOBID"

# --- 5. Single-trial merged array job -- depends on the alpha ROI ---
cat > tmp_stmerged.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%A_%a.out
#SBATCH --cpus-per-task 2
#SBATCH --time 2:0:0
#SBATCH --mem 40gb
#SBATCH --array=1-${NSUBJ}
#SBATCH --job-name RSA_IspecBB_singleTrialMerged
#SBATCH --dependency=afterok:${ALPHAROI_JOBID}
module load matlab/R2017b
SUBJ=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $SUBJLIST)
matlab -r "config_and_run_RSA_singleTrial_merged('\$SUBJ','$FREQLO','$FREQHI'); exit"
EOF
ST_JOBID=$(sbatch --parsable tmp_stmerged.batch)
echo "  single-trial merged array job ($NSUBJ tasks): $ST_JOBID"

# --- 6. Aggregation job -- afterany (not afterok) so one failed subject
#        doesn't block aggregation; missing subjects are warned about
#        and skipped by aggregate_singleTrial_merged.m ---
cat > tmp_stagg.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 0:30:0
#SBATCH --mem 20gb
#SBATCH --job-name RSA_IspecBB_singleTrialAggregate
#SBATCH --dependency=afterany:${ST_JOBID}
module load matlab/R2017b
matlab -r "aggregate_singleTrial_merged($FREQLO,$FREQHI); exit"
EOF
AGG_JOBID=$(sbatch --parsable tmp_stagg.batch)
echo "  single-trial aggregation job: $AGG_JOBID"

rm -f tmp_plots.batch tmp_alpharoi.batch tmp_stmerged.batch tmp_stagg.batch
echo ""
echo "Submitted. Check progress with: squeue -u \$USER"

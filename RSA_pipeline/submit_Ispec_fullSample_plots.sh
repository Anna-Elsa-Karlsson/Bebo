#!/bin/bash
# submit_Ispec_fullSample_plots.sh
#
# Runs ONLY:
#   1. Ispec fullSample, baseline-corrected (2nd-level cluster stat +
#      3rd-level cluster means, YA vs OA)
#   2. Ispec plots (Fig3/Fig4 + stats file) -- starts when 1 has succeeded
#
# Assumes the per-age jobs (YA, OA) have finished. The plots also read the
# per-subject alpha (Results/SingleTrialMerged/.../Ispec_30to100Hz_alphaIspec_fullSample.mat)
# from the last single-trial run -- that is NOT recomputed here.
#
# Run with: bash submit_Ispec_fullSample_plots.sh   (from inside RSA_pipeline/)

cd "$(dirname "$0")" || exit 1
mkdir -p "$HOME/logs"

cat > tmp_full.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 22:0:0
#SBATCH --mem 85gb
#SBATCH --job-name RSA_IspecBB_full
module load matlab/R2017b
matlab -r "config_and_run_RSA_group_fullSample('Ispec','30','100'); exit"
EOF
FULL=$(sbatch --parsable tmp_full.batch) || { echo "sbatch FAILED for fullSample"; exit 1; }
echo "  Ispec fullSample (corrected): $FULL"

cat > tmp_plots.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 22:0:0
#SBATCH --mem 20gb
#SBATCH --job-name RSA_IspecBB_plots
#SBATCH --dependency=afterok:${FULL}
module load matlab/R2017b
matlab -r "config_and_run_RSA_plots('Ispec','30','100'); exit"
EOF
PLOTS=$(sbatch --parsable tmp_plots.batch) || { echo "sbatch FAILED for plots"; exit 1; }
echo "  Ispec plots (after $FULL): $PLOTS"

rm -f tmp_full.batch tmp_plots.batch
echo "Submitted. Check progress with: squeue -u \$USER"

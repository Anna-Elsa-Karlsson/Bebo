#!/bin/bash
# submit_Ispec_plots.sh
#
# Runs ONLY the Ispec plot/stats job (rsa_plots.m: Fig3, Fig4, Fig4b +
# stats file). Assumes fullSample, alpha ROI, single-trial aggregation and
# the alpha_pre2post per-subject files already exist.
#
# Run with: bash submit_Ispec_plots.sh   (from inside RSA_pipeline/)

cd "$(dirname "$0")" || exit 1
mkdir -p "$HOME/logs"
cat > tmp_plots.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 4:0:0
#SBATCH --mem 20gb
#SBATCH --job-name RSA_IspecBB_plots
module load matlab/R2017b
matlab -r "config_and_run_RSA_plots('Ispec','30','100'); exit"
EOF
JID=$(sbatch --parsable tmp_plots.batch) && echo "Submitted Ispec plots: $JID (log: ~/logs/slurm-$JID.out)"
rm -f tmp_plots.batch

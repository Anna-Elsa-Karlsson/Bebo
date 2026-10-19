#!/bin/bash
# submit_lag_effects_fullSample.sh
#
# Runs test_lag_effects_onesubj.m as a 93-task array job (one per
# subject), then submits aggregate_lag_effects.m to run once afterward.
# Diagnostic only -- does not touch the production pipeline.
#
# Usage: bash submit_lag_effects_fullSample.sh   (from inside RSA_pipeline/)

SUBJLIST="subj_list_all93.txt"
NSUBJ=$(wc -l < "$SUBJLIST")
mkdir -p "$HOME/logs"

cat > tmp_lagfull.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%A_%a.out
#SBATCH --cpus-per-task 2
#SBATCH --time 8:0:0
#SBATCH --mem 40gb
#SBATCH --array=1-${NSUBJ}
#SBATCH --job-name RSA_lagEffects_full
module load matlab/R2017b
SUBJ=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $SUBJLIST)
matlab -r "test_lag_effects_onesubj('\$SUBJ','30','100','legacy'); exit"
EOF
ARRAY_JOBID=$(sbatch --parsable tmp_lagfull.batch)
echo "Lag effects array job ($NSUBJ tasks): $ARRAY_JOBID"

cat > tmp_lagagg.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 0:30:0
#SBATCH --mem 20gb
#SBATCH --job-name RSA_lagEffects_aggregate
#SBATCH --dependency=afterok:${ARRAY_JOBID}
module load matlab/R2017b
matlab -r "aggregate_lag_effects; exit"
EOF
AGG_JOBID=$(sbatch --parsable tmp_lagagg.batch)
echo "Aggregation job (depends on array): $AGG_JOBID"

rm -f tmp_lagfull.batch tmp_lagagg.batch

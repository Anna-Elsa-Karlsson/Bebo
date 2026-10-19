#!/bin/bash
# submit_lag_effects_test.sh
#
# Runs test_lag_effects_onesubj.m for one subject (default 3109) on the
# grid. Diagnostic only. See test_lag_effects_onesubj.m header.
#
# Usage: bash submit_lag_effects_test.sh [subj]   (from inside RSA_pipeline/)

SUBJ="${1:-3109}"
mkdir -p "$HOME/logs"

cat > tmp_lagtest.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 2
#SBATCH --time 2:0:0
#SBATCH --mem 40gb
#SBATCH --job-name RSA_lagEffectsTest
module load matlab/R2017b
matlab -r "test_lag_effects_onesubj('$SUBJ','30','100','legacy'); exit"
EOF
JOBID=$(sbatch --parsable tmp_lagtest.batch)
echo "Lag effects diagnostic job for subject $SUBJ: $JOBID"
echo "Log: $HOME/logs/slurm-${JOBID}.out"

rm -f tmp_lagtest.batch

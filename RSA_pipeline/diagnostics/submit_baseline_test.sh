#!/bin/bash
# submit_baseline_test.sh
#
# Runs test_baseline_correction_onesubj.m for one subject (default 3109)
# on the grid. Diagnostic only -- does not touch production pipeline
# output. See test_baseline_correction_onesubj.m header for what it
# compares.
#
# Usage: bash submit_baseline_test.sh [subj]   (from inside RSA_pipeline/)

SUBJ="${1:-3109}"
mkdir -p "$HOME/logs"

cat > tmp_baselinetest.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 2
#SBATCH --time 2:0:0
#SBATCH --mem 80gb
#SBATCH --job-name RSA_baselineTest
module load matlab/R2017b
matlab -r "test_baseline_correction_onesubj('$SUBJ','30','100','legacy'); exit"
EOF
JOBID=$(sbatch --parsable tmp_baselinetest.batch)
echo "Baseline diagnostic job for subject $SUBJ: $JOBID"
echo "Log: $HOME/logs/slurm-${JOBID}.out"

rm -f tmp_baselinetest.batch

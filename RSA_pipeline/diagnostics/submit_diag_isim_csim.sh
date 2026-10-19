#!/bin/bash
# submit_diag_isim_csim.sh
#
# Runs diag_isim_csim_per_person.m on the grid: Isim and Csim of every
# participant (mean +/- SE over pairs, number of pairs), summary panels and
# the consistency check single-trial tables vs group-level values.
# Diagnostic only -- does not touch production pipeline output. Needs the
# finished pipeline (single-trial tables + Ispec_derive_alphaROI output).
#
# Usage (from anywhere): bash ~/BEBO2/RSA_pipeline/diagnostics/submit_diag_isim_csim.sh
# Output: ~/BEBO2/Results/Plots/Ispec_30to100Hz/Diag_IsimCsim_perPerson.png/.pdf,
#         Diag_IsimCsim_summary.png/.pdf;
#         ~/BEBO2/Results/Stats/Ispec_30to100Hz/Diag_IsimCsim_perPerson.csv;
#         numbers at the end of the log.

DIAGDIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$HOME/logs"

cat > "$DIAGDIR/tmp_diagIsimCsim.batch" << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 0:30:0
#SBATCH --mem 8gb
#SBATCH --job-name Diag_IsimCsim
module load matlab/R2017b
matlab -nodisplay -r "try, addpath('$DIAGDIR'); diag_isim_csim_per_person; catch e, disp(getReport(e)); exit(1); end; exit(0)"
EOF
JOBID=$(sbatch --parsable "$DIAGDIR/tmp_diagIsimCsim.batch")
echo "Isim/Csim diagnostic job: $JOBID"
echo "Log: $HOME/logs/slurm-${JOBID}.out"

rm -f "$DIAGDIR/tmp_diagIsimCsim.batch"

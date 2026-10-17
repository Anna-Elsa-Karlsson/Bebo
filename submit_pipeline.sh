#!/bin/bash
# submit_pipeline.sh -- submits the whole BEBO rerun on the grid (SLURM).
# Plan: BEBO/analysis_plan_2026-10-17_FINAL.txt. Settings: bebo_config.m.
#
#   S1  gamma RSA per person (array)        s1_gamma_subject.m
#   S2  low-frequency power per person      s2_alpha_subject.m
#   G1  Ispec cluster test + per person     g1_ispec_group.m      <- S1
#   G2  alpha SME cluster test + per person g2_sme_group.m        <- S2
#   G3  single-trial tables + split half    g3_extract.m          <- G1, S2
#   G4  person-level statistics             g4_person_stats.m     <- G2, G3
#   P1  group-level figures (plots only)    p1_figures.m          <- G4
#   B1  behaviour (figure + stats)          b1_behavior.m         (no dependencies)
# Each job starts only if the jobs it needs FINISHED WITHOUT ERROR.
# The trial-level models (R/models_results.R) run locally on the tables.
#
# Run from inside this folder:
#   bash submit_pipeline.sh one      TEST: S1 and S2 for the first person only
#                                    (check trial_log/ before running everything)
#   bash submit_pipeline.sh          everything
#   bash submit_pipeline.sh group    only G1-G4 (S1 and S2 already done)
#   bash submit_pipeline.sh g3       only G3, G4, P1
#   bash submit_pipeline.sh figs     only P1 (figures)
#   bash submit_pipeline.sh behavior only B1 (behaviour)
# Output: ../Results_v2/   Logs: ~/logs/

cd "$(dirname "$0")" || exit 1
REPO=$(pwd)
SUBJLIST="$REPO/subj_list_all93.txt"
NSUBJ=$(grep -c . "$SUBJLIST")   # non-empty lines
START=${1:-all}
mkdir -p "$HOME/logs"

case "$START" in all|one|group|g3|figs|behavior) ;; *) echo "Unknown start point '$START' (all, one, group, g3, figs or behavior)"; exit 1 ;; esac
if [ ! -d "$REPO/../PreprocData" ]; then echo "PreprocData not found next to this folder"; exit 1; fi

# MATLAB call that exits with an error code when the function fails
mcall() { echo "matlab -nodisplay -nosplash -r \"try, addpath('$REPO'); $1; catch e, disp(getReport(e)); exit(1); end; exit(0)\""; }

# submit a job; arguments: name, time, memory, dependency job IDs, command [, array]
submit() {
  local name=$1 time=$2 mem=$3 deps=$4 cmd=$5 array=$6 f out dep jid
  f=$(mktemp "$REPO/tmp_${name}_XXXXXX")
  out="$HOME/logs/${name}_%j.out"; [ -n "$array" ] && out="$HOME/logs/${name}_%A_%a.out"
  {
    echo "#!/bin/bash"
    echo "#SBATCH --job-name $name"
    echo "#SBATCH --output $out"
    echo "#SBATCH --time $time"
    echo "#SBATCH --mem $mem"
    echo "#SBATCH --cpus-per-task 2"
    [ -n "$array" ] && echo "#SBATCH --array=$array"
    dep=$(echo $deps | xargs | tr ' ' ':')
    [ -n "$dep" ] && echo "#SBATCH --dependency=afterok:$dep"
    echo "module load matlab/R2017b"
    [ -n "$array" ] && echo "SUBJ=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" $SUBJLIST)"
    echo "$cmd"
  } > "$f"
  jid=$(sbatch --parsable "$f") || { echo "sbatch FAILED for $name -- stopping" >&2; rm -f "$f"; exit 1; }
  rm -f "$f"
  echo "$jid"
}

S1=""; S2=""; G1=""
if [ "$START" = behavior ]; then
  B1=$(submit B1_behavior 0:30:0 8gb "" "$(mcall 'b1_behavior')") || exit 1
  echo "B1 behaviour: $B1"
  exit 0
fi
if [ "$START" = one ]; then
  S1=$(submit S1_gamma 12:0:0 40gb "" "$(mcall 's1_gamma_subject($SUBJ)')" "1") || exit 1
  S2=$(submit S2_alpha 4:0:0 20gb "" "$(mcall 's2_alpha_subject($SUBJ)')" "1") || exit 1
  echo "Test jobs for the first person: S1 $S1, S2 $S2. Then check ../Results_v2/trial_log/"
  exit 0
fi
if [ "$START" = all ]; then
  S1=$(submit S1_gamma 12:0:0 40gb "" "$(mcall 's1_gamma_subject($SUBJ)')" "1-$NSUBJ") || exit 1
  echo "S1 gamma RSA per person:   $S1"
  S2=$(submit S2_alpha 4:0:0 20gb "" "$(mcall 's2_alpha_subject($SUBJ)')" "1-$NSUBJ") || exit 1
  echo "S2 alpha power per person: $S2"
fi
if [ "$START" = all ] || [ "$START" = group ]; then
  G1=$(submit G1_ispec 22:0:0 60gb "$S1" "$(mcall 'g1_ispec_group')") || exit 1
  echo "G1 Ispec cluster test:     $G1"
  G2=$(submit G2_sme 12:0:0 30gb "$S2" "$(mcall 'g2_sme_group')") || exit 1
  echo "G2 SME cluster test:       $G2"
fi
G3=""; G4=""
if [ "$START" != figs ]; then
  G3=$(submit G3_extract 6:0:0 30gb "$G1 $S2" "$(mcall 'g3_extract')") || exit 1
  echo "G3 single-trial tables:    $G3"
  G4=$(submit G4_stats 1:0:0 8gb "$G2 $G3" "$(mcall 'g4_person_stats')") || exit 1
  echo "G4 person statistics:      $G4"
fi
P1=$(submit P1_figures 1:0:0 16gb "$G4" "$(mcall 'p1_figures')") || exit 1
echo "P1 figures:                $P1"
if [ "$START" != figs ]; then
  B1=$(submit B1_behavior 0:30:0 8gb "" "$(mcall 'b1_behavior')") || exit 1
  echo "B1 behaviour:              $B1"
fi
echo "All submitted. Check with: squeue -u \$USER"

#!/bin/bash
# submit_IspecAll.sh -- Ispec WITHOUT the response filter (2026-10-18):
# all encoding trials of exemplars 1-2, presentations 1-2, whatever the
# later response ('old', 'similar', 'new' or not tested); four experimental
# conditions only. Same settings as the main Ispec analysis otherwise
# (30-100 Hz, BOSC 'legacy', baseline correction, indepsamplesT t-maps,
# cluster test of the full sample).
#
#   C1  subject array   config_and_run_RSA_subject(<ID>,'IspecAll',...)
#   C2  per age group   config_and_run_RSA_group_perAge(YA/OA,'IspecAll',...)   <- C1
#   C3  full sample     config_and_run_RSA_group_fullSample('IspecAll',...)     <- C2
#
# Output ONLY in the IspecAll_30to100Hz folders (Results/Similarity_Data,
# Results/GA, Results/Stats) -- the main Ispec_30to100Hz results are not
# touched. Cluster result + YA vs OA test: end of the C3 log in ~/logs/.
#
# Run from inside RSA_pipeline/:
#   bash submit_IspecAll.sh            everything (C1 -> C2 -> C3)
#   bash submit_IspecAll.sh perAge     C2 -> C3 only (C1 already done)
#   bash submit_IspecAll.sh fullSample C3 only (C1, C2 already done)

cd "$(dirname "$0")" || exit 1

FREQLO=30
FREQHI=100
BGMODE=legacy
CTYPE=IspecAll
SUBJLIST="subj_list_all93.txt" # must match ids_ALL.m
NSUBJ=$(wc -l < "$SUBJLIST")
mkdir -p "$HOME/logs"

START=${1:-all}
case "$START" in
  all|perAge|fullSample) ;;
  *) echo "Unknown start point '$START' -- use: all, perAge or fullSample"; exit 1 ;;
esac

dep() {
  local type=$1; shift
  local ids
  ids=$(echo "$@" | xargs | tr ' ' ':')
  [ -n "$ids" ] && echo "#SBATCH --dependency=${type}:${ids}"
}
submit() {
  local jid
  jid=$(sbatch --parsable "$1") || { echo "sbatch FAILED for $1 -- aborting" >&2; exit 1; }
  echo "$jid"
}
HDR="#!/bin/bash
#SBATCH --mail-type NONE"

echo "=== Submitting IspecAll ($NSUBJ subjects, ${FREQLO}-${FREQHI}Hz, start: $START) ==="
C1=""; C2=()

if [ "$START" = all ]; then
cat > tmp_C1.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%A_%a.out
#SBATCH --cpus-per-task 2
#SBATCH --time 22:0:0
#SBATCH --mem 100gb
#SBATCH --array=1-${NSUBJ}
#SBATCH --job-name RSA_IspecAll_subj
module load matlab/R2017b
SUBJ=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $SUBJLIST)
matlab -r "config_and_run_RSA_subject('\$SUBJ','$CTYPE','$FREQLO','$FREQHI','$BGMODE'); exit"
EOF
C1=$(submit tmp_C1.batch) || exit 1
echo "  C1 IspecAll subject array: $C1"
fi

if [ "$START" = all ] || [ "$START" = perAge ]; then
for GROUP in YA OA; do
cat > tmp_C2.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 2
#SBATCH --time 22:0:0
#SBATCH --mem 60gb
#SBATCH --job-name RSA_IspecAll_${GROUP}
$(dep afterok $C1)
module load matlab/R2017b
matlab -r "config_and_run_RSA_group_perAge('$GROUP','$CTYPE','$FREQLO','$FREQHI'); exit"
EOF
  JID=$(submit tmp_C2.batch) || exit 1
  C2+=("$JID"); echo "  C2 IspecAll per-age ($GROUP): $JID"
done
fi

cat > tmp_C3.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 22:0:0
#SBATCH --mem 85gb
#SBATCH --job-name RSA_IspecAll_full
$(dep afterok ${C2[@]})
module load matlab/R2017b
matlab -r "config_and_run_RSA_group_fullSample('$CTYPE','$FREQLO','$FREQHI'); exit"
EOF
C3=$(submit tmp_C3.batch) || exit 1
echo "  C3 IspecAll fullSample: $C3  -> cluster + YA vs OA at the end of ~/logs/slurm-$C3.out"

rm -f tmp_C1.batch tmp_C2.batch tmp_C3.batch
echo "All submitted. Check with: squeue -u \$USER"

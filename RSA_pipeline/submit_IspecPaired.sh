#!/bin/bash
# submit_IspecPaired.sh -- Ispec with a PAIRED 1st level (2026-10-19):
# per item X, Isim(X) = X presentation 1 vs 2 and an item-specific
# Csim(X) = X presentation 1 vs presentation 2 of the other exemplar of X's
# category; per person a depsamplesT across items (paired by item), then
# the usual group test of the t-maps (full sample, cluster permutation).
# Same decomposition, background removal, baseline correction, trial rule
# for Isim and ROI-free group test as the main Ispec analysis.
#
#   P1  subject array   config_and_run_RSA_subject(<ID>,'IspecPaired',...)
#   P2  per age group   config_and_run_RSA_group_perAge(YA/OA,'IspecPaired',...) <- P1
#   P3  full sample     config_and_run_RSA_group_fullSample('IspecPaired',...)   <- P2
#
# Output ONLY in the IspecPaired_30to100Hz folders (Results/Similarity_Data,
# Results/GA, Results/Stats) -- nothing of the main run is touched.
# Cluster result + YA vs OA test: end of the P3 log in ~/logs/.
#
# Run from inside RSA_pipeline/:
#   bash submit_IspecPaired.sh             start now
#   bash submit_IspecPaired.sh 1234567     wait until job 1234567 (the main
#                                          run's Ispec subject array,
#                                          RSA_IspecBB_subj) has finished, so
#                                          the cached decomposition is reused
#   bash submit_IspecPaired.sh perAge      P2 -> P3 only (P1 already done)
#   bash submit_IspecPaired.sh fullSample  P3 only

cd "$(dirname "$0")" || exit 1

FREQLO=30
FREQHI=100
BGMODE=legacy
CTYPE=IspecPaired
SUBJLIST="subj_list_all96.txt" # must match ids_ALL.m
NSUBJ=$(wc -l < "$SUBJLIST")
mkdir -p "$HOME/logs"

ARG=${1:-all}
WAITFOR=""
case "$ARG" in
  all|perAge|fullSample) START=$ARG ;;
  *[!0-9]*) echo "Unknown argument '$ARG' -- use: (nothing), a job number, perAge or fullSample"; exit 1 ;;
  *) START=all; WAITFOR=$ARG ;;
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

echo "=== Submitting IspecPaired ($NSUBJ subjects, ${FREQLO}-${FREQHI}Hz, start: $START${WAITFOR:+, after job $WAITFOR}) ==="
P1=""; P2=()

if [ "$START" = all ]; then
cat > tmp_P1.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%A_%a.out
#SBATCH --cpus-per-task 2
#SBATCH --time 22:0:0
#SBATCH --mem 100gb
#SBATCH --array=1-${NSUBJ}
#SBATCH --job-name RSA_IspecPaired_subj
$(dep afterany $WAITFOR)
module load matlab/R2017b
SUBJ=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $SUBJLIST)
matlab -r "config_and_run_RSA_subject('\$SUBJ','$CTYPE','$FREQLO','$FREQHI','$BGMODE'); exit"
EOF
P1=$(submit tmp_P1.batch) || exit 1
echo "  P1 IspecPaired subject array: $P1"
fi

if [ "$START" = all ] || [ "$START" = perAge ]; then
for GROUP in YA OA; do
cat > tmp_P2.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 2
#SBATCH --time 22:0:0
#SBATCH --mem 60gb
#SBATCH --job-name RSA_IspecPaired_${GROUP}
$(dep afterok $P1)
module load matlab/R2017b
matlab -r "config_and_run_RSA_group_perAge('$GROUP','$CTYPE','$FREQLO','$FREQHI'); exit"
EOF
  JID=$(submit tmp_P2.batch) || exit 1
  P2+=("$JID"); echo "  P2 IspecPaired per-age ($GROUP): $JID"
done
fi

cat > tmp_P3.batch << EOF
$HDR
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 22:0:0
#SBATCH --mem 85gb
#SBATCH --job-name RSA_IspecPaired_full
$(dep afterok ${P2[@]})
module load matlab/R2017b
matlab -r "config_and_run_RSA_group_fullSample('$CTYPE','$FREQLO','$FREQHI'); exit"
EOF
P3=$(submit tmp_P3.batch) || exit 1
echo "  P3 IspecPaired fullSample: $P3  -> cluster + YA vs OA at the end of ~/logs/slurm-$P3.out"

rm -f tmp_P1.batch tmp_P2.batch tmp_P3.batch
echo "All submitted. Check with: squeue -u \$USER"

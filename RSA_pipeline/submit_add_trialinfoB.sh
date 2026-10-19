#!/bin/bash
# submit_add_trialinfoB.sh
#
# One grid job per subject (array job over subj_list_all93.txt): adds
# trialinfo2_B to the ORIGINAL sim_mat files in
# Results/Similarity_Data/Ispec_30to100Hz/ and re-saves them uncompressed
# (see add_trialinfoB_one_subject.m). One-off conversion step.
#
# Run with: bash submit_add_trialinfoB.sh   (from inside RSA_pipeline/)

cd "$(dirname "$0")" || exit 1
SUBJLIST="subj_list_all93.txt"
NSUBJ=$(wc -l < "$SUBJLIST")
mkdir -p "$HOME/logs"

cat > tmp_addB.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/addB_%A_%a.out
#SBATCH --cpus-per-task 1
#SBATCH --time 4:0:0
#SBATCH --mem 60gb
#SBATCH --array=1-${NSUBJ}
#SBATCH --job-name addTrialinfoB
module load matlab/R2017b
SUBJ=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $SUBJLIST)
matlab -nodisplay -nosplash -r "add_trialinfoB_one_subject('\$SUBJ'); exit"
EOF
JID=$(sbatch --parsable tmp_addB.batch) && echo "Submitted array job $JID ($NSUBJ subjects). Logs: ~/logs/addB_${JID}_<n>.out"
rm -f tmp_addB.batch

#!/bin/bash
# submit_filter_isim_exemplar12.sh
#
# One job per subject: drops exemplar 3/4 pairs from the within-item
# sim_mat files (see filter_isim_exemplar12_one_subject.m). One-off step,
# run AFTER the trialinfo2_B conversion and BEFORE step 2.
#
# Run with: bash submit_filter_isim_exemplar12.sh   (from inside RSA_pipeline/)

cd "$(dirname "$0")" || exit 1
SUBJLIST="subj_list_all93.txt"
NSUBJ=$(wc -l < "$SUBJLIST")
mkdir -p "$HOME/logs"

cat > tmp_filt.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/filtIsim_%A_%a.out
#SBATCH --cpus-per-task 1
#SBATCH --time 4:0:0
#SBATCH --mem 60gb
#SBATCH --array=1-${NSUBJ}
#SBATCH --job-name filterIsim12
module load matlab/R2017b
SUBJ=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $SUBJLIST)
matlab -nodisplay -nosplash -r "filter_isim_exemplar12_one_subject('\$SUBJ'); exit"
EOF
JID=$(sbatch --parsable tmp_filt.batch) && echo "Submitted array job $JID ($NSUBJ subjects). Logs: ~/logs/filtIsim_${JID}_<n>.out"
rm -f tmp_filt.batch

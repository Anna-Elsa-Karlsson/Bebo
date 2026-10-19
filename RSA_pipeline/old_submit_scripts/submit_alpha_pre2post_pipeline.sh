#!/bin/bash
# submit_alpha_pre2post_pipeline.sh
#
# Grid submission for the alpha activation-vs-zero (pre-to-post) analysis
# -- the non-circular alpha frequency-range/ROI-existence test (replaces
# the old SME-based approach, see session notes). Rerun from scratch here
# specifically to pick up the 3 newly-added subjects (3213, 3111, 4210)
# now included in ids_ALL.m / subj_list_all96.txt.
#
# Chain:
#   1. subject-level array job (96 tasks) -- per-subject wavelet TFR,
#      baseline correction, response-collapsed trial average
#      (alpha_pre2post_subject_level.m)
#   2. group-level job -- grand averages (all/YA/OA) + depsamplesT-vs-zero
#      cluster stat (0.4-1.7s, 6-30Hz, clusteralpha=0.01, minnbchan=4,
#      500 randomizations -- the settings that produced the settled
#      8-15Hz effect) (alpha_pre2post_group_level.m)
#   3. PF_derive_roi job -- derives PF's channel ROI/window/freq band
#      directly from stat_alpha.mat (replaces the old SME-derived,
#      circular ROI) (PF_derive_roi.m)
#   4. PF subject-level array job (96 tasks) -- per-subject raw spectral
#      estimation using the derived ROI (PF_subject_level.m)
#   5. PF group-level job -- background correction, peak detection
#      (Corcoran et al. 2018), group stats, plots (PF_group_level.m)
#
# Output: BEBO2/Alpha_pre2post/ (sibling of RSA_pipeline/) -- per-subject
# <id>_alphaPow.mat files, plus GApow.mat, GApow_YA.mat, GApow_OA.mat,
# stat_alpha.mat. PF output: BEBO2/Results_PF/ (PF_roi_config.mat,
# PF_raw/<id>_PFraw.mat, BEBO_PF_estimates.csv/.mat, BEBO_PF_stats.txt)
# and BEBO2/Plots_PF/ (FigPF_Raincloud, FigPF_GroupSpectra).
#
# Time/mem budgets are placeholder estimates (this decomposition is
# lighter than the gamma DPSS one -- 15 freq bins, wavelet not
# multitaper) -- adjust after the first run if needed.
#
# Run with: bash submit_alpha_pre2post_pipeline.sh   (from inside RSA_pipeline/)

SUBJLIST="subj_list_all96.txt"
NSUBJ=$(wc -l < "$SUBJLIST")
mkdir -p "$HOME/logs"

echo "=== Submitting Alpha_pre2post pipeline ($NSUBJ subjects) ==="

# --- 1. Subject-level array job ---
cat > tmp_alpha_subj.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%A_%a.out
#SBATCH --cpus-per-task 2
#SBATCH --time 2:0:0
#SBATCH --mem 30gb
#SBATCH --array=1-${NSUBJ}
#SBATCH --job-name AlphaPre2Post_subj
module load matlab/R2017b
SUBJ=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $SUBJLIST)
matlab -r "config_and_run_alpha_pre2post_subject('\$SUBJ'); exit"
EOF
SUBJ_JOBID=$(sbatch --parsable tmp_alpha_subj.batch)
echo "  subject-level array job ($NSUBJ tasks): $SUBJ_JOBID"

# --- 2. Group-level job -- depends on the full array finishing ---
cat > tmp_alpha_group.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 2
#SBATCH --time 6:0:0
#SBATCH --mem 40gb
#SBATCH --job-name AlphaPre2Post_group
#SBATCH --dependency=afterok:${SUBJ_JOBID}
module load matlab/R2017b
matlab -r "config_and_run_alpha_pre2post_group(); exit"
EOF
GROUP_JOBID=$(sbatch --parsable tmp_alpha_group.batch)
echo "  group-level job: $GROUP_JOBID"

# --- 3. PF_derive_roi -- depends on the alpha_pre2post group-level job ---
cat > tmp_pf_roi.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 0:15:0
#SBATCH --mem 10gb
#SBATCH --job-name PF_deriveROI
#SBATCH --dependency=afterok:${GROUP_JOBID}
module load matlab/R2017b
matlab -r "config_and_run_PF_derive_roi(); exit"
EOF
PFROI_JOBID=$(sbatch --parsable tmp_pf_roi.batch)
echo "  PF_derive_roi job: $PFROI_JOBID"

# --- 4. PF subject-level array job -- depends on PF_derive_roi ---
cat > tmp_pf_subj.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%A_%a.out
#SBATCH --cpus-per-task 2
#SBATCH --time 1:0:0
#SBATCH --mem 20gb
#SBATCH --array=1-${NSUBJ}
#SBATCH --job-name PF_subj
#SBATCH --dependency=afterok:${PFROI_JOBID}
module load matlab/R2017b
SUBJ=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $SUBJLIST)
matlab -r "config_and_run_PF_subject('\$SUBJ'); exit"
EOF
PFSUBJ_JOBID=$(sbatch --parsable tmp_pf_subj.batch)
echo "  PF subject-level array job ($NSUBJ tasks): $PFSUBJ_JOBID"

# --- 5. PF group-level job -- depends on the PF subject-level array ---
cat > tmp_pf_group.batch << EOF
#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --output $HOME/logs/slurm-%j.out
#SBATCH --cpus-per-task 1
#SBATCH --time 1:0:0
#SBATCH --mem 20gb
#SBATCH --job-name PF_group
#SBATCH --dependency=afterok:${PFSUBJ_JOBID}
module load matlab/R2017b
matlab -r "config_and_run_PF_group(); exit"
EOF
PFGROUP_JOBID=$(sbatch --parsable tmp_pf_group.batch)
echo "  PF group-level job: $PFGROUP_JOBID"

rm -f tmp_alpha_subj.batch tmp_alpha_group.batch tmp_pf_roi.batch tmp_pf_subj.batch tmp_pf_group.batch
echo ""
echo "Alpha_pre2post + PF pipeline submitted end-to-end. Check progress with: squeue -u \$USER"

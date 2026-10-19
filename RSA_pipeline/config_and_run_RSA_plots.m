function config_and_run_RSA_plots(contrastType, freqLoStr, freqHiStr)
% Grid entry point for the plotting step (rsa_plots.m).
% Called as: matlab -r "config_and_run_RSA_plots('Ispec','30','100'); exit"
% Run once per variant, after config_and_run_RSA_group_fullSample
% completes for that variant.
%
% SELF-CONTAINED: see config_and_run_RSA_subject.m header for the
% portable-path convention used throughout this folder.
%
% Needs, as siblings of this RSA_pipeline/ folder (confirmed present on
% the grid at /home/mpib/karlsson/BEBO2, 2026-08-25):
%   Performance_All.txt (matched by ID to ids_ALL; extra rows ignored),
%   Elec/layout_EEG_cap_20161018.mat, and the aggregated per-subject
%   alpha from aggregate_singleTrial_merged.m (Results/SingleTrialMerged/)

pipelineRoot = fileparts(mfilename('fullpath'));
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/'); % Confirmed 2026-08-25
ft_defaults

pdat.ga   = fullfile(pipelineRoot, '..', 'Results','GA');
pdat.stat = fullfile(pipelineRoot, '..', 'Results','Stats');
pdat.fig  = fullfile(pipelineRoot, '..', 'Results','Plots');
pdat.stm  = fullfile(pipelineRoot, '..', 'Results','SingleTrialMerged'); % per-subject alphaIspec (aggregate_singleTrial_merged.m) + Ispec_alphaROI_config.mat
pdat.alpha = fullfile(pipelineRoot, '..', 'Alpha_pre2post');            % per-subject baseline-corrected alpha (<ID>_alphaPow.mat)

baseDir = fullfile(pipelineRoot, '..'); % expects ../Performance_All.txt (behavioural data only)
cfg_layout = fullfile(pipelineRoot, '..', 'Elec', 'layout_EEG_cap_20161018.mat');

rsa_plots(contrastType, str2double(freqLoStr), str2double(freqHiStr), pdat, baseDir, cfg_layout);

end

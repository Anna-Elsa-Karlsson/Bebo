function config_and_run_alpha_pre2post_plots()
% Grid entry point for alpha_pre2post_plots.m. Run after the alpha
% group-level job.
% Called as: matlab -r "config_and_run_alpha_pre2post_plots(); exit"
%
% SELF-CONTAINED: see config_and_run_RSA_subject.m header for the
% portable-path convention used throughout this folder.

pipelineRoot = fileparts(mfilename('fullpath'));
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/');
ft_defaults

pdat.alpha = fullfile(pipelineRoot, '..', 'Alpha_pre2post');
pdat.fig   = fullfile(pipelineRoot, '..', 'Plots_Alpha');
pdat.perfFile = fullfile(pipelineRoot, '..', 'Performance_All.txt'); % IM/LDI, matched to ids_ALL by ID
cfg_layout = fullfile(pipelineRoot, '..', 'Elec', 'layout_EEG_cap_20161018.mat');

alpha_pre2post_plots(pdat, cfg_layout);

end

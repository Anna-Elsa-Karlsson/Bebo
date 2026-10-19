function config_and_run_alpha_pre2post_group()
% Grid entry point for the group-level step (alpha_pre2post_group_level.m).
% Called as: matlab -r "config_and_run_alpha_pre2post_group(); exit"
% Run once, after the full subject-level array job completes.
%
% SELF-CONTAINED: see config_and_run_RSA_subject.m header for the
% portable-path convention used throughout this folder.

pipelineRoot = fileparts(mfilename('fullpath'));
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/');
ft_defaults

pdat.out  = fullfile(pipelineRoot, '..', 'Alpha_pre2post');
pdat.elec = fullfile(pipelineRoot, '..', 'Elec'); % elec.mat lives here (confirmed 2026-09-XX, NOT under TFRdata/ like the RSA pipeline's copy)

alpha_pre2post_group_level(pdat);

end

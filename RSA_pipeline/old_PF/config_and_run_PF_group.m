function config_and_run_PF_group()
% Grid entry point for the PF group-level step (PF_group_level.m).
% Called as: matlab -r "config_and_run_PF_group(); exit"
% Run once, after the full PF subject-level array job completes.

pipelineRoot = fileparts(mfilename('fullpath'));
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/');
ft_defaults

pdat.out     = fullfile(pipelineRoot, '..', 'Results_PF', 'PF_raw');
pdat.roiDir  = fullfile(pipelineRoot, '..', 'Results_PF');
pdat.plotDir = fullfile(pipelineRoot, '..', 'Plots_PF');

PF_group_level(pdat);

end

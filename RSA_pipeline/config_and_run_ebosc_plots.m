function config_and_run_ebosc_plots()
% Grid entry point for ebosc_plots.m (plotting only).
% Called as: matlab -r "config_and_run_ebosc_plots(); exit"
pipelineRoot = fileparts(mfilename('fullpath'));
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/'); % ft_topoplotER (PF topographies)
ft_defaults
pdat.res = fullfile(pipelineRoot, '..', 'Results_eBOSC');
pdat.fig = fullfile(pipelineRoot, '..', 'Plots_eBOSC');
pdat.layout = fullfile(pipelineRoot, '..', 'Elec', 'layout_EEG_cap_20161018.mat');
ebosc_plots(pdat);
end

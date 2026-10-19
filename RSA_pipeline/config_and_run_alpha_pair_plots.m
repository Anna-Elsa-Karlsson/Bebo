function config_and_run_alpha_pair_plots()
% Grid entry point for alpha_pair_plots.m (alpha x pair type: does alpha
% desync go with more item-specific representations?). Called as:
%   matlab -r "config_and_run_alpha_pair_plots(); exit"

pipelineRoot = fileparts(mfilename('fullpath'));
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/');
ft_defaults

pdat.variant = 'Ispec_30to100Hz';
pdat.stm = fullfile(pipelineRoot, '..', 'Results', 'SingleTrialMerged', pdat.variant); % T_pair (aggregate_singleTrial_merged.m)
pdat.out = fullfile(pipelineRoot, '..', 'Results', 'Stats', pdat.variant);            % alpha_pair_results.mat + alpha_pair_stats.txt
pdat.fig = fullfile(pipelineRoot, '..', 'Results', 'Plots', pdat.variant);            % Fig_AlphaPair

alpha_pair_plots(pdat);

end

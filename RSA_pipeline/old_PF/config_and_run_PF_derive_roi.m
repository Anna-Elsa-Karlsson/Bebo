function config_and_run_PF_derive_roi()
% Grid entry point for PF_derive_roi.m. Run once, after
% alpha_pre2post_group_level.m has produced stat_alpha.mat.
% Called as: matlab -r "config_and_run_PF_derive_roi(); exit"

pipelineRoot = fileparts(mfilename('fullpath'));
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/');
ft_defaults

pdat.alphaDir = fullfile(pipelineRoot, '..', 'Alpha_pre2post');
pdat.out      = fullfile(pipelineRoot, '..', 'Results_PF');

PF_derive_roi(pdat);

end

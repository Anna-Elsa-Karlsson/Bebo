function config_and_run_ebosc_subject(subjectIdStr)
% Grid entry point for ebosc_subject.m (one subject).
% Called as: matlab -r "config_and_run_ebosc_subject('3109'); exit"
pipelineRoot = fileparts(mfilename('fullpath'));
addpath(pipelineRoot);                                             % alpha_band.m
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/');
ft_defaults
pdat.rawEnc = fullfile(pipelineRoot, '..', 'PreprocData');
pdat.fun    = fullfile(pipelineRoot, 'ebosc');                     % eBOSC functions (from CFC_CMY_bosc)
pdat.out    = fullfile(pipelineRoot, '..', 'Results_eBOSC', 'subjects');
ebosc_subject(subjectIdStr, pdat);
end

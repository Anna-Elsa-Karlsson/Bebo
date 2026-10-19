function config_and_run_PF_subject(subjectIdStr)
% Grid entry point for the PF per-subject step (PF_subject_level.m).
% Called as: matlab -r "config_and_run_PF_subject('3109'); exit"
% Run only after config_and_run_PF_derive_roi has produced PF_roi_config.mat.

pipelineRoot = fileparts(mfilename('fullpath'));
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/');
ft_defaults

pdat.rawEnc = fullfile(pipelineRoot, '..', 'PreprocData');
pdat.roiDir = fullfile(pipelineRoot, '..', 'Results_PF');
pdat.out    = fullfile(pipelineRoot, '..', 'Results_PF', 'PF_raw');

PF_subject_level(subjectIdStr, pdat);

end

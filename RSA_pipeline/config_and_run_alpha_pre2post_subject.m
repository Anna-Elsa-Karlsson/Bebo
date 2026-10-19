function config_and_run_alpha_pre2post_subject(subjectIdStr)
% Grid entry point for the per-subject step (alpha_pre2post_subject_level.m).
% Called as: matlab -r "config_and_run_alpha_pre2post_subject('3109'); exit"
%
% SELF-CONTAINED: see config_and_run_RSA_subject.m header for the
% portable-path convention used throughout this folder.

pipelineRoot = fileparts(mfilename('fullpath'));
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/');
ft_defaults

pdat.rawEnc = fullfile(pipelineRoot, '..', 'PreprocData');
pdat.out    = fullfile(pipelineRoot, '..', 'Alpha_pre2post');

alpha_pre2post_subject_level(subjectIdStr, pdat);

end

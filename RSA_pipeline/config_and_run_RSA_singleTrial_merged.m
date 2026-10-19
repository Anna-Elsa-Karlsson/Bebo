function config_and_run_RSA_singleTrial_merged(subjectIdStr, freqLoStr, freqHiStr)
% Grid entry point for rsa_singleTrial_merged.m -- REPLACES
% config_and_run_RSA_singleTrial_alpha.m + config_and_run_RSA_
% singleTrial_extract.m. Runs once per subject per Ispec variant.
% Called as: matlab -r "config_and_run_RSA_singleTrial_merged('3109','30','100'); exit"
%
% SELF-CONTAINED: see config_and_run_RSA_subject.m header for the
% portable-path convention used throughout this folder.

pipelineRoot = fileparts(mfilename('fullpath'));
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/');
ft_defaults

pdat.rsa = fullfile(pipelineRoot, '..', 'Results', 'Similarity_Data');
pdat.stat = fullfile(pipelineRoot, '..', 'Results', 'Stats');
% Non-circular alpha ROI (2026-09-XX): channels + time window derived
% from the Ispec cluster (Ispec_derive_alphaROI.m), NOT the old SME-based
% Masks_IvC.mat (Apow/) -- see rsa_singleTrial_merged.m header.
pdat.ispecAlphaROI = fullfile(pipelineRoot, '..', 'Results', 'SingleTrialMerged');
pdat.rawEnc = fullfile(pipelineRoot, '..', 'PreprocData');
pdat.out = fullfile(pipelineRoot, '..', 'Results', 'SingleTrialMerged');

rsa_singleTrial_merged(subjectIdStr, str2double(freqLoStr), str2double(freqHiStr), pdat);

end

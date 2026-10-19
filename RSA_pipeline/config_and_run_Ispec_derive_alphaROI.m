function config_and_run_Ispec_derive_alphaROI(freqLoStr, freqHiStr)
% Grid entry point for Ispec_derive_alphaROI.m. Run once, after the Ispec
% fullSample job (config_and_run_RSA_group_fullSample) has produced the
% 2nd-level cluster stat.
% Called as: matlab -r "config_and_run_Ispec_derive_alphaROI('30','100'); exit"

pipelineRoot = fileparts(mfilename('fullpath'));
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/');
ft_defaults

pdat.stat = fullfile(pipelineRoot, '..', 'Results', 'Stats');
pdat.ga   = fullfile(pipelineRoot, '..', 'Results', 'GA');     % grand averages for the per-subject ROI values
pdat.out  = fullfile(pipelineRoot, '..', 'Results', 'SingleTrialMerged'); % same dir rsa_singleTrial_merged.m reads from

Ispec_derive_alphaROI(str2double(freqLoStr), str2double(freqHiStr), pdat);

end

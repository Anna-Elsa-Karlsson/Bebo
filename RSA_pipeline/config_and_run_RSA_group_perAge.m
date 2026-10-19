function config_and_run_RSA_group_perAge(group, contrastType, freqLoStr, freqHiStr, simField)
% Grid entry point for the per-age-group step (rsa_group_level_perAge.m).
% Called as: matlab -r "config_and_run_RSA_group_perAge('YA','Ispec','30','100'); exit"
% Run once per {group} x {variant} -- e.g. 2 age groups x 7 variants = 14 jobs.
% Must complete (both YA and OA) before config_and_run_RSA_group_fullSample
% runs for the same variant.
%
% simField: (optional, default 'powspctrm') 'powspctrm' (baseline-
% corrected, primary) or 'powspctrmOrig' (uncorrected) -- 2026-09-XX
% addition. Pass 'powspctrmOrig' as a 5th arg to run the uncorrected
% comparison; must be run separately (a second call) alongside the
% default corrected run, not instead of it.
%
% SELF-CONTAINED: see config_and_run_RSA_subject.m header for the
% portable-path convention used throughout this folder.

pipelineRoot = fileparts(mfilename('fullpath'));
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/'); % Confirmed 2026-08-25
ft_defaults

pdat.rsa  = fullfile(pipelineRoot, '..', 'Results','Similarity_Data');
pdat.ga   = fullfile(pipelineRoot, '..', 'Results','GA');
pdat.stat = fullfile(pipelineRoot, '..', 'Results','Stats');

if nargin < 5 || isempty(simField), simField = 'powspctrm'; end
rsa_group_level_perAge(group, contrastType, str2double(freqLoStr), str2double(freqHiStr), pdat, simField);

end

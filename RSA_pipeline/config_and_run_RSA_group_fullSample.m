function config_and_run_RSA_group_fullSample(contrastType, freqLoStr, freqHiStr, simField)
% Grid entry point for the fullSample step (rsa_group_level_fullSample.m).
% Called as: matlab -r "config_and_run_RSA_group_fullSample('Ispec','30','100'); exit"
% Run once per variant (7 total) -- ONLY after both YA and OA
% config_and_run_RSA_group_perAge jobs for this same variant have
% completed (explicit SLURM --dependency in submit_rsa_all_variants.sh).
%
% This is the expensive step (montecarlo cluster permutation, 500
% randomizations, on the full 60x213x213 time-time-channel space) -- up
% to ~22h per the original pipeline's SBATCH --time setting.
%
% simField: (optional, default 'powspctrm') must match whatever
% config_and_run_RSA_group_perAge was run with for these subjects --
% 'powspctrm' (corrected) or 'powspctrmOrig' (uncorrected).
%
% SELF-CONTAINED: see config_and_run_RSA_subject.m header for the
% portable-path convention used throughout this folder.

pipelineRoot = fileparts(mfilename('fullpath'));
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/'); % Confirmed 2026-08-25
ft_defaults

pdat.elec = fullfile(pipelineRoot, '..', 'Elec'); % elec.mat -- same file as the alpha_pre2post group step
pdat.ga   = fullfile(pipelineRoot, '..', 'Results','GA');
pdat.stat = fullfile(pipelineRoot, '..', 'Results','Stats');

if nargin < 4 || isempty(simField), simField = 'powspctrm'; end
rsa_group_level_fullSample(contrastType, str2double(freqLoStr), str2double(freqHiStr), pdat, simField);

end

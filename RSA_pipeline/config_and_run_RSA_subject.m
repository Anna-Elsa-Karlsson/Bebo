function config_and_run_RSA_subject(subjectIdStr, contrastType, freqLoStr, freqHiStr, bgFitMode)
% Grid entry point for the per-subject step (rsa_subject_level.m).
% Called as: matlab -r "config_and_run_RSA_subject('3109','Ispec','30','100','wide'); exit"
% bgFitMode: 'wide' (default, used for 6 of 7 variants) or 'legacy'
% (Ispec broadband only -- see rsa_subject_level.m header for why).
%
% SELF-CONTAINED: this whole RSA_pipeline/ folder can be copied anywhere
% on the grid and run as-is -- no path edits needed. All paths below are
% built relative to this file's own location (assumed to sit as a direct
% subfolder of the BEBO2 root, alongside TFRdata/, TFRdat_BL/), not
% hardcoded absolute paths. FieldTrip itself is the one exception (a
% large external toolbox, not something to vendor into this folder) --
% verify that path is still correct for your grid.

pipelineRoot = fileparts(mfilename('fullpath'));
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/'); % Confirmed 2026-08-25
ft_defaults

% Raw TFR data -- lives one level up (siblings of this folder), NOT
% copied in (large, shared, already-existing grid resource). Only used
% for theta/low-freq variants (freqHi<=20) -- see rsa_subject_level.m's
% load_bosc_correct_restrict_legacyFile.
pdat.tfr   = fullfile(pipelineRoot, '..', 'TFRdata', 'Time-frequencyEEGdata');
pdat.tfrBL = fullfile(pipelineRoot, '..', 'TFRdat_BL');

% Preprocessed data + raw-decomposition cache -- used for gamma variants
% (freqHi>20), which now decompose directly from Enc_allcond instead of
% reading a pre-made file (2026-09-04, see rsa_subject_level.m's
% load_bosc_correct_restrict_rawGamma). Same .rawEnc convention as
% config_and_run_RSA_singleTrial_merged.m's pdat.rawEnc.
pdat.rawEnc = fullfile(pipelineRoot, '..', 'PreprocData');
pdat.tfrOut = fullfile(pipelineRoot, '..', 'Results', 'GammaTFR_Cache'); % moved out of RSA_pipeline/ (2026-09-05) -- deleting/recopying the code folder must never wipe the cache/results again

% Output -- lives as a SIBLING of this folder (BEBO2/Results/, not
% RSA_pipeline/Results/), deliberately decoupled from the code (2026-09-05)
% so deleting/recopying RSA_pipeline/ to deploy a code update can never
% wipe the cache or any computed results again.
pdat.rsa = fullfile(pipelineRoot, '..', 'Results', 'Similarity_Data');

if nargin < 5, bgFitMode = 'wide'; end
rsa_subject_level(subjectIdStr, contrastType, str2double(freqLoStr), str2double(freqHiStr), pdat, bgFitMode);

end

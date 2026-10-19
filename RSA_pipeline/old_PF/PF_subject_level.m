function PF_subject_level(subj, pdat)
% Per-subject raw spectral estimation for peak alpha frequency (PF).
% Ported from BEBO_PF_analysis.m Section 1, unchanged in method (FFT/
% multitaper, single Hanning taper, trial-averaged then ROI-averaged) --
% only the ROI/window source has changed (loaded from PF_roi_config.mat,
% produced by PF_derive_roi.m from the alpha_pre2post cluster, instead of
% hardcoded SME-derived values).
%
%   subj : subject ID (numeric or string)
%   pdat : struct with fields .rawEnc (PreprocData root), .roiDir
%          (contains PF_roi_config.mat), .out (per-subject output dir)
%
% Output: saves <subj>_PFraw.mat (spec_raw row vector over foiGrid,
% nTrials) into pdat.out.

subj = num2str(subj);

L = load(fullfile(pdat.roiDir, 'PF_roi_config.mat'), 'roiChannels', 'winLat');
roiChannels = L.roiChannels;
winLat = L.winLat;

foiGrid = 2:0.25:20; % broad enough for QC + background flanks regardless of the derived search band
padDur = 4;          % seconds -- zero-pad for a smooth 0.25Hz-sampled spectrum

loadFile = fullfile(pdat.rawEnc, subj, 'eeg', [subj '_Enc_allcond.mat']);
L = load(loadFile, 'data');
data_preprocessed = L.data; clear L

% Trial selection: presentations 1-2, exemplars 1-2, response in
% {old, similar} -- exactly the trials that enter the Ispec analysis,
% same rule as every other analysis (2026-09-26; was presentation 1 only).
trlInfo = data_preprocessed.trialinfo2;
trlIdx  = find(cell2mat(trlInfo(:,9)) <= 2 & cell2mat(trlInfo(:,8)) < 3 & ...
               (strcmp(trlInfo(:,11),'old') | strcmp(trlInfo(:,11),'similar')));
nTrials = length(trlIdx);

cfg_sel = [];
cfg_sel.trials  = trlIdx;
cfg_sel.channel = roiChannels;
cfg_sel.latency = winLat;
data_win = ft_selectdata(cfg_sel, data_preprocessed);
clear data_preprocessed

cfg_freq = [];
cfg_freq.method     = 'mtmfft';
cfg_freq.taper      = 'hanning';
cfg_freq.output     = 'pow';
cfg_freq.keeptrials = 'yes';
cfg_freq.foi        = foiGrid;
cfg_freq.pad        = padDur;
tfr = ft_freqanalysis(cfg_freq, data_win);
clear data_win

pow_trialavg = squeeze(mean(tfr.powspctrm, 1, 'omitnan')); % chan x freq
spec_raw = mean(pow_trialavg, 1, 'omitnan');                % 1 x freq

if ~exist(pdat.out, 'dir'), mkdir(pdat.out); end
save(fullfile(pdat.out, [subj '_PFraw.mat']), 'spec_raw', 'nTrials', 'foiGrid', '-v7.3')
fprintf('Saved PF raw spectrum for %s -> %s\n', subj, pdat.out);

end

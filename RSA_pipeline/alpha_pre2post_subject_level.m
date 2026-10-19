function alpha_pre2post_subject_level(subj, pdat)
% Per-subject step for the alpha activation-vs-zero (pre-to-post) analysis.
% Faithful port of the user's local alpha-pre2post script's per-subject
% loop body (wavelet TFR -> channel/latency select -> log10 -> baseline
% correct -> response-collapsed trial average), parameterized for the
% grid array job instead of looping over all subjects in one process.
%
%   subj : subject ID (numeric or string), e.g. 3109
%   pdat : struct with fields .rawEnc (PreprocData root) and .out
%          (per-subject output dir)
%
% Output: saves <subj>_alphaPow.mat (one response-collapsed,
% baseline-corrected TFR descriptive per subject, 'tmp' -> saved as
% 'pow') into pdat.out.

subj = num2str(subj);

ChanSelection = {'all','-A1','-IOL','-LHEOG','-RHEOG'};

%% Load preprocessed data
load(fullfile(pdat.rawEnc, [subj '_Enc_allcond.mat']), 'data')
data_preprocessed = data; clear data

%% TFR decomposition -- unchanged from the validated local script
cfg = [];
cfg.output     = 'pow';
cfg.channel    = 'all';
cfg.width      = 7;
cfg.keeptrials = 'yes';
cfg.method     = 'wavelet';
cfg.foi        = 2:2:30;
cfg.toi        = -2:0.04:3.5;
cfg.pad        = ceil(max(cellfun(@numel, data_preprocessed.time)) / data_preprocessed.fsample);

TFR_data_WV = ft_freqanalysis(cfg, data_preprocessed);

%% Select channels + latency
cfg = [];
cfg.channel = ChanSelection;
cfg.latency = [-0.6 1.7];
TFR_data_WV = ft_selectdata(cfg, TFR_data_WV);
TFR_data_WV.trialinfo2 = data_preprocessed.trialinfo2;

%% log10 transform
cfg = []; cfg.parameter = 'powspctrm'; cfg.operation = 'log10';
TFR_data_WV = ft_math(cfg, TFR_data_WV);
TFR_data_WV.trialinfo2 = data_preprocessed.trialinfo2;

%% Baseline correct (absolute, -0.6 to -0.2s)
cfg = [];
cfg.baseline = [-0.6 -0.2];
cfg.baselinetype = 'absolute';
cfg.keeptrials = 'yes';
TFR_data_WV = ft_freqbaseline(cfg, TFR_data_WV);
TFR_data_WV.trialinfo2 = data_preprocessed.trialinfo2;
clear data_preprocessed

%% Select trials: presentations 1-2, exemplars 1-2, valid (old/similar)
%  response -- exactly the trials that enter the Ispec analysis (Isim:
%  presentations 1+2 of each item; Csim: presentation 1 of exemplars 1+2),
%  same rule as every other analysis (2026-09-26; was presentation 1
%  only). Response-collapsed (both memory outcomes combined), so
%  non-circular w.r.t. the old-vs-similar contrast used elsewhere.
trlInfo = TFR_data_WV.trialinfo2;
trlIdx_I = find(cell2mat(trlInfo(:,9)) <= 2 & cell2mat(trlInfo(:,8)) < 3 & strcmp(trlInfo(:,11),'old'));
trlIdx_C = find(cell2mat(trlInfo(:,9)) <= 2 & cell2mat(trlInfo(:,8)) < 3 & strcmp(trlInfo(:,11),'similar'));

cfg = [];
cfg.keeptrials = 'no';
cfg.trials = [trlIdx_I; trlIdx_C];
pow = ft_freqdescriptives(cfg, TFR_data_WV);
cfg.trials = trlIdx_I;
pow_old = ft_freqdescriptives(cfg, TFR_data_WV);
cfg.trials = trlIdx_C;
pow_sim = ft_freqdescriptives(cfg, TFR_data_WV);
nTrl_old = numel(trlIdx_I);
nTrl_sim = numel(trlIdx_C);
clear TFR_data_WV

%% Save
if ~exist(pdat.out, 'dir'), mkdir(pdat.out); end
save(fullfile(pdat.out, [subj '_alphaPow.mat']), 'pow', '-v7.3')
save(fullfile(pdat.out, [subj '_alphaPow_byResp.mat']), 'pow_old', 'pow_sim', 'nTrl_old', 'nTrl_sim', '-v7.3')
fprintf('Saved alpha pre2post subject-level pow for %s -> %s\n', subj, pdat.out);

end

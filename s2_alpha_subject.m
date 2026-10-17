function s2_alpha_subject(subj)
% S2_ALPHA_SUBJECT  Per-person low-frequency power (analysis plan 1b and 2).
% Computes everything from the preprocessed data.
%
%   1. load PreprocData/<ID>_Enc_allcond.mat, keep the analysed trials
%      (bebo_trials.m, the same trials as the gamma RSA) and the 60 EEG
%      channels
%   2. TFR: Morlet wavelets, 7 cycles, 2:2:30 Hz, 40 ms steps,
%      -0.28 to 1.68 s; log10 per trial; NO baseline correction
%   3. condition averages for the SME: presentation-1 trials of items later
%      called 'old' and of items later called 'similar' (one trial per item)
%
% Output (C.dir.alphaSubj/<ID>_alpha.mat):
%   pow        single-trial log10 power [trial x chan x freq x time], single,
%              rows = rows of T (used for the trial-level power, plan 2)
%   powOld, powSim   mean log10 power, later 'old' / later 'similar',
%                    presentation 1 [chan x freq x time]
%   nOld, nSim       number of items (presentation-1 trials) in each
%   T, label, freq, time
% Trial log (C.dir.trialLog/<ID>_alpha_log.csv).

if isnumeric(subj), subj = num2str(subj); end
C = bebo_config();
bebo_init_fieldtrip(C);
for d = {C.dir.alphaSubj, C.dir.trialLog}, if ~exist(d{1}, 'dir'), mkdir(d{1}); end, end

%% 1. Data and trial rule
fprintf('[%s] loading preprocessed data\n', subj);
L = load(fullfile(C.preproc, [subj '_Enc_allcond.mat']), 'data');
data = L.data; clear L
T = bebo_trials(data.trialinfo2, C);
[selRows, ~, dataRow] = unique(T.row);
T.dataRow = dataRow;
cfg = []; cfg.trials = selRows; cfg.channel = C.channel;
data = ft_selectdata(cfg, data);

%% 2. TFR
cfg = [];
cfg.method     = 'wavelet';
cfg.output     = 'pow';
cfg.width      = C.alpha.width;
cfg.foi        = C.alpha.foi;
cfg.toi        = C.alpha.toi;
cfg.keeptrials = 'yes';
cfg.pad        = C.alpha.pad;     % 6 s: frequency grid exactly 2:2:30 Hz
fprintf('[%s] TFR (Morlet %d cycles, %d-%d Hz)\n', subj, cfg.width, cfg.foi(1), cfg.foi(end));
tfr = ft_freqanalysis(cfg, data);
clear data
P = log10(tfr.powspctrm);
label = tfr.label; freq = tfr.freq; time = tfr.time;
clear tfr
if any(~isfinite(P(:)))
    error('s2_alpha_subject:nan', '[%s] non-finite power (check that the wavelets fit into the epochs)', subj);
end

% rows in T order
pow = single(P(T.dataRow, :, :, :));
clear P

%% 3. Condition averages (SME): PRESENTATION 1 only (one trial per item;
%  decided 2026-10-17). Single-trial power is kept for presentations 1 and 2
%  (the pair model needs the power of both trials of an Isim pair).
isOld = T.Acc == 1 & T.pres == 1;
isSim = T.Acc == 0 & T.pres == 1;
nOld = sum(isOld); nSim = sum(isSim);
if nOld < 2 || nSim < 2
    error('s2_alpha_subject:cond', '[%s] too few trials (old %d, similar %d)', subj, nOld, nSim);
end
powOld = squeeze(mean(double(pow(isOld,:,:,:)), 1));
powSim = squeeze(mean(double(pow(isSim,:,:,:)), 1));

%% Save
save(fullfile(C.dir.alphaSubj, [subj '_alpha.mat']), 'pow', 'powOld', 'powSim', 'nOld', 'nSim', ...
    'T', 'label', 'freq', 'time', '-v7.3');
tlog = table(string(subj), height(T), nOld, nSim, 'VariableNames', {'ID','nTrials','nItemsOld','nItemsSimilar'});
writetable(tlog, fullfile(C.dir.trialLog, [subj '_alpha_log.csv']));
fprintf('[%s] done: %d trials; SME items (presentation 1): later old %d, later similar %d\n', subj, height(T), nOld, nSim);
end

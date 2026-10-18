function s1_gamma_subject(subj)
% S1_GAMMA_SUBJECT  Per-person gamma RSA with the original pre-stimulus
% baseline correction (analysis plan 1a; decided 2026-10-18). Computes
% everything from the preprocessed data.
%
%   1. load PreprocData/<ID>_Enc_allcond.mat, keep the analysed trials
%      (bebo_trials.m) and the 60 EEG channels
%   2. line noise: DFT filter at 50 and 100 Hz
%   3. TFR: DPSS mtmconvol, 35:2:100 Hz, +/-5 Hz, 0.4 s, 20 ms, 0-1.7 s;
%      pre-stimulus baseline TFR with the same settings, window centres
%      -1.0 to -0.7 s (data -1.2 to -0.5 s)
%   4. aperiodic background removed (BOSC fit on 35-100 Hz) separately for
%      activation and baseline, log10; baseline averaged over time
%   5. similarity of every Isim and Csim pair: Pearson r between the two
%      trials' spectra at every channel x time1 x time2, Fisher z (earlier
%      trial = time1); same values as spectral_rsa.m, computed faster
%   6. baseline for every Isim and Csim pair (the original method,
%      rsa_subject_level.m correlate_pairs_and_package, "average first,
%      then correlate"): corr(A activation at time1, B's time-averaged
%      baseline spectrum), Fisher z, the same for every time2
%   7. corrected similarity = pair similarity - its baseline
%   8. per-person t-map, corrected Isim vs Csim pairs: FieldTrip
%      ft_freqstatistics, analytic indepsamplesT (as the old pipeline;
%      only used to FIND the cluster at the group level)
%
% Output (C.dir.gammaSubj/<ID>_gamma.mat):
%   simI, simC     corrected similarity per pair  [pair x chan x t1 x t2], single
%   tmap           t-map Isim vs Csim             [chan x t1 x t2]
%   meanI, meanC   mean corrected Isim / Csim     [chan x t1 x t2]
%   rawMeanI, rawMeanC, refMeanI, refMeanC   uncorrected means and the mean
%                  baseline subtracted (corrected = raw - baseline; the drift
%                  check uses refMeanI - refMeanC)
%   PI, PC         pair tables (trials, item/category, condition, Acc, lag,
%                  pair midpoint in the session, mean baseline similarity)
%   T              analysed trials (bebo_trials.m)
%   label, time    channel labels, time axis
% Trial log (C.dir.trialLog/<ID>_gamma_log.csv): counts, lags.

if isnumeric(subj), subj = num2str(subj); end
C = bebo_config();
bebo_init_fieldtrip(C);
for d = {C.dir.gammaSubj, C.dir.trialLog}, if ~exist(d{1}, 'dir'), mkdir(d{1}); end, end

%% 1. Data and trial rule
fprintf('[%s] loading preprocessed data\n', subj);
L = load(fullfile(C.preproc, [subj '_Enc_allcond.mat']), 'data');
data = L.data; clear L
[T, pI, pC] = bebo_trials(data.trialinfo2, C);
fprintf('[%s] %d analysed trials, %d Isim pairs, %d Csim pairs\n', subj, height(T), size(pI,1), size(pC,1));

% keep only the analysed trials; T(k) -> row dataRow(k) of the reduced data
[selRows, ~, dataRow] = unique(T.row);   % sorted rows (ft_selectdata keeps data order)
T.dataRow = dataRow;
cfg = []; cfg.trials = selRows; cfg.channel = C.channel;
data = ft_selectdata(cfg, data);

%% 2. Line noise
cfg = []; cfg.dftfilter = 'yes'; cfg.dftfreq = C.gamma.dftfreq;
data = ft_preprocessing(cfg, data);

%% 3. TFR
cfg = [];
cfg.method     = 'mtmconvol';
cfg.output     = 'pow';
cfg.taper      = 'dpss';
cfg.foi        = C.gamma.foi;
cfg.t_ftimwin  = C.gamma.t_ftimwin * ones(size(cfg.foi));
cfg.tapsmofrq  = C.gamma.tapsmofrq * ones(size(cfg.foi));
cfg.toi        = C.gamma.toi;
cfg.keeptrials = 'yes';
cfg.pad        = C.gamma.pad;
fprintf('[%s] TFR (DPSS %d-%d Hz, +/-%d Hz, %.1f s, %d ms)\n', subj, cfg.foi(1), cfg.foi(end), C.gamma.tapsmofrq, C.gamma.t_ftimwin, round(1000*mean(diff(cfg.toi))));
tfr = ft_freqanalysis(cfg, data);
label = tfr.label; time = tfr.time;
cfg.toi = C.gamma.blToi;                       % pre-stimulus baseline, same settings
tfrBL = ft_freqanalysis(cfg, data);
clear data

%% 4. Aperiodic background (BOSC fit on the analysis band), log10
cfgb = []; cfgb.parameter = 'powspctrm';
warning('off', 'stats:statrobustfit:IterationLimit')
fsel = tfr.freq >= C.gamma.boscBand(1) - 1e-6 & tfr.freq <= C.gamma.boscBand(2) + 1e-6;
tfr.powspctrm = tfr.powspctrm(:, :, fsel, :); tfr.freq = tfr.freq(fsel);
tfr = subtract_mean_noise_spectrum(cfgb, tfr);
P = log10(tfr.powspctrm);                       % trial x chan x freq x time
clear tfr
if any(~isfinite(P(:)))
    error('s1_gamma_subject:nan', '[%s] non-finite power after BOSC/log10', subj);
end

% baseline: BOSC fitted separately (as rsa_subject_level.m), log10, then
% averaged over the baseline time points (trial x chan x freq)
tfrBL.powspctrm = tfrBL.powspctrm(:, :, fsel, :); tfrBL.freq = tfrBL.freq(fsel);
tfrBL = subtract_mean_noise_spectrum(cfgb, tfrBL);
PB = mean(log10(tfrBL.powspctrm), 4);
clear tfrBL
if any(~isfinite(PB(:)))
    error('s1_gamma_subject:nanBL', '[%s] non-finite baseline power after BOSC/log10', subj);
end

% z-score every spectrum over frequency once: then Pearson r = Z1'*Z2/(nf-1)
% (identical to spectral_rsa.m -- corr + atanh -- to rounding, checked on
% simulated data 2026-10-17; faster)
nf = size(P, 3);
Z = (P - mean(P, 3)) ./ std(P, 0, 3);
ZB = (PB - mean(PB, 3)) ./ std(PB, 0, 3);
clear P PB
nCh = numel(label); nT = numel(time);

%% 5. Isim and Csim pairs (uncorrected)
lagI = T.trialNr(pI(:,2)) - T.trialNr(pI(:,1));
lagC = T.trialNr(pC(:,2)) - T.trialNr(pC(:,1));
rawI = zeros([size(pI,1) nCh nT nT], 'single');
for k = 1:size(pI,1), rawI(k,:,:,:) = sim_pair(Z, T.dataRow(pI(k,1)), T.dataRow(pI(k,2))); end
rawC = zeros([size(pC,1) nCh nT nT], 'single');
for k = 1:size(pC,1), rawC(k,:,:,:) = sim_pair(Z, T.dataRow(pC(k,1)), T.dataRow(pC(k,2))); end

%% 6. Baseline (original method): A's activation vs B's baseline spectrum
%  [pair x chan x time1]; identical for every time2
blI = zeros([size(pI,1) nCh nT], 'single');
for k = 1:size(pI,1), blI(k,:,:) = sim_bl(Z, ZB, T.dataRow(pI(k,1)), T.dataRow(pI(k,2))); end
blC = zeros([size(pC,1) nCh nT], 'single');
for k = 1:size(pC,1), blC(k,:,:) = sim_bl(Z, ZB, T.dataRow(pC(k,1)), T.dataRow(pC(k,2))); end
refI = repmat(blI, [1 1 1 nT]); refC = repmat(blC, [1 1 1 nT]);
clear blI blC
infoI = mean(reshape(refI, size(refI,1), []), 2);   % mean baseline similarity per pair
infoC = mean(reshape(refC, size(refC,1), []), 2);

%% 7. Corrected similarity
simI = rawI - refI;
simC = rawC - refC;
rawMeanI = squeeze(mean(rawI, 1)); rawMeanC = squeeze(mean(rawC, 1));
refMeanI = squeeze(mean(refI, 1));  refMeanC = squeeze(mean(refC, 1));
clear rawI rawC refI refC Z ZB

%% 8. Per-person t-map: FieldTrip ft_freqstatistics, analytic indepsamplesT,
%  corrected Isim pairs (condition 1) vs corrected Csim pairs (condition 2),
%  set up as in the old pipeline (rsa_group_level_perAge.m,
%  independent_samples_ttest). The 'freq' axis carries time1 (FieldTrip
%  dimord convention, as before).
nI = size(simI,1); nC = size(simC,1);
if nI < 2 || nC < 2, error('s1_gamma_subject:pairs', '[%s] too few pairs (Isim %d, Csim %d)', subj, nI, nC); end
mk = @(X) struct('dimord', 'rpt_chan_freq_time', 'label', {label(:)}, 'freq', time(:)', ...
    'time', time(:)', 'powspctrm', double(X));
cfg = [];
cfg.channel = 'all'; cfg.latency = 'all'; cfg.frequency = 'all';
cfg.avgoverchan = 'no'; cfg.avgovertime = 'no'; cfg.avgoverfreq = 'no';
cfg.parameter = 'powspctrm';
cfg.method    = 'analytic';
cfg.statistic = 'indepsamplesT';
cfg.tail      = 0;
cfg.design    = [ones(1, nI) 2*ones(1, nC)];
cfg.ivar      = 1;
ft_warning off 'FieldTrip:getdimord:warning_dimord_could_not_be_determined'
st = ft_freqstatistics(cfg, mk(simI), mk(simC));
tmap = st.stat;
clear st
meanI = squeeze(mean(double(simI), 1)); meanC = squeeze(mean(double(simC), 1));

%% Pair tables
PI = table(pI(:,1), pI(:,2), T.item(pI(:,1)), T.category(pI(:,1)), T.exemplar(pI(:,1)), ...
    T.catSize(pI(:,1)), T.itemRep(pI(:,1)), T.Acc(pI(:,1)), lagI, ...
    (T.trialNr(pI(:,1)) + T.trialNr(pI(:,2))) / 2, infoI, ...
    'VariableNames', {'tA','tB','item','category','exemplar','catSize','itemRep','Acc','lag','pos','refSim'});
PC = table(pC(:,1), pC(:,2), T.item(pC(:,1)), T.item(pC(:,2)), T.category(pC(:,1)), ...
    T.exemplar(pC(:,1)), T.exemplar(pC(:,2)), T.catSize(pC(:,1)), T.Acc(pC(:,1)), T.Acc(pC(:,2)), ...
    lagC, (T.trialNr(pC(:,1)) + T.trialNr(pC(:,2))) / 2, infoC, ...
    'VariableNames', {'tA','tB','itemA','itemB','category','exemplarA','exemplarB','catSize','AccA','AccB', ...
    'lag','pos','refSim'});

%% Save
save(fullfile(C.dir.gammaSubj, [subj '_gamma.mat']), 'simI', 'simC', 'tmap', 'meanI', 'meanC', ...
    'rawMeanI', 'rawMeanC', 'refMeanI', 'refMeanC', 'PI', 'PC', 'T', 'label', 'time', '-v7.3');

%% Trial log
tlog = table(string(subj), height(T), sum(T.pres==1), sum(T.Acc(T.pres==1)), sum(1-T.Acc(T.pres==1)), nI, nC, ...
    median(lagI), min(lagI), max(lagI), median(lagC), min(lagC), max(lagC), ...
    'VariableNames', {'ID','nTrials','nItems','nOld','nSimilar','nIsim','nCsim', ...
    'lagI_median','lagI_min','lagI_max','lagC_median','lagC_min','lagC_max'});
writetable(tlog, fullfile(C.dir.trialLog, [subj '_gamma_log.csv']));
fprintf('[%s] done: Isim %d, Csim %d pairs; lags Isim %d-%d, Csim %d-%d; mean baseline similarity Isim %.4f, Csim %.4f\n', ...
    subj, nI, nC, min(lagI), max(lagI), min(lagC), max(lagC), mean(infoI), mean(infoC));
end

function S = sim_bl(Z, ZB, a, b)
% baseline similarity of pair (a, b): a's activation spectrum at every time
% point vs b's time-averaged baseline spectrum -> chan x t1, Fisher z
% (= spectral_rsa(tf1, repmat(mean(tf2bl,2), ...)) of the original, whose
% columns are all equal)
[~, nCh, nf, nT] = size(Z);
S = zeros(nCh, nT, 'single');
for ch = 1:nCh
    r = reshape(Z(a, ch, :, :), nf, nT)' * reshape(ZB(b, ch, :), nf, 1) / (nf - 1);
    r = min(max(r, -0.999999), 0.999999);
    S(ch, :) = single(atanh(r));
end
end

function S = sim_pair(Z, a, b)
% similarity of data rows a and b: chan x t1 x t2, Fisher z. Z holds the
% spectra z-scored over frequency, so r = Za' * Zb / (nf - 1) -- the same
% Pearson r as spectral_rsa.m (corr + atanh), computed faster.
[~, nCh, nf, nT] = size(Z);
S = zeros(nCh, nT, nT, 'single');
for ch = 1:nCh
    r = reshape(Z(a, ch, :, :), nf, nT)' * reshape(Z(b, ch, :, :), nf, nT) / (nf - 1);
    r = min(max(r, -0.999999), 0.999999);
    S(ch, :, :) = single(atanh(r));
end
end

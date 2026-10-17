function s1_gamma_subject(subj)
% S1_GAMMA_SUBJECT  Per-person gamma RSA with lag-matched reference
% (analysis plan 1a). Computes everything from the preprocessed data.
%
%   1. load PreprocData/<ID>_Enc_allcond.mat, keep the analysed trials
%      (bebo_trials.m) and the 60 EEG channels
%   2. line noise: DFT filter at 50 and 100 Hz
%   3. TFR: DPSS mtmconvol, 35:2:100 Hz, +/-5 Hz, 0.4 s, 20 ms, 0-1.7 s
%   4. aperiodic background removed (BOSC fit on 35-100 Hz), log10
%   5. similarity of every Isim and Csim pair: Pearson r between the two
%      trials' spectra at every channel x time1 x time2, Fisher z
%      (earlier trial = time1)
%   6. reference for every Isim and Csim pair: mean similarity of the 20
%      reference pairs (two analysed trials of DIFFERENT categories) with
%      the same lag whose midpoint is closest in the session (lags +/-1,
%      +/-2 added if fewer than 20 exist) -- matched on lag AND position
%   7. corrected similarity = pair similarity - its reference
%   8. per-person independent-samples t-map, corrected Isim vs Csim pairs
%      (only used to FIND the cluster at the group level)
%
% Output (C.dir.gammaSubj/<ID>_gamma.mat):
%   simI, simC     corrected similarity per pair  [pair x chan x t1 x t2], single
%   tmap           t-map Isim vs Csim             [chan x t1 x t2]
%   meanI, meanC   mean corrected Isim / Csim     [chan x t1 x t2]
%   rawMeanI, rawMeanC, refMeanI, refMeanC   uncorrected means and the mean
%                  reference subtracted (corrected = raw - ref; the drift
%                  check uses refMeanI - refMeanC)
%   PI, PC         pair tables (trials, item/category, condition, Acc, lag,
%                  position, n reference pairs, lag widening, position
%                  distance of the farthest reference, mean reference similarity)
%   T              analysed trials (bebo_trials.m)
%   label, time    channel labels, time axis
% Trial log (C.dir.trialLog/<ID>_gamma_log.csv): counts, lags, references.

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

% z-score every spectrum over frequency once: then r = Z1'*Z2/(nf-1)
nf = size(P, 3);
Z = (P - mean(P, 3)) ./ std(P, 0, 3);
clear P
nCh = numel(label); nT = numel(time);

%% 5. Isim and Csim pairs (uncorrected)
lagI = T.trialNr(pI(:,2)) - T.trialNr(pI(:,1));
lagC = T.trialNr(pC(:,2)) - T.trialNr(pC(:,1));
rawI = zeros([size(pI,1) nCh nT nT], 'single');
for k = 1:size(pI,1), rawI(k,:,:,:) = sim_pair(Z, T.dataRow(pI(k,1)), T.dataRow(pI(k,2))); end
rawC = zeros([size(pC,1) nCh nT nT], 'single');
for k = 1:size(pC,1), rawC(k,:,:,:) = sim_pair(Z, T.dataRow(pC(k,1)), T.dataRow(pC(k,2))); end

%% 6. Reference matched on lag AND session position (per pair)
% Candidate reference pairs: all pairs of analysed trials from different
% categories, earlier trial first. For each Isim/Csim pair: the nRef
% candidates with the same lag (widening to +/-1, +/-2 if fewer exist)
% whose midpoint in the session is CLOSEST to the pair's midpoint.
% Matching position as well as lag is needed because slow changes over a
% session (fatigue, impedance, ...) are not uniform: a lag-only reference
% left a bias that followed session position in the simulations
% (test/check_lag_correction.m). Deterministic -- no random draw.
nTr = height(T);
[ia, ib] = find(triu(true(nTr), 1));
okCat = T.category(ia) ~= T.category(ib);
ia = ia(okCat); ib = ib(okCat);
refLag = T.trialNr(ib) - T.trialNr(ia);
refMid = (T.trialNr(ia) + T.trialNr(ib)) / 2;

    function [ref, info] = local_reference(pairs, lagsP)
        ref = zeros([size(pairs,1) nCh nT nT], 'single');
        info = zeros(size(pairs,1), 4);   % nRef, widen, max position distance, mean ref similarity
        for k = 1:size(pairs,1)
            mid = (T.trialNr(pairs(k,1)) + T.trialNr(pairs(k,2))) / 2;
            for w = 0:C.lag.maxWiden
                cand = find(abs(refLag - lagsP(k)) <= w);
                if numel(cand) >= C.lag.nRef, break, end
            end
            if isempty(cand)
                error('s1_gamma_subject:noRef', '[%s] no reference pair for lag %d', subj, lagsP(k));
            end
            [dist, o] = sort(abs(refMid(cand) - mid));
            n = min(C.lag.nRef, numel(cand));
            sel = cand(o(1:n));
            acc = zeros(nCh, nT, nT);
            for q = 1:n, acc = acc + double(sim_pair(Z, T.dataRow(ia(sel(q))), T.dataRow(ib(sel(q))))); end
            ref(k,:,:,:) = single(acc / n);
            info(k,:) = [n, w, dist(n), mean(acc(:)) / n];
        end
    end

[refI, infoI] = local_reference(pI, lagI);
[refC, infoC] = local_reference(pC, lagC);
if any([infoI(:,1); infoC(:,1)] < C.lag.nRef)
    warning('s1_gamma_subject:fewRef', '[%s] %d pair(s) with fewer than %d reference pairs', subj, sum([infoI(:,1); infoC(:,1)] < C.lag.nRef), C.lag.nRef);
end

%% 7. Corrected similarity
simI = rawI - refI;
simC = rawC - refC;
rawMeanI = squeeze(mean(rawI, 1)); rawMeanC = squeeze(mean(rawC, 1));
refMeanI = squeeze(mean(refI, 1));  refMeanC = squeeze(mean(refC, 1));
clear rawI rawC refI refC Z

%% 8. Per-person t-map (pooled-variance independent-samples t, = FieldTrip indepsamplesT)
nI = size(simI,1); nC = size(simC,1);
if nI < 2 || nC < 2, error('s1_gamma_subject:pairs', '[%s] too few pairs (Isim %d, Csim %d)', subj, nI, nC); end
meanI = squeeze(mean(double(simI), 1)); meanC = squeeze(mean(double(simC), 1));
vI = squeeze(var(double(simI), 0, 1)); vC = squeeze(var(double(simC), 0, 1));
sp = sqrt(((nI-1)*vI + (nC-1)*vC) / (nI+nC-2));
tmap = (meanI - meanC) ./ (sp * sqrt(1/nI + 1/nC));

%% Pair tables
PI = table(pI(:,1), pI(:,2), T.item(pI(:,1)), T.category(pI(:,1)), T.exemplar(pI(:,1)), ...
    T.catSize(pI(:,1)), T.itemRep(pI(:,1)), T.Acc(pI(:,1)), lagI, T.trialNr(pI(:,1)), ...
    infoI(:,1), infoI(:,2), infoI(:,3), infoI(:,4), ...
    'VariableNames', {'tA','tB','item','category','exemplar','catSize','itemRep','Acc','lag','posA', ...
    'nRef','widen','refPosDist','refSim'});
PC = table(pC(:,1), pC(:,2), T.item(pC(:,1)), T.item(pC(:,2)), T.category(pC(:,1)), ...
    T.exemplar(pC(:,1)), T.exemplar(pC(:,2)), T.catSize(pC(:,1)), T.Acc(pC(:,1)), T.Acc(pC(:,2)), ...
    lagC, T.trialNr(pC(:,1)), infoC(:,1), infoC(:,2), infoC(:,3), infoC(:,4), ...
    'VariableNames', {'tA','tB','itemA','itemB','category','exemplarA','exemplarB','catSize','AccA','AccB', ...
    'lag','posA','nRef','widen','refPosDist','refSim'});

%% Save
save(fullfile(C.dir.gammaSubj, [subj '_gamma.mat']), 'simI', 'simC', 'tmap', 'meanI', 'meanC', ...
    'rawMeanI', 'rawMeanC', 'refMeanI', 'refMeanC', 'PI', 'PC', 'T', 'label', 'time', '-v7.3');

%% Trial log
tlog = table(string(subj), height(T), sum(T.pres==1), sum(T.Acc(T.pres==1)), sum(1-T.Acc(T.pres==1)), nI, nC, ...
    median(lagI), min(lagI), max(lagI), median(lagC), min(lagC), max(lagC), ...
    min([infoI(:,1); infoC(:,1)]), sum([infoI(:,2); infoC(:,2)] > 0), median([infoI(:,3); infoC(:,3)]), max([infoI(:,3); infoC(:,3)]), ...
    'VariableNames', {'ID','nTrials','nItems','nOld','nSimilar','nIsim','nCsim', ...
    'lagI_median','lagI_min','lagI_max','lagC_median','lagC_min','lagC_max', ...
    'minRefPerPair','nPairsWidened','refPosDist_median','refPosDist_max'});
writetable(tlog, fullfile(C.dir.trialLog, [subj '_gamma_log.csv']));
fprintf('[%s] done: Isim %d, Csim %d pairs; lags Isim %d-%d, Csim %d-%d; references: min %d per pair, position distance median %.0f, max %.0f\n', ...
    subj, nI, nC, min(lagI), max(lagI), min(lagC), max(lagC), min([infoI(:,1); infoC(:,1)]), median([infoI(:,3); infoC(:,3)]), max([infoI(:,3); infoC(:,3)]));
end

function S = sim_pair(Z, a, b)
% similarity of data rows a and b: chan x t1 x t2, Fisher z. Z holds the
% spectra z-scored over frequency, so r = Za' * Zb / (nf - 1).
[~, nCh, nf, nT] = size(Z);
S = zeros(nCh, nT, nT, 'single');
for ch = 1:nCh
    r = reshape(Z(a, ch, :, :), nf, nT)' * reshape(Z(b, ch, :, :), nf, nT) / (nf - 1);
    r = min(max(r, -0.999999), 0.999999);
    S(ch, :, :) = single(atanh(r));
end
end

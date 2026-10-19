function test_lag_effects_onesubj(subj, freqLo, freqHi, bgFitMode)
% DIAGNOSTIC ONLY -- does not touch the production pipeline. Characterizes
% whether trial lag (how far apart two trials are in the sequence)
% predicts corrected similarity, independent of condition (2026-08-27,
% see session notes: is item sim's apparent specificity partly a lag
% artifact, since item pairs are inherently close together and between-
% cat pairs inherently far apart?).
%
% Uses ONLY trimmed, shape-matched data -- activation [0 1.5]s and
% baseline [-1.5 0]s, both the same duration -- so every correlation here
% is a literal, single-timepoint-vs-single-timepoint comparison with NO
% averaging-order machinery (no OLD/NEW vector-vs-correlation averaging
% question at all; see test_baseline_correction_onesubj.m for that
% separate question). Per user request 2026-08-27: keep this one simple.
%
% For all three pairing types (within-item, within-cat, between-cat --
% between-cat pairs kept individually, NOT averaged across partners like
% production does, so we get more (lag, similarity) data points), computes
% per pair, averaged over channels and the full trimmed time-time grid
% into ONE scalar per pair:
%   lag        = |trialinfo2(HatA,1) - trialinfo2(HatB,1)| -- column 1 of
%                trialinfo2 is the actual trial number (confirmed
%                2026-08-27), so this is the true trial-sequence lag, not
%                a row-index proxy.
%   real       = mean corr(activ HatA, activ HatB)
%   actVsBL    = mean corr(activ HatA, baseline HatB)     [current-style correction term]
%   blVsBl     = mean corr(baseline HatA, baseline HatB)  [proposed correction term]
%   corrected_actVsBL = real - actVsBL
%   corrected_blVsBl  = real - blVsBl
%
% Reports BOTH the RAW (uncorrected) and corrected lag relationships,
% within condition and pooled -- the raw one is the actual "is there a
% lag effect at all" question; the corrected ones show whether either
% correction term removes it. Saves the per-pair table to
% Results/BaselineDiagnostic/<subj>_lag_effects.mat/.txt.
%
% ALSO builds a full-sample (all valid trials, any condition), trial-
% order-sorted raw similarity matrix and a binned lag-decay curve, both
% saved as .png/.mat. SIMPLIFIED for speed: each trial is collapsed to
% one frequency-vector (averaged over channels AND time) before
% correlating, rather than the single-timepoint, no-averaging approach
% used for the condition-specific pairs above -- this reintroduces the
% same "average-then-correlate" concern from the OLD/NEW baseline-method
% debate, so treat this as a fast first-look visual/quantitative check,
% not a rigorous measure (agreed with user 2026-08-27).
%
% Usage: matlab -r "test_lag_effects_onesubj('3109','30','100','legacy'); exit"

if nargin < 1 || isempty(subj), error('subj is required, e.g. ''3109'''); end
if nargin < 2 || isempty(freqLo), freqLo = 30; end
if nargin < 3 || isempty(freqHi), freqHi = 100; end
if nargin < 4 || isempty(bgFitMode), bgFitMode = 'legacy'; end
if ischar(freqLo), freqLo = str2double(freqLo); end
if ischar(freqHi), freqHi = str2double(freqHi); end
subj = num2str(subj);

pipelineRoot = fileparts(mfilename('fullpath'));
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/');
ft_defaults

pdat.tfr   = fullfile(pipelineRoot, '..', 'TFRdata', 'Time-frequencyEEGdata');
pdat.tfrBL = fullfile(pipelineRoot, '..', 'TFRdat_BL');

if strcmp(bgFitMode, 'legacy')
    bgFitBand = [freqLo freqHi];
elseif freqHi <= 20
    bgFitBand = [2 20];
else
    bgFitBand = [25 125];
end
analysisBand = [freqLo freqHi];

fprintf('Loading + BOSC-correcting subject %s, TRIMMED windows only (%d-%dHz, bgFitMode=%s)...\n', subj, freqLo, freqHi, bgFitMode);
Data_trim = load_activation_trimmed(subj, pdat, bgFitBand, analysisBand);
Data_bl   = load_baseline_trimmed(subj, pdat, bgFitBand, analysisBand);

nTimeTrim = length(Data_trim.time);
nTimeBL = length(Data_bl.time);
if nTimeTrim ~= nTimeBL
    error(['Trimmed activation [0 1.5]s has %d timepoints but wide baseline ' ...
        '[-1.5 0]s has %d -- different sampling grids, literal comparisons not valid.'], nTimeTrim, nTimeBL);
end
fprintf('Trimmed windows: %d timepoints each, %d channels\n', nTimeTrim, size(Data_trim.powspctrm,2));

%% Build pairs for all three conditions

pairsItem = pairs_within_item(Data_trim);
pairsCat  = pairs_within_cat(Data_trim);
[pairsBw1, pairsBw2] = pairs_between_cat(Data_trim);

fprintf('Pairs: within-item=%d, within-cat=%d, between-cat=%d\n', ...
    size(pairsItem,1), size(pairsCat,1), length(pairsBw1));

%% Compute per-pair stats for each condition

T_item = compute_condition_table(Data_trim, Data_bl, pairsItem(:,1), pairsItem(:,2), 'within-item');
T_cat  = compute_condition_table(Data_trim, Data_bl, pairsCat(:,1),  pairsCat(:,2),  'within-cat');
T_bw   = compute_condition_table(Data_trim, Data_bl, pairsBw1,       pairsBw2,       'between-cat');

T = [T_item; T_cat; T_bw];

%% Report lag-vs-corrected-similarity relationship, within condition and pooled

fprintf('\n=== LAG EFFECTS ON RAW (uncorrected) SIMILARITY -- subject %s, %d-%dHz, %s ===\n', subj, freqLo, freqHi, bgFitMode);
fprintf('(this is the actual question: does lag predict similarity BEFORE any correction is applied?)\n');
conds = {'within-item','within-cat','between-cat'};
for c = 1:length(conds)
    rows = strcmp(T.condition, conds{c});
    lag_c = T.lag(rows);
    r_real = corr(lag_c, T.realVal(rows), 'rows','pairwise');
    r_actRaw = corr(lag_c, T.actVsBL(rows), 'rows','pairwise');
    r_blRaw  = corr(lag_c, T.blVsBl(rows),  'rows','pairwise');
    fprintf('%-14s n=%4d  lag range [%4d,%4d]  corr(lag,real)=%6.3f  corr(lag,actVsBL)=%6.3f  corr(lag,blVsBl)=%6.3f\n', ...
        conds{c}, sum(rows), min(lag_c), max(lag_c), r_real, r_actRaw, r_blRaw);
end
r_real_all = corr(T.lag, T.realVal, 'rows','pairwise');
r_actRaw_all = corr(T.lag, T.actVsBL, 'rows','pairwise');
r_blRaw_all  = corr(T.lag, T.blVsBl,  'rows','pairwise');
fprintf('%-14s n=%4d  lag range [%4d,%4d]  corr(lag,real)=%6.3f  corr(lag,actVsBL)=%6.3f  corr(lag,blVsBl)=%6.3f\n', ...
    'POOLED (all)', height(T), min(T.lag), max(T.lag), r_real_all, r_actRaw_all, r_blRaw_all);

fprintf('\nMean RAW (uncorrected) values by condition:\n');
for c = 1:length(conds)
    rows = strcmp(T.condition, conds{c});
    fprintf('%-14s real mean=%.4f  actVsBL mean=%.4f  blVsBl mean=%.4f\n', ...
        conds{c}, mean(T.realVal(rows),'omitnan'), mean(T.actVsBL(rows),'omitnan'), mean(T.blVsBl(rows),'omitnan'));
end

fprintf('\n=== LAG EFFECTS ON CORRECTED SIMILARITY (after subtracting each baseline term) ===\n');
for c = 1:length(conds)
    rows = strcmp(T.condition, conds{c});
    lag_c = T.lag(rows);
    ca = T.corrected_actVsBL(rows);
    cb = T.corrected_blVsBl(rows);
    r_act = corr(lag_c, ca, 'rows','pairwise');
    r_bl  = corr(lag_c, cb, 'rows','pairwise');
    fprintf('%-14s n=%4d  lag range [%4d,%4d]  corr(lag,corrected_actVsBL)=%6.3f  corr(lag,corrected_blVsBl)=%6.3f\n', ...
        conds{c}, sum(rows), min(lag_c), max(lag_c), r_act, r_bl);
end
r_act_all = corr(T.lag, T.corrected_actVsBL, 'rows','pairwise');
r_bl_all  = corr(T.lag, T.corrected_blVsBl,  'rows','pairwise');
fprintf('%-14s n=%4d  lag range [%4d,%4d]  corr(lag,corrected_actVsBL)=%6.3f  corr(lag,corrected_blVsBl)=%6.3f\n', ...
    'POOLED (all)', height(T), min(T.lag), max(T.lag), r_act_all, r_bl_all);

fprintf('\nMean corrected similarity by condition:\n');
for c = 1:length(conds)
    rows = strcmp(T.condition, conds{c});
    fprintf('%-14s corrected_actVsBL mean=%.4f  corrected_blVsBl mean=%.4f\n', ...
        conds{c}, mean(T.corrected_actVsBL(rows),'omitnan'), mean(T.corrected_blVsBl(rows),'omitnan'));
end

%% Full-sample (all valid trials) raw similarity matrix, ordered by true
%  trial number -- a quick visual check for a lag/diagonal-banding
%  pattern in RAW similarity, independent of any pairing condition.
%  SIMPLIFIED for speed: each trial collapsed to ONE frequency-vector
%  (averaged over time AND channels) rather than the full per-channel
%  time-time treatment used above -- this is a visual diagnostic, not a
%  final measure, so the full N-trials-squared pairwise comparison is
%  done as a single corr() call instead of looping.

fprintf('\nBuilding full trial-order raw similarity matrix...\n');
validTrials = find(strcmp(Data_trim.trialinfo2(:,11),'old') | strcmp(Data_trim.trialinfo2(:,11),'similar'));
trialNrAll = cell2mat(Data_trim.trialinfo2(validTrials,1));
[trialNrSorted, sortIdx] = sort(trialNrAll);
validTrialsSorted = validTrials(sortIdx);

nValid = length(validTrialsSorted);
vecs = nan(size(Data_trim.powspctrm,3), nValid); % [nFreq x nValidTrials]
for t = 1:nValid
    tfr = squeeze(Data_trim.powspctrm(validTrialsSorted(t),:,:,:)); % chan x freq x time
    vecs(:,t) = squeeze(mean(mean(tfr,1,'omitnan'),3,'omitnan'));   % avg over chan and time -> [nFreq x 1]
end
fullSimMat = corr(vecs, 'rows','pairwise');

fig = figure('Visible','off','Units','centimeters','Position',[2 2 20 18]);
imagesc(fullSimMat); axis square; colorbar;
xlabel('Trial (sorted by trial number)'); ylabel('Trial (sorted by trial number)');
title(sprintf('Subject %s: raw activation similarity, trial-order sorted (n=%d valid trials)', subj, nValid));
outDir = fullfile(pipelineRoot, 'Results', 'BaselineDiagnostic');
if ~exist(outDir,'dir'), mkdir(outDir); end
print(fig, fullfile(outDir, [subj '_fullTrialOrder_simMatrix.png']), '-dpng', '-r200');
fprintf('Saved trial-order similarity matrix plot -> %s\n', fullfile(outDir, [subj '_fullTrialOrder_simMatrix.png']));
save(fullfile(outDir, [subj '_fullTrialOrder_simMatrix.mat']), 'fullSimMat', 'trialNrSorted', 'validTrialsSorted', '-v7.3');

%% Quantify similarity as a function of lag, using the same full matrix --
%  every unique trial pair (not just condition-matched pairs), so this is
%  the most representative, condition-agnostic look at whether lag
%  predicts raw similarity at all.

[I,J] = find(triu(true(nValid),1)); % every unique pair, excluding the trivial diagonal
lagFull = abs(trialNrSorted(I) - trialNrSorted(J));
simFull = fullSimMat(sub2ind([nValid nValid], I, J));

r_lag_full = corr(lagFull, simFull, 'rows','pairwise');
fprintf('\nFull-sample (all %d valid-trial pairs) corr(lag, raw similarity) = %.3f\n', length(lagFull), r_lag_full);

edges = [0 5 10 20 40 80 160 320 640 max(lagFull)+1];
binIdx = discretize(lagFull, edges);
nBins = max(binIdx);
binMean = nan(nBins,1); binSE = nan(nBins,1); binN = nan(nBins,1); binCenter = nan(nBins,1);
for b = 1:nBins
    inBin = binIdx==b;
    binMean(b) = mean(simFull(inBin),'omitnan');
    binSE(b) = std(simFull(inBin),'omitnan')/sqrt(sum(inBin));
    binN(b) = sum(inBin);
    binCenter(b) = mean(lagFull(inBin),'omitnan');
end

fprintf('\nMean raw similarity by lag bin:\n');
fprintf('%12s %8s %10s %10s\n', 'lag~center', 'n', 'mean', 'SE');
for b = 1:nBins
    fprintf('%12.1f %8d %10.4f %10.4f\n', binCenter(b), binN(b), binMean(b), binSE(b));
end

fig2 = figure('Visible','off','Units','centimeters','Position',[2 2 20 14]);
errorbar(binCenter, binMean, binSE, '-o','LineWidth',1.2,'MarkerFaceColor','auto');
xlabel('Trial lag'); ylabel('Mean raw similarity (\pm SE)');
title(sprintf('Subject %s: raw similarity vs. trial lag (n=%d pairs)', subj, length(lagFull)));
box off;
print(fig2, fullfile(outDir, [subj '_lagDecayCurve.png']), '-dpng', '-r200');
fprintf('Saved lag decay curve plot -> %s\n', fullfile(outDir, [subj '_lagDecayCurve.png']));

save(fullfile(outDir, [subj '_lagDecayCurve.mat']), ...
    'lagFull','simFull','r_lag_full','binMean','binSE','binN','binCenter','-v7.3');

save(fullfile(outDir, [subj '_lag_effects.mat']), 'T', '-v7.3');
writetable(T, fullfile(outDir, [subj '_lag_effects.txt']));
fprintf('\nSaved to %s\n', fullfile(outDir, [subj '_lag_effects.mat']));

end


%% =========================================================
%% LOCAL FUNCTIONS
%% =========================================================

function T = compute_condition_table(Data_trim, Data_bl, stim1list, stim2list, condName)
n = length(stim1list);
lag = nan(n,1); realVal = nan(n,1); actVsBL = nan(n,1); blVsBl = nan(n,1);
nChan = size(Data_trim.powspctrm,2);
trialNr = cell2mat(Data_trim.trialinfo2(:,1)); % actual trial number, confirmed 2026-08-27

for p = 1:n
    s1 = stim1list(p); s2 = stim2list(p);
    lag(p) = abs(trialNr(s1) - trialNr(s2));

    realC = nan(nChan,1); actBlC = nan(nChan,1); blBlC = nan(nChan,1);
    for chan = 1:nChan
        tf1 = squeeze(Data_trim.powspctrm(s1,chan,:,:));
        tf2 = squeeze(Data_trim.powspctrm(s2,chan,:,:));
        r = spectral_rsa([],tf1,tf2,Data_trim.time);
        realC(chan) = mean(r(:),'omitnan');

        tf2bl = squeeze(Data_bl.powspctrm(s2,chan,:,:));
        r2 = spectral_rsa([],tf1,tf2bl,Data_trim.time);
        actBlC(chan) = mean(r2(:),'omitnan');

        tf1bl = squeeze(Data_bl.powspctrm(s1,chan,:,:));
        r3 = spectral_rsa([],tf1bl,tf2bl,Data_bl.time);
        blBlC(chan) = mean(r3(:),'omitnan');
    end
    realVal(p) = mean(realC,'omitnan');
    actVsBL(p) = mean(actBlC,'omitnan');
    blVsBl(p)  = mean(blBlC,'omitnan');

    if mod(p,25)==0, fprintf('  [%s] processed %d/%d pairs\n', condName, p, n); end
end

condition = repmat({condName}, n, 1);
corrected_actVsBL = realVal - actVsBL;
corrected_blVsBl  = realVal - blVsBl;
T = table(condition, lag, realVal, actVsBL, blVsBl, corrected_actVsBL, corrected_blVsBl);

end


function pairs = pairs_within_item(Data_trim)
% Identical logic to rsa_subject_level.m's within-item pairing.
pres1 = 1;
pairs = [];
for stim1 = 1:size(Data_trim.powspctrm,1)
    if ~ismember(stim1,pairs) && Data_trim.trialinfo2{stim1,9} == pres1 && ~isempty(find(strcmp(Data_trim.trialinfo2{stim1,11},'old') | strcmp(Data_trim.trialinfo2{stim1,11},'similar')))
        curr_pair = Data_trim.trialinfo2{stim1,7};
        if strcmp(Data_trim.trialinfo2(stim1,4),'normal') || strcmp(Data_trim.trialinfo2(stim1,4),'high-exemplar')
            if find(strcmp(Data_trim.trialinfo2(:,7),curr_pair) & cell2mat(Data_trim.trialinfo2(:,9)) == 2)
                pairs(stim1,2) = find(strcmp(Data_trim.trialinfo2(:,7),curr_pair) & cell2mat(Data_trim.trialinfo2(:,9)) == 2);
                pairs(stim1,1) = stim1;
            end
        elseif strcmp(Data_trim.trialinfo2(stim1,4),'high-repetition+exemplar') || strcmp(Data_trim.trialinfo2(stim1,4),'high-repetition')
            if find(strcmp(Data_trim.trialinfo2(:,7),curr_pair) & cell2mat(Data_trim.trialinfo2(:,9)) == 2)
                pairs(stim1,1) = stim1;
                pairs(stim1,2) = find(strcmp(Data_trim.trialinfo2(:,7),curr_pair) & cell2mat(Data_trim.trialinfo2(:,9)) == 2);
            end
        end
    end
end
pairs(~all(pairs,2),:) = [];
end


function pairs = pairs_within_cat(Data_trim)
% Identical logic to rsa_subject_level.m's pairs_within_cat_exemplar1v2.
pairs = [];
for stim1 = 1:size(Data_trim.powspctrm,1)
    if ~ismember(stim1,pairs) && cell2mat(Data_trim.trialinfo2(stim1,8)) == 1 && cell2mat(Data_trim.trialinfo2(stim1,9)) == 1 && ~isempty(find(strcmp(Data_trim.trialinfo2{stim1,11},'old') | strcmp(Data_trim.trialinfo2{stim1,11},'similar')))
        curr_pair = Data_trim.trialinfo2{stim1,5};
        stim2 = find(strcmp(Data_trim.trialinfo2(:,5),curr_pair) & cell2mat(Data_trim.trialinfo2(:,8)) == 2 & cell2mat(Data_trim.trialinfo2(:,9)) == 1);
        if ~isempty(stim2)
            pairs(stim1,1) = stim1;
            pairs(stim1,2) = stim2(1);
        end
    end
end
pairs(~all(pairs,2),:) = [];
end


function [pairs1, pairs2] = pairs_between_cat(Data_trim)
% Same candidate-finding logic as rsa_subject_level.m's between-cat
% section, but kept as INDIVIDUAL pairs (not averaged across partners
% per reference item) -- gives more (lag, similarity) data points for
% this lag-characterization diagnostic.
pairs1 = [];
pairs2 = [];
for stim1 = 1:size(Data_trim.powspctrm,1)
    if ~ismember(stim1,pairs1) && cell2mat(Data_trim.trialinfo2(stim1,8)) == 1 && cell2mat(Data_trim.trialinfo2(stim1,9)) == 1 && ~isempty(find(strcmp(Data_trim.trialinfo2{stim1,11},'old') | strcmp(Data_trim.trialinfo2{stim1,11},'similar')))
        curr_pair = Data_trim.trialinfo2{stim1,5};
        stim2 = find(~strcmp(Data_trim.trialinfo2(:,5),curr_pair) & cell2mat(Data_trim.trialinfo2(:,8)) == 2 & cell2mat(Data_trim.trialinfo2(:,9)) == 1 & (strcmp(Data_trim.trialinfo2(:,11),'old') | strcmp(Data_trim.trialinfo2(:,11),'similar')));
        if ~isempty(stim2)
            for s = 1:length(stim2)
                pairs1(end+1,1) = stim1; %#ok<AGROW>
                pairs2(end+1,1) = stim2(s); %#ok<AGROW>
            end
        end
    end
end
end


function Data_logpow = load_activation_trimmed(subj, pdat, bgFitBand, analysisBand)
% Simplified vs. rsa_subject_level.m's load_activation: fits BOSC and
% restricts latency directly to [0 1.5]s from the start (no [0 1.7]-then-
% trim two-step), since this diagnostic only ever needs the trimmed
% window. Minor deviation from production fidelity, acceptable here per
% user request to keep this simple.

disp('loading activation TFR...')
load(fullfile(pdat.tfr, [subj '_Enc_tfr_2_125Hz_ST.mat']), 'tfr')

cfg = [];
cfg.latency = [0 1.5];
cfg.frequency = bgFitBand;
tfr = ft_selectdata(cfg, tfr);

cfg = []; cfg.parameter = 'powspctrm'; cfg.operation = 'log10';
warning('off', 'stats:statrobustfit:IterationLimit')
disp('running BOSC (activation)...')
Data_c = subtract_mean_noise_spectrum(cfg, tfr);
clear tfr

ft_warning off 'FieldTrip:getdimord:warning_dimord_could_not_be_determined'
Data_logpow_wide = ft_math(cfg, Data_c);
Data_logpow_wide.trialinfo2 = Data_c.trialinfo2;
clear Data_c

cfg = []; cfg.frequency = analysisBand;
Data_logpow = ft_selectdata(cfg, Data_logpow_wide);
Data_logpow.trialinfo2 = Data_logpow_wide.trialinfo2;

end


function Data_logpow_bl = load_baseline_trimmed(subj, pdat, bgFitBand, analysisBand)
% Same simplification for baseline: fits BOSC and restricts latency
% directly to the wide [-1.5 0]s window.

disp('loading baseline TFR...')
load(fullfile(pdat.tfrBL, [subj '_Enc_tfrBL_2_125Hz_ST.mat']), 'tfr')

cfg = [];
cfg.latency = [-1.5 0];
cfg.frequency = bgFitBand;
tfr = ft_selectdata(cfg, tfr);

cfg = []; cfg.parameter = 'powspctrm'; cfg.operation = 'log10';
warning('off', 'stats:statrobustfit:IterationLimit')
disp('running BOSC (baseline)...')
Data_c = subtract_mean_noise_spectrum(cfg, tfr);
clear tfr

ft_warning off 'FieldTrip:getdimord:warning_dimord_could_not_be_determined'
Data_logpow_bl_wide = ft_math(cfg, Data_c);
Data_logpow_bl_wide.trialinfo2 = Data_c.trialinfo2;
clear Data_c

cfg = []; cfg.frequency = analysisBand;
Data_logpow_bl = ft_selectdata(cfg, Data_logpow_bl_wide);
Data_logpow_bl.trialinfo2 = Data_logpow_bl_wide.trialinfo2;

end

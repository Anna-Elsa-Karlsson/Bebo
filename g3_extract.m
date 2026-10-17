function g3_extract()
% G3_EXTRACT  Single-trial tables and split-half values (analysis plan 1e, 2).
% Needs g1 (Ispec cluster). Uses ONLY the Ispec cluster for the trial-level
% measures (never the SME cluster: that would make power -> memory circular).
%
% Per person:
%   similarity per pair = mean corrected similarity over all Ispec cluster
%     points
%   trial power = mean log10 power over the Ispec cluster channels x 8-25 Hz
%     x the Ispec cluster's time extent (alpha/beta, a priori band)
%   split half: categories split at random (fixed seed per person); Ispec in
%     each half = mean Isim pairs - mean Csim pairs of that half's categories
%
% Output (C.dir.tables):
%   items.csv   one row per item with an Isim pair (memory model, condition
%               models): ID, YA, category, item, exemplar, catSize, itemRep,
%               Acc, Isim, Csim (mean over the item's Csim pairs; NaN if
%               none), nCsimPairs, Power1st, Power2nd, lagIsim
%   pairs.csv   one row per Isim or Csim pair (pair model): ID, YA, Type,
%               Sim, PowerPair (mean of the two trials), PowerA, PowerB,
%               category, catSize, itemRep, lag
%   person_split.csv   ID, YA, Ispec_half1, Ispec_half2, Ispec_pairs (check)
%   lags.csv    one row per Isim/Csim pair: lag, position, reference position
%               distance, mean reference similarity, lag widening

C = bebo_config();
if ~exist(C.dir.tables, 'dir'), mkdir(C.dir.tables); end
ids = bebo_ids('all'); nYA = numel(bebo_ids('ya')); n = numel(ids);
G = load(fullfile(C.dir.group, 'ispec_cluster.mat'), 'mask', 'clusterChannels', 'timeExtent', 'P');
fprintf('Ispec cluster: %d channels, %.2f-%.2f s; power band %d-%d Hz\n', numel(G.clusterChannels), G.timeExtent, C.alpha.glmBand);

items = cell(n,1); pairs = cell(n,1); lags = cell(n,1);
split = table(ids(:), (1:n)' <= nYA, nan(n,1), nan(n,1), nan(n,1), ...
    'VariableNames', {'ID','YA','Ispec_half1','Ispec_half2','Ispec_pairs'});

for s = 1:n
    id = ids(s); isYA = s <= nYA;
    Lg = load(fullfile(C.dir.gammaSubj, sprintf('%d_gamma.mat', id)), 'simI', 'simC', 'PI', 'PC', 'T');
    La = load(fullfile(C.dir.alphaSubj, sprintf('%d_alpha.mat', id)), 'pow', 'T', 'label', 'freq', 'time');
    if ~isequal(Lg.T.trialNr, La.T.trialNr)
        error('g3_extract:trials', 'Subject %d: gamma and alpha steps have different trials', id);
    end
    T = Lg.T;

    % similarity per pair in the Ispec cluster
    sI = cluster_mean(Lg.simI, G.mask);
    sC = cluster_mean(Lg.simC, G.mask);

    % trial power: Ispec cluster channels x 8-25 Hz x cluster time extent
    ch = ismember(La.label, G.clusterChannels);
    fr = La.freq >= C.alpha.glmBand(1) - 1e-6 & La.freq <= C.alpha.glmBand(2) + 1e-6;
    tm = La.time >= G.timeExtent(1) - 1e-6 & La.time <= G.timeExtent(2) + 1e-6;
    if ~any(ch) || ~any(fr) || ~any(tm)
        error('g3_extract:power', 'Subject %d: empty power window', id);
    end
    pw = La.pow(:, ch, fr, tm);
    pw = mean(reshape(double(pw), size(pw,1), []), 2);      % one value per trial (row of T)

    % items (one per Isim pair)
    PI = Lg.PI; PC = Lg.PC;
    csimItem = nan(height(PI),1); nCs = zeros(height(PI),1);
    for k = 1:height(PI)
        in = PC.tA == PI.tA(k) | PC.tB == PI.tA(k);       % Csim pairs containing this item's presentation 1
        nCs(k) = sum(in);
        if nCs(k) > 0, csimItem(k) = mean(sC(in)); end
    end
    items{s} = table(repmat(id, height(PI), 1), repmat(isYA, height(PI), 1), PI.category, PI.item, PI.exemplar, ...
        PI.catSize, PI.itemRep, PI.Acc, sI, csimItem, nCs, pw(PI.tA), pw(PI.tB), PI.lag, ...
        'VariableNames', {'ID','YA','category','item','exemplar','catSize','itemRep','Acc','Isim','Csim', ...
        'nCsimPairs','Power1st','Power2nd','lagIsim'});

    % pairs (Isim and Csim stacked)
    nI = height(PI); nC = height(PC);
    pairs{s} = table(repmat(id, nI+nC, 1), repmat(isYA, nI+nC, 1), [repmat("Isim", nI, 1); repmat("Csim", nC, 1)], ...
        [sI; sC], [(pw(PI.tA) + pw(PI.tB))/2; (pw(PC.tA) + pw(PC.tB))/2], [pw(PI.tA); pw(PC.tA)], [pw(PI.tB); pw(PC.tB)], ...
        [PI.category; PC.category], [PI.catSize; PC.catSize], [PI.itemRep; T.itemRep(PC.tA)], [PI.lag; PC.lag], ...
        'VariableNames', {'ID','YA','Type','Sim','PowerPair','PowerA','PowerB','category','catSize','itemRep','lag'});

    % lags and references (for the lag figure)
    lags{s} = table(repmat(id, nI+nC, 1), repmat(isYA, nI+nC, 1), [repmat("Isim", nI, 1); repmat("Csim", nC, 1)], ...
        [PI.lag; PC.lag], [PI.posA; PC.posA], [PI.refPosDist; PC.refPosDist], [PI.refSim; PC.refSim], [PI.widen; PC.widen], ...
        'VariableNames', {'ID','YA','Type','lag','posA','refPosDist','refSim','widen'});

    % split half (categories)
    rng(C.seed.splitHalf + id, 'twister');
    cats = unique([PI.category; PC.category]);
    h1 = cats(randperm(numel(cats), floor(numel(cats)/2)));
    inH1I = ismember(PI.category, h1); inH1C = ismember(PC.category, h1);
    split.Ispec_half1(s) = mean(sI(inH1I)) - mean(sC(inH1C));
    split.Ispec_half2(s) = mean(sI(~inH1I)) - mean(sC(~inH1C));
    split.Ispec_pairs(s) = mean(sI) - mean(sC);
    fprintf('%d: %d items, %d Csim pairs, Ispec %.4f (halves %.4f / %.4f)\n', id, nI, nC, split.Ispec_pairs(s), split.Ispec_half1(s), split.Ispec_half2(s));
end

% consistency: per-pair cluster means average to the per-person value of g1
dev = max(abs(split.Ispec_pairs - G.P.Ispec));
if dev > 1e-4
    error('g3_extract:check', 'Ispec from pairs differs from g1 by up to %.2g', dev);
end

writetable(vertcat(items{:}), fullfile(C.dir.tables, 'items.csv'));
writetable(vertcat(pairs{:}), fullfile(C.dir.tables, 'pairs.csv'));
writetable(split, fullfile(C.dir.tables, 'person_split.csv'));
writetable(vertcat(lags{:}), fullfile(C.dir.tables, 'lags.csv'));
fprintf('Saved items.csv, pairs.csv, person_split.csv (check vs g1: max deviation %.2g)\n', dev);
end

function v = cluster_mean(X, mask)
% mean over all cluster points per row; X = [pair x chan x t1 x t2]
X = reshape(X, size(X,1), []);
v = mean(double(X(:, mask(:))), 2);
end

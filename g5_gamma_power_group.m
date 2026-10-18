function g5_gamma_power_group()
% G5_GAMMA_POWER_GROUP  Gamma power at the group level (analysis plan 1f, 1g;
% added 2026-10-18). Two cluster tests on the gamma power of s1 (the RSA's
% DPSS TFR, log10, no BOSC):
%
%   1f GAMMA SYNC: per person gPow = mean over all analysed trials of
%      log10 power minus the trial's own baseline (s1); group: gPow vs 0
%      (depsamplesT).
%      Cluster = all significant POSITIVE clusters (power above baseline).
%      gSync = mean gPow (log10 power - baseline) over the cluster.
%   1g GAMMA SME (as the alpha SME, g2): per person mean log10 power of
%      later 'old' and later 'similar' trials (s1; presentation 1, no
%      baseline); group: paired test old vs similar (depsamplesT).
%   Decided 2026-10-18: averages for every cluster test (s1 still saves the
%   t-maps tSync, tSME; not used).
%      Cluster = all significant POSITIVE clusters (expected: more gamma
%      for later 'old'). gSME = mean (old - similar) over the cluster.
%   In both, significant clusters in the other direction are reported, not
%   used. 35-100 Hz, 0-1.7 s. Cluster settings for both = the standard ones
%   (C.clus: clusteralpha .05, minnbchan 1), as Ispec and the alpha SME.
%
% Output:
%   C.dir.group/gsync_cluster.mat   stat, mask, sigPos, sigNeg, extent, grand averages
%   C.dir.group/gsme_cluster.mat    stat, mask, sigPos, sigNeg, extent, grand averages
%   C.dir.tables/person_gpower.csv  ID, YA, gSync, gSME, nTrials, nOld, nSim

C = bebo_config();
bebo_init_fieldtrip(C);
for d = {C.dir.group, C.dir.tables}, if ~exist(d{1}, 'dir'), mkdir(d{1}); end, end
ids = bebo_ids('all'); nYA = numel(bebo_ids('ya')); n = numel(ids);

%% Load
for s = 1:n
    L = load(fullfile(C.dir.gammaSubj, sprintf('%d_gamma.mat', ids(s))), ...
        'gPow', 'nPow', 'gOld', 'gSim', 'gnOld', 'gnSim', 'label', 'freq', 'time');
    if s == 1
        label = L.label; freq = L.freq; time = L.time;
        A = zeros([n size(L.gPow)]); Old = A; Sim = A; nPow = zeros(n,1); nOld = nPow; nSim = nPow;
    elseif ~isequal(L.label, label) || ~isequal(L.freq, freq) || ~isequal(L.time, time)
        error('g5_gamma_power_group:axes', 'Subject %d: axes differ from subject %d', ids(s), ids(1));
    end
    A(s,:,:,:) = L.gPow; Old(s,:,:,:) = L.gOld; Sim(s,:,:,:) = L.gSim;
    nPow(s) = L.nPow; nOld(s) = L.gnOld; nSim(s) = L.gnSim;
end
P = table(ids(:), (1:n)' <= nYA, nan(n,1), nan(n,1), nPow, nOld, nSim, ...
    'VariableNames', {'ID','YA','gSync','gSME','nTrials','nOld','nSim'});

%% 1f. Gamma sync: per-person mean (post - baseline) vs 0
sel = struct('frequency', C.gsync.freq, 'latency', C.gsync.latency);
fprintf('\n=== GAMMA SYNC (per-person means post - baseline, vs 0): %d persons, %d-%d Hz, %.2f-%.2f s, clusteralpha %.2f, minnbchan %d\n', ...
    n, sel.frequency, sel.latency, C.clus.clusteralpha, C.clus.minnbchan);
stat = bebo_cluster_test(A, zeros(size(A)), label, freq, time, C, C.seed.clusterG, sel);
[mask, sigPos, sigNeg, extent, P.gSync] = pos_cluster(stat, A, label, freq, time, C.clus.alpha, 'gamma sync');
GA = squeeze(mean(A, 1)); GA_YA = squeeze(mean(A(1:nYA,:,:,:), 1)); GA_OA = squeeze(mean(A(nYA+1:end,:,:,:), 1));
save(fullfile(C.dir.group, 'gsync_cluster.mat'), 'stat', 'mask', 'sigPos', 'sigNeg', 'extent', ...
    'GA', 'GA_YA', 'GA_OA', 'label', 'freq', 'time', '-v7.3');

%% 1g. Gamma SME: mean later old vs mean later similar (paired)
sel = struct('frequency', C.gsme.freq, 'latency', C.gsme.latency);
fprintf('\n=== GAMMA SME (per-person means, old vs similar): %d persons, %d-%d Hz, %.2f-%.2f s, clusteralpha %.2f, minnbchan %d\n', ...
    n, sel.frequency, sel.latency, C.clus.clusteralpha, C.clus.minnbchan);
stat = bebo_cluster_test(Old, Sim, label, freq, time, C, C.seed.clusterGS, sel);
D = Old - Sim;
[mask, sigPos, sigNeg, extent, P.gSME] = pos_cluster(stat, D, label, freq, time, C.clus.alpha, 'gamma SME');
GA = squeeze(mean(D, 1)); GA_YA = squeeze(mean(D(1:nYA,:,:,:), 1)); GA_OA = squeeze(mean(D(nYA+1:end,:,:,:), 1));
save(fullfile(C.dir.group, 'gsme_cluster.mat'), 'stat', 'mask', 'sigPos', 'sigNeg', 'extent', ...
    'GA', 'GA_YA', 'GA_OA', 'label', 'freq', 'time', '-v7.3');

writetable(P, fullfile(C.dir.tables, 'person_gpower.csv'));
fprintf('Saved gsync_cluster.mat, gsme_cluster.mat and person_gpower.csv\n');
end

function [mask, sigPos, sigNeg, extent, val] = pos_cluster(stat, X, label, freq, time, alpha, name)
% mask of all significant POSITIVE clusters; per-person mean of X over it
n = size(X, 1);
mask = false(size(stat.stat)); sigPos = []; sigNeg = [];
if isfield(stat, 'posclusters') && ~isempty(stat.posclusters)
    sigPos = find([stat.posclusters.prob] < alpha);
    for k = sigPos, mask = mask | stat.posclusterslabelmat == k; end
end
if isfield(stat, 'negclusters') && ~isempty(stat.negclusters)
    sigNeg = find([stat.negclusters.prob] < alpha);
end
fprintf('Significant positive clusters: %d; negative: %d (reported, not used)\n', numel(sigPos), numel(sigNeg));
for k = sigPos, fprintf('  pos cluster %d: p = %.4f, sum t = %.1f\n', k, stat.posclusters(k).prob, stat.posclusters(k).clusterstat); end
for k = sigNeg, fprintf('  neg cluster %d: p = %.4f, sum t = %.1f\n', k, stat.negclusters(k).prob, stat.negclusters(k).clusterstat); end
val = nan(n, 1); extent = struct();
if ~any(mask(:))
    warning('g5_gamma_power_group:noCluster', 'No significant positive %s cluster -- per-person value is NaN.', name);
    return
end
[okC, iC] = ismember(stat.label, label);
iF = arrayfun(@(f) find(abs(freq - f) < 1e-6, 1), stat.freq);
iT = arrayfun(@(t) find(abs(time - t) < 1e-6, 1), stat.time);
if ~all(okC), error('g5_gamma_power_group:chan', 'stat channels missing in the data'); end
D = reshape(X(:, iC, iF, iT), n, []);
val = mean(D(:, mask(:)), 2);
extent.channels = stat.label(any(reshape(mask, size(mask,1), []), 2));
extent.freq = stat.freq(squeeze(any(any(mask, 1), 3)));
extent.time = stat.time(squeeze(any(any(mask, 1), 2)));
fprintf('%s cluster: %d channels, %g-%g Hz, %.2f-%.2f s\n', name, numel(extent.channels), ...
    min(extent.freq), max(extent.freq), min(extent.time), max(extent.time));
end

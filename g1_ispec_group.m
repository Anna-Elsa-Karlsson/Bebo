function g1_ispec_group()
% G1_ISPEC_GROUP  Ispec cluster test and per-person Ispec (analysis plan 1a, 1c).
%
%   1. per-person t-maps (s1) tested against 0 across the full sample
%      (depsamplesT, cluster permutation, settings in bebo_config)
%   2. Ispec cluster = all significant POSITIVE clusters
%   3. per person, mean over all cluster points of the corrected Isim and
%      Csim; Ispec = Isim - Csim (MEAN DIFFERENCES, not t). Also the drift
%      check: baseline(Isim) - baseline(Csim) in the cluster, and the
%      uncorrected values.
%
% Output:
%   C.dir.group/ispec_cluster.mat   stat, mask [chan x t1 x t2], clusterChannels,
%                                   timeExtent [s], P (person table)
%   C.dir.group/ispec_group_means.mat   GA.YA/OA.Isim, .Csim (corrected) and
%                                   .rawIsim, .rawCsim: mean over persons of the
%                                   per-person mean maps [chan x t1 x t2]
%                                   (= the old grand averages; for Fig3)
%   C.dir.tables/person_ispec.csv   ID, YA, Isim, Csim, Ispec, rawIsim, rawCsim,
%                                   refIsim, refCsim (baseline), drift

C = bebo_config();
bebo_init_fieldtrip(C);
for d = {C.dir.group, C.dir.tables}, if ~exist(d{1}, 'dir'), mkdir(d{1}); end, end
ids = bebo_ids('all'); nYA = numel(bebo_ids('ya')); n = numel(ids);

%% 1. Load the per-person maps
vars = {'tmap','meanI','meanC','rawMeanI','rawMeanC','refMeanI','refMeanC','label','time'};
for s = 1:n
    L = load(fullfile(C.dir.gammaSubj, sprintf('%d_gamma.mat', ids(s))), vars{:});
    if s == 1
        label = L.label; time = L.time;
        sz = size(L.tmap);
        M = struct(); for v = vars(1:7), M.(v{1}) = zeros([n sz], 'single'); end
    elseif ~isequal(L.label, label) || ~isequal(L.time, time)
        error('g1_ispec_group:axes', 'Subject %d: channels or time axis differ from subject %d', ids(s), ids(1));
    end
    for v = vars(1:7), M.(v{1})(s,:,:,:) = single(L.(v{1})); end
end

%% 2. Cluster test: t-maps vs 0
fprintf('Ispec cluster test: %d persons, %d channels x %d x %d time points\n', n, sz(1), sz(2), sz(3));
stat = bebo_cluster_test(M.tmap, zeros(size(M.tmap), 'single'), label, time, time, C, C.seed.clusterI);

mask = false(size(stat.stat));
sigPos = [];
if isfield(stat, 'posclusters') && ~isempty(stat.posclusters)
    sigPos = find([stat.posclusters.prob] < C.clus.alpha);
    for k = sigPos, mask = mask | stat.posclusterslabelmat == k; end
end
sigNeg = [];
if isfield(stat, 'negclusters') && ~isempty(stat.negclusters)
    sigNeg = find([stat.negclusters.prob] < C.clus.alpha);
end
fprintf('Significant positive clusters: %d (p = %s); negative: %d (not used)\n', numel(sigPos), ...
    mat2str(round([stat.posclusters(sigPos).prob], 4)), numel(sigNeg));
if ~any(mask(:))
    save(fullfile(C.dir.group, 'ispec_cluster.mat'), 'stat', 'mask', '-v7.3');
    error('g1_ispec_group:noCluster', 'No significant positive Ispec cluster -- the steps that need the Ispec cluster cannot run.');
end
clusterChannels = label(any(reshape(mask, size(mask,1), []), 2));
[~, i2, i3] = ind2sub(size(mask), find(mask));
timeExtent = [min([time(i2) time(i3)]) max([time(i2) time(i3)])];
fprintf('Ispec cluster: %d channels, time extent %.2f-%.2f s, %d points\n', numel(clusterChannels), timeExtent, nnz(mask));

%% 3. Per-person values in the cluster (mean differences)
inCl = @(X) cluster_mean(X, mask);
P = table(ids(:), (1:n)' <= nYA, 'VariableNames', {'ID','YA'});
P.Isim    = inCl(M.meanI);    P.Csim    = inCl(M.meanC);    P.Ispec = P.Isim - P.Csim;
P.rawIsim = inCl(M.rawMeanI); P.rawCsim = inCl(M.rawMeanC);
P.refIsim = inCl(M.refMeanI); P.refCsim = inCl(M.refMeanC); P.drift = P.refIsim - P.refCsim;

%% 4. Group-average similarity matrices (= the old grand-average files,
%  for the Fig3 similarity-matrix figure): mean over persons per age group
%  of the per-person mean corrected Isim and Csim maps (and uncorrected)
GA = struct();
grp = struct('YA', (1:n)' <= nYA, 'OA', (1:n)' > nYA);
for g = {'YA', 'OA'}
    ix = grp.(g{1});
    GA.(g{1}).Isim    = squeeze(mean(M.meanI(ix,:,:,:), 1));
    GA.(g{1}).Csim    = squeeze(mean(M.meanC(ix,:,:,:), 1));
    GA.(g{1}).rawIsim = squeeze(mean(M.rawMeanI(ix,:,:,:), 1));
    GA.(g{1}).rawCsim = squeeze(mean(M.rawMeanC(ix,:,:,:), 1));
end
save(fullfile(C.dir.group, 'ispec_group_means.mat'), 'GA', 'label', 'time', '-v7.3');

save(fullfile(C.dir.group, 'ispec_cluster.mat'), 'stat', 'mask', 'clusterChannels', 'timeExtent', ...
    'sigPos', 'sigNeg', 'P', 'label', 'time', '-v7.3');
writetable(P, fullfile(C.dir.tables, 'person_ispec.csv'));
fprintf('Saved ispec_cluster.mat and person_ispec.csv\n');
end

function v = cluster_mean(X, mask)
% per-person mean over all cluster points; X = [subject x chan x t1 x t2]
X = reshape(double(X), size(X,1), []);
v = mean(X(:, mask(:)), 2);
end

function g2_sme_group()
% G2_SME_GROUP  Alpha SME cluster test and per-person SME (analysis plan 1b, 1d).
%
%   1. per person: 1st-level t-map, later 'old' vs later 'similar' trials
%      (s2: indepsamplesT, log10 power, no baseline; as the original SME)
%   2. group: the t-maps vs 0 (depsamplesT), full sample, 2-30 Hz, -0.3 to
%      1.7 s, cluster permutation (bebo_config)
%   3. SME cluster = all significant NEGATIVE clusters (expected direction:
%      more desynchronisation for later 'old'); a significant positive
%      cluster is reported, not used
%   4. per person: mean over the SME cluster of (old - similar) = SME
%      (MEAN DIFFERENCE, not t)
%
% Output:
%   C.dir.group/sme_cluster.mat   stat, mask [chan x freq x time on the stat axes],
%                                 P (person table), cluster extent
%   C.dir.tables/person_sme.csv   ID, YA, SME, nOld, nSim

C = bebo_config();
bebo_init_fieldtrip(C);
for d = {C.dir.group, C.dir.tables}, if ~exist(d{1}, 'dir'), mkdir(d{1}); end, end
ids = bebo_ids('all'); nYA = numel(bebo_ids('ya')); n = numel(ids);

%% 1. Load
for s = 1:n
    L = load(fullfile(C.dir.alphaSubj, sprintf('%d_alpha.mat', ids(s))), 'tmap', 'powOld', 'powSim', 'nOld', 'nSim', 'label', 'freq', 'time');
    if s == 1
        label = L.label; freq = L.freq; time = L.time;
        A = zeros([n size(L.powOld)]); B = A; TM = A; nOld = zeros(n,1); nSim = nOld;
    elseif ~isequal(L.label, label) || ~isequal(L.freq, freq) || ~isequal(L.time, time)
        error('g2_sme_group:axes', 'Subject %d: axes differ from subject %d', ids(s), ids(1));
    end
    A(s,:,:,:) = L.powOld; B(s,:,:,:) = L.powSim; TM(s,:,:,:) = L.tmap; nOld(s) = L.nOld; nSim(s) = L.nSim;
end

%% 2. Cluster test: 1st-level t-maps (later old vs later similar) vs 0
sel = struct('frequency', C.clus.smeFreq, 'latency', C.clus.smeLatency);
fprintf('SME cluster test (t-maps vs 0): %d persons, %d-%d Hz, %.2f-%.2f s\n', n, sel.frequency, sel.latency);
stat = bebo_cluster_test(TM, zeros(size(TM)), label, freq, time, C, C.seed.clusterS, sel);

%% 3. SME cluster = significant negative clusters
mask = false(size(stat.stat));
sigNeg = []; sigPos = [];
if isfield(stat, 'negclusters') && ~isempty(stat.negclusters)
    sigNeg = find([stat.negclusters.prob] < C.clus.alpha);
    for k = sigNeg, mask = mask | stat.negclusterslabelmat == k; end
end
if isfield(stat, 'posclusters') && ~isempty(stat.posclusters)
    sigPos = find([stat.posclusters.prob] < C.clus.alpha);
end
fprintf('Significant negative clusters: %d; positive: %d (reported, not used)\n', numel(sigNeg), numel(sigPos));

P = table(ids(:), (1:n)' <= nYA, nan(n,1), nOld, nSim, 'VariableNames', {'ID','YA','SME','nOld','nSim'});
extent = struct();
if any(mask(:))
    % map the stat axes back onto the full arrays
    [okC, iC] = ismember(stat.label, label);
    iF = arrayfun(@(f) find(abs(freq - f) < 1e-6, 1), stat.freq);
    iT = arrayfun(@(t) find(abs(time - t) < 1e-6, 1), stat.time);
    if ~all(okC), error('g2_sme_group:chan', 'stat channels missing in the data'); end
    D = A(:, iC, iF, iT) - B(:, iC, iF, iT);
    D = reshape(D, n, []);
    P.SME = mean(D(:, mask(:)), 2);
    extent.channels = stat.label(any(reshape(mask, size(mask,1), []), 2));
    extent.freq = stat.freq(squeeze(any(any(mask, 1), 3)));
    extent.time = stat.time(squeeze(any(any(mask, 1), 2)));
    fprintf('SME cluster: %d channels, %g-%g Hz, %.2f-%.2f s\n', numel(extent.channels), ...
        min(extent.freq), max(extent.freq), min(extent.time), max(extent.time));
else
    warning('g2_sme_group:noCluster', 'No significant negative SME cluster -- SME per person is NaN.');
end

save(fullfile(C.dir.group, 'sme_cluster.mat'), 'stat', 'mask', 'sigNeg', 'sigPos', 'extent', 'P', '-v7.3');
writetable(P, fullfile(C.dir.tables, 'person_sme.csv'));
fprintf('Saved sme_cluster.mat and person_sme.csv\n');
end

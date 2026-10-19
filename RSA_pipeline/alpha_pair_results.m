function alpha_pair_results(pdat)
% RESULTS step (2026-10-04; control block 2026-10-05): does alpha desync
% shape the representation? Within each person, across similarity pairs
% (T_pair, built by aggregate_singleTrial_merged.m), per pair type the OLS
% slope of similarity on alpha (z-scored within person; LOWER alpha =
% stronger desync, so a NEGATIVE slope = more similarity with stronger
% desync). Within-person, so skull / conductance / absolute power / age
% differences in power cannot drive it. Then the per-person slopes and
% their differences across people: one-sample t vs 0 (all, YA, OA) and YA
% vs OA (t-test, equal variances, Cohen d).
%
% Two blocks:
%  MAIN     Isim vs Csim; alpha = mean of the pair's two trials (PowerZ)
%           -> slopes + difference Isim - Csim (item specificity)
%  CONTROL  Isim, Csim and Bsim (between-category: exemplar 1 of A vs
%           exemplar 2 of every other category, averaged); alpha = the
%           1ST pair member for all three (PowerRefZ: presentation 1 /
%           exemplar 1 of A), so the types are on the same footing ->
%           slopes + differences Isim - Bsim, Csim - Bsim, Isim - Csim.
%           If desync made ALL patterns more reproducible (a general state
%           effect), Bsim would follow and the differences would be ~0.
%           Only run if Bsim rows exist.
% People need >= minPairs rows of each type in a block. Also saved: group-
% mean similarity per alpha quintile and type (figures). Mixed model:
% Scripts/BEBO_pair_lmer.R.
%
%   pdat : .stm (Results/SingleTrialMerged/<variant>), .out (stats dir), .variant

minPairs = 5;
L = load(fullfile(pdat.stm, [pdat.variant '_singleTrial_pairs_fullSample.mat']), 'T_pair');
T = L.T_pair;

ID_ya = ids_ALL('ya'); ID_oa = ids_ALL('oa');
ids = [ID_ya(:); ID_oa(:)];
nYA = numel(ID_ya); nOA = numel(ID_oa);
idx_YA = (1:nYA)'; idx_OA = (nYA+1:nYA+nOA)';

blocks = struct();
blocks.main = run_block(T, ids, idx_YA, idx_OA, {'Isim','Csim'}, 'PowerZ', {'Isim','Csim'}, minPairs);
if any(strcmp(T.Type, 'Bsim'))
    blocks.control = run_block(T, ids, idx_YA, idx_OA, {'Isim','Csim','Bsim'}, 'PowerRefZ', ...
        {'Isim','Bsim'; 'Csim','Bsim'; 'Isim','Csim'}, minPairs);
end

if ~exist(pdat.out, 'dir'), mkdir(pdat.out); end
save(fullfile(pdat.out, 'alpha_pair_results.mat'), 'ids', 'idx_YA', 'idx_OA', 'nYA', 'nOA', 'blocks', 'minPairs');

fid = fopen(fullfile(pdat.out, 'alpha_pair_stats.txt'), 'w');
fprintf(fid, 'ALPHA x PAIR TYPE: does alpha desync shape the representation? -- %s\n', datestr(now));
fprintf(fid, 'Within person, per pair type: slope of similarity on alpha (z-scored within person; LOWER alpha = stronger\n');
fprintf(fid, 'desync, so a NEGATIVE slope = more similarity with stronger desync). Then across people.\n\n');
print_block(fid, blocks.main, 'MAIN: Isim vs Csim, alpha = mean of the pair''s two trials', minPairs, idx_YA, idx_OA);
if isfield(blocks, 'control')
    print_block(fid, blocks.control, ['CONTROL: Isim, Csim, Bsim (between-category, all other categories averaged), ' ...
        'alpha = 1st pair member (presentation 1 / exemplar 1)'], minPairs, idx_YA, idx_OA);
else
    fprintf(fid, 'CONTROL block not run: no between-category (Bsim) rows in T_pair.\n');
end
fclose(fid);
fprintf('Saved %s\n', fullfile(pdat.out, 'alpha_pair_stats.txt'));
end


function B = run_block(T, ids, idx_YA, idx_OA, types, powVar, dpairs, minPairs)
nBin = 5; nAll = numel(ids); nT = numel(types);
slope = nan(nAll, nT); nPairs = zeros(nAll, nT);
binSim = nan(nAll, nBin, nT); binPow = nan(nAll, nBin, nT);
for s = 1:nAll
    r = T.ID == ids(s) & ~isnan(T.(powVar)) & ~isnan(T.Sim);
    for ty = 1:nT
        rr = find(r & strcmp(T.Type, types{ty}));
        nPairs(s,ty) = numel(rr);
        if numel(rr) < minPairs, continue; end
        b = [ones(numel(rr),1), T.(powVar)(rr)] \ T.Sim(rr);
        slope(s,ty) = b(2);
        [~, ord] = sort(T.(powVar)(rr));
        edges = round(linspace(0, numel(rr), nBin+1));
        for k = 1:nBin
            ii = rr(ord(edges(k)+1:edges(k+1)));
            binSim(s,k,ty) = mean(T.Sim(ii)); binPow(s,k,ty) = mean(T.(powVar)(ii));
        end
    end
end
included = all(~isnan(slope), 2);
slope(~included,:) = NaN;   % the same people for every comparison in the block

meas = struct('name', {}, 'label', {}, 'v', {});
for ty = 1:nT
    meas(end+1) = struct('name', ['slope_' types{ty}], 'label', sprintf('slope alpha -> %s', types{ty}), 'v', slope(:,ty)); %#ok<AGROW>
end
if ~isempty(dpairs) && size(dpairs,2) == 2
    for k = 1:size(dpairs,1)
        a = strcmp(types, dpairs{k,1}); b = strcmp(types, dpairs{k,2});
        meas(end+1) = struct('name', sprintf('diff_%s_%s', dpairs{k,1}, dpairs{k,2}), ...
            'label', sprintf('slope %s - slope %s', dpairs{k,1}, dpairs{k,2}), 'v', slope(:,a) - slope(:,b)); %#ok<AGROW>
    end
end
G = struct();
for m = 1:numel(meas)
    v = meas(m).v;
    G.(meas(m).name) = struct('all', t0(v), 'YA', t0(v(idx_YA)), 'OA', t0(v(idx_OA)), 'age', t2(v(idx_YA), v(idx_OA)));
end
grpBin = struct();
for g = {'YA','OA'}
    ix = idx_YA; if strcmp(g{1}, 'OA'), ix = idx_OA; end
    ix = ix(included(ix));
    grpBin.(g{1}).sim_M  = reshape(mean(binSim(ix,:,:), 1, 'omitnan'), nBin, nT);
    grpBin.(g{1}).sim_SE = reshape(std(binSim(ix,:,:), 0, 1, 'omitnan'), nBin, nT) / sqrt(numel(ix));
    grpBin.(g{1}).pow_M  = reshape(mean(binPow(ix,:,:), 1, 'omitnan'), nBin, nT);
end
B = struct('types', {types}, 'powVar', powVar, 'slope', slope, 'nPairs', nPairs, 'included', included, ...
    'meas', meas, 'G', G, 'grpBin', grpBin);
end

function print_block(fid, B, title_str, minPairs, idx_YA, idx_OA)
fprintf(fid, '==================== %s ====================\n', title_str);
fprintf(fid, 'Included: %d people (>= %d rows of each type): YA %d, OA %d. Median rows per person: %s\n\n', ...
    sum(B.included), minPairs, sum(B.included(idx_YA)), sum(B.included(idx_OA)), ...
    strjoin(arrayfun(@(k) sprintf('%s %.0f', B.types{k}, median(B.nPairs(B.included,k))), 1:numel(B.types), 'uni', 0), ', '));
for m = 1:numel(B.meas)
    g = B.G.(B.meas(m).name);
    fprintf(fid, '--- %s ---\n', B.meas(m).label);
    pr0(fid, 'All', g.all); pr0(fid, 'YA ', g.YA); pr0(fid, 'OA ', g.OA);
    fprintf(fid, 'YA vs OA: t(%d)=%.3f, p=%.4f, d=%.3f\n\n', g.age.df, g.age.t, g.age.p, g.age.d);
end
end

function R = t0(v)
v = v(~isnan(v));
[~, p, ~, st] = ttest(v);
R = struct('M', mean(v), 'SD', std(v), 'n', numel(v), 't', st.tstat, 'df', st.df, 'p', p, 'dz', mean(v)/std(v));
end

function R = t2(a, b)
a = a(~isnan(a)); b = b(~isnan(b));
[~, p, ~, st] = ttest2(a, b);
d = (mean(a) - mean(b)) / sqrt(((numel(a)-1)*var(a) + (numel(b)-1)*var(b)) / (numel(a)+numel(b)-2));
R = struct('t', st.tstat, 'df', st.df, 'p', p, 'd', d);
end

function pr0(fid, lab, g)
fprintf(fid, '%s: M=%.4f SD=%.4f (n=%d), t(%d)=%.3f, p=%.4f, dz=%.3f\n', lab, g.M, g.SD, g.n, g.df, g.t, g.p, g.dz);
end

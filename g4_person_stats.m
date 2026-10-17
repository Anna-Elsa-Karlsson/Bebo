function g4_person_stats()
% G4_PERSON_STATS  Between-person statistics (analysis plan 1c, 1d, 1e, 3).
% Reads the person tables of g1-g3 and the behavioural file; writes
% C.dir.stats/person_stats.txt (+ person_stats.mat for the figure script)
% and C.dir.tables/person_all.csv.
%
%   Age (YA vs OA, Student t, Cohen's d): Ispec, Isim, Csim, SME
%   Drift check: reference(Isim) - reference(Csim) in the Ispec cluster;
%     20% of its mean = the largest Ispec the remaining drift could produce
%   Position check: per person, r between session position and corrected
%     similarity (Isim, Csim separately), tested against 0 across persons
%   Correlations: partial r controlling age (full sample), r within YA and
%     OA, Fisher z YA vs OA. BH-FDR within each family, separately for the
%     partial r and the Fisher z:
%       family 1  Ispec x IM, Ispec x LDI
%       family 2  SME x IM,   SME x LDI
%       single    Ispec x SME
%   Split half: the age test in each half; split-half reliability
%     (Pearson r between halves, Spearman-Brown corrected)

C = bebo_config();
if ~exist(C.dir.stats, 'dir'), mkdir(C.dir.stats); end
ids = bebo_ids('all'); n = numel(ids);

I = readtable(fullfile(C.dir.tables, 'person_ispec.csv'));
S = readtable(fullfile(C.dir.tables, 'person_sme.csv'));
H = readtable(fullfile(C.dir.tables, 'person_split.csv'));
B = readtable(C.perfFile);
if ~isequal(I.ID, ids(:)) || ~isequal(S.ID, ids(:)) || ~isequal(H.ID, ids(:))
    error('g4_person_stats:ids', 'Person tables are not in bebo_ids order');
end
[ok, ib] = ismember(ids(:), B.ID);
if ~all(ok), error('g4_person_stats:perf', 'Performance file is missing: %s', mat2str(ids(~ok))); end
P = I; P.SME = S.SME; P.nOld = S.nOld; P.nSim = S.nSim;
P.Ispec_half1 = H.Ispec_half1; P.Ispec_half2 = H.Ispec_half2;
P.IM = B.IM(ib); P.LDI = B.LDI(ib);
writetable(P, fullfile(C.dir.tables, 'person_all.csv'));

ya = P.YA == 1; age = double(~ya) + 1;   % 1 = YA, 2 = OA (covariate for partialcorr)
fid = fopen(fullfile(C.dir.stats, 'person_stats.txt'), 'w');
out = @(varargin) fprintf(fid, varargin{:});
out('BEBO person-level statistics -- %s\nN = %d (%d YA, %d OA). Per-person values = mean differences in the cluster.\n\n', ...
    datestr(now), n, sum(ya), sum(~ya));

%% Age differences
R = struct();
out('=== AGE DIFFERENCES (YA vs OA; Student t, Cohen''s d) ===\n');
for v = {'Ispec', 'Isim', 'Csim', 'SME'}
    x = P.(v{1}); R.age.(v{1}) = age_test(out, v{1}, x(ya), x(~ya));
end

%% Drift check
out('\n=== DRIFT CHECK (Ispec cluster) ===\n');
out('Drift effect = reference(Isim) - reference(Csim): M = %.4f (SD %.4f); uncorrected Isim - Csim: M = %.4f\n', ...
    mean(P.drift), std(P.drift), mean(P.rawIsim - P.rawCsim));
out('Largest Ispec the remaining drift could produce (20%% of the drift effect): %.4f | observed Ispec M = %.4f\n', ...
    0.2 * mean(P.drift), mean(P.Ispec));

%% Position check (lag correction): do the corrected similarities still
%  depend on where in the session a pair lies? Per person and pair type:
%  Pearson r between pair position (midpoint) and corrected similarity;
%  Fisher-z r tested against 0 across persons (one-sample t).
out('\n=== POSITION CHECK (lag correction; corrected similarity vs session position) ===\n');
PR = readtable(fullfile(C.dir.tables, 'pairs.csv'));
for ty = {'Isim', 'Csim'}
    zr = nan(n, 1);
    for s = 1:n
        ix = PR.ID == ids(s) & strcmp(PR.Type, ty{1});
        if sum(ix) >= 4, zr(s) = atanh(corr(PR.pos(ix), PR.Sim(ix))); end
    end
    [~, p, ~, st] = ttest(zr);
    out('%s: mean r(position, corrected similarity) = %.3f; t(%d) = %.2f, p = %.4f\n', ...
        ty{1}, tanh(mean(zr, 'omitnan')), st.df, st.tstat, p);
    R.position.(ty{1}) = struct('r', tanh(mean(zr, 'omitnan')), 't', st.tstat, 'df', st.df, 'p', p);
end

%% Correlations
out('\n=== CORRELATIONS (partial r controlling age; BH-FDR within family) ===\n');
R.Ispec_beh = corr_family(out, 'Family 1: Ispec', P.Ispec, {P.IM, P.LDI}, {'IM', 'LDI'}, age, ya);
R.SME_beh   = corr_family(out, 'Family 2: SME',   P.SME,   {P.IM, P.LDI}, {'IM', 'LDI'}, age, ya);
R.Ispec_SME = corr_family(out, 'Single test',     P.Ispec, {P.SME},       {'SME'},       age, ya);

%% Split half
out('\n=== SPLIT-HALF CHECK (categories split at random per person) ===\n');
age_test(out, 'Ispec, half 1', P.Ispec_half1(ya), P.Ispec_half1(~ya));
age_test(out, 'Ispec, half 2', P.Ispec_half2(ya), P.Ispec_half2(~ya));
r = corr(P.Ispec_half1, P.Ispec_half2);
out('Split-half reliability: r = %.3f, Spearman-Brown = %.3f\n', r, 2*r/(1+r));
fclose(fid);
save(fullfile(C.dir.stats, 'person_stats.mat'), 'R', 'P');   % numbers for the figure script
fprintf('%s', fileread(fullfile(C.dir.stats, 'person_stats.txt')));
end

function res = age_test(out, name, a, b)
[~, p, ~, st] = ttest2(a, b);
sp = sqrt(((numel(a)-1)*var(a) + (numel(b)-1)*var(b)) / (numel(a)+numel(b)-2));
out('%-14s YA M=%.4f SD=%.4f | OA M=%.4f SD=%.4f | t(%d)=%.3f, p=%.4f, d=%.3f\n', ...
    name, mean(a), std(a), mean(b), std(b), st.df, st.tstat, p, (mean(a)-mean(b))/sp);
res = struct('t', st.tstat, 'df', st.df, 'p', p, 'd', (mean(a)-mean(b))/sp);
end

function res = corr_family(out, title, x, ys, names, age, ya)
k = numel(ys); r = nan(1,k); p = r; rY = r; rO = r; z = r; pz = r;
nY = sum(ya); nO = sum(~ya);
for j = 1:k
    y = ys{j};
    [r(j), p(j)] = partialcorr(x, y, age);
    rY(j) = corr(x(ya), y(ya)); rO(j) = corr(x(~ya), y(~ya));
    z(j) = (atanh(rY(j)) - atanh(rO(j))) / sqrt(1/(nY-3) + 1/(nO-3));
    pz(j) = 2 * (1 - normcdf(abs(z(j))));
end
q = bh(p); qz = bh(pz);
out('%s (%d test%s)\n', title, k, repmat('s', 1, k > 1));
for j = 1:k
    out('  x %-4s partial r=%.3f, p=%.4f, p_FDR=%.4f | YA r=%.3f, OA r=%.3f | Fisher z=%.2f, p=%.4f, p_FDR=%.4f\n', ...
        names{j}, r(j), p(j), q(j), rY(j), rO(j), z(j), pz(j), qz(j));
end
res = struct('names', {names}, 'r', r, 'p', p, 'pFDR', q, 'rYA', rY, 'rOA', rO, 'z', z, 'pz', pz, 'pzFDR', qz);
end

function q = bh(p)
% Benjamini-Hochberg adjusted p-values
[ps, o] = sort(p); m = numel(p);
qs = min(1, cummin(ps .* m ./ (1:m), 'reverse'));
q = zeros(size(p)); q(o) = qs;
end

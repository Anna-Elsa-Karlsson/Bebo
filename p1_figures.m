function p1_figures()
% P1_FIGURES  Group-level figures (PLOT step only; no statistics computed).
% Reads the outputs of g1-g4 and draws:
%   Fig3_SimMatrices  group-average Isim and Csim matrices, YA and OA (old Fig3)
%   Fig_Ispec   t-map (mean over cluster channels) with cluster outline,
%               cluster topography, Isim/Csim and Ispec per age group,
%               split halves
%   Fig_SME     t-map (mean over SME cluster channels) with outline,
%               topography, SME per age group
%   Fig_corr    between-person correlations (numbers from g4)
%   Fig_lags    lags of Isim / Csim pairs, reference similarity by lag,
%               drift effect per person, position distance of references
% R2017b-compatible (subplot, print). Output: C.dir.figs (png + pdf).

C = bebo_config();
if ~exist(C.dir.figs, 'dir'), mkdir(C.dir.figs); end
haveFT = exist(C.layoutFile, 'file') == 2;
if haveFT, bebo_init_fieldtrip(C); end
col = struct('YA', [30 144 255]/255, 'OA', [80 80 80]/255);

GI = load(fullfile(C.dir.group, 'ispec_cluster.mat'));
GS = load(fullfile(C.dir.group, 'sme_cluster.mat'));
S  = load(fullfile(C.dir.stats, 'person_stats.mat'));   % R (numbers), P (person table)
P  = S.P; R = S.R; ya = P.YA == 1;
lagT = readtable(fullfile(C.dir.tables, 'lags.csv'));

%% Fig_Ispec
f = figure('Units', 'centimeters', 'Position', [2 2 34 20], 'Color', 'w');
chI = any(reshape(GI.mask, size(GI.mask,1), []), 2);
ax = subplot(2,3,1);
tm = squeeze(mean(GI.stat.stat(chI,:,:), 1));
imagesc(ax, GI.time, GI.time, tm'); axis(ax, 'xy', 'square'); hold(ax, 'on');
contour(ax, GI.time, GI.time, double(squeeze(any(GI.mask(chI,:,:), 1)))', [0.5 0.5], 'k', 'LineWidth', 1.2);
colormap(ax, bwr(256)); c = max(abs(tm(:))); caxis(ax, [-c c]); colorbar(ax);
xlabel(ax, 'Presentation 1 / exemplar A time (s)'); ylabel(ax, 'Presentation 2 / exemplar B time (s)');
title(ax, sprintf('Ispec t (mean over %d cluster channels)', sum(chI)));
ax = subplot(2,3,2);
topo(ax, GI.label, sum(reshape(GI.mask, size(GI.mask,1), []), 2), C, haveFT, 'cluster points per channel');
ax = subplot(2,3,3);
group_dots(ax, {P.Isim(ya), P.Isim(~ya), P.Csim(ya), P.Csim(~ya)}, [1 2 4 5], {col.YA, col.OA, col.YA, col.OA});
set(ax, 'XTick', [1.5 4.5], 'XTickLabel', {'Item similarity', 'Category similarity'});
ylabel(ax, 'Corrected similarity (Fisher z)'); title(ax, 'Isim and Csim (blue YA, grey OA)');
ax = subplot(2,3,4);
group_dots(ax, {P.Ispec(ya), P.Ispec(~ya)}, [1 2], {col.YA, col.OA});
set(ax, 'XTick', [1 2], 'XTickLabel', {'Younger', 'Older'}); ylabel(ax, 'Ispec (Isim - Csim)');
a = R.age.Ispec; title(ax, sprintf('t(%d) = %.2f, p = %.3f, d = %.2f', a.df, a.t, a.p, a.d));
ax = subplot(2,3,5);
scatter(ax, P.Ispec_half1(ya), P.Ispec_half2(ya), 25, col.YA, 'filled'); hold(ax, 'on');
scatter(ax, P.Ispec_half1(~ya), P.Ispec_half2(~ya), 25, col.OA, 'filled');
xlabel(ax, 'Ispec, half 1'); ylabel(ax, 'Ispec, half 2'); title(ax, 'Split half (random categories)');
save_fig(f, fullfile(C.dir.figs, 'Fig_Ispec'));

%% Fig3_SimMatrices (group-average similarity, as the old Fig3)
GM = load(fullfile(C.dir.group, 'ispec_group_means.mat'));
f = figure('Units', 'centimeters', 'Position', [2 2 26 22], 'Color', 'w');
outline = double(squeeze(any(GI.mask(chI,:,:), 1)))';
panels = {'YA','Isim','Younger - item similarity'; 'YA','Csim','Younger - category similarity'; ...
          'OA','Isim','Older - item similarity';   'OA','Csim','Older - category similarity'};
mats = cell(4,1);
for k = 1:4, mats{k} = squeeze(mean(GM.GA.(panels{k,1}).(panels{k,2})(chI,:,:), 1))'; end
lim = [min(cellfun(@(m) min(m(:)), mats)) max(cellfun(@(m) max(m(:)), mats))];
for k = 1:4
    ax = subplot(2,2,k);
    imagesc(ax, GM.time, GM.time, mats{k}); axis(ax, 'xy', 'square'); hold(ax, 'on');
    contour(ax, GM.time, GM.time, outline, [0.5 0.5], 'k', 'LineWidth', 1);
    colormap(ax, flipud(autumn(256))); caxis(ax, lim); colorbar(ax);
    xlabel(ax, 'Presentation 1 / exemplar A time (s)'); ylabel(ax, 'Presentation 2 / exemplar B time (s)');
    title(ax, panels{k,3});
end
annotation(f, 'textbox', [0 0.95 1 0.05], 'String', sprintf('Corrected similarity (Fisher z), mean over the %d Ispec cluster channels; outline = Ispec cluster', sum(chI)), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
save_fig(f, fullfile(C.dir.figs, 'Fig3_SimMatrices'));

%% Fig_SME
f = figure('Units', 'centimeters', 'Position', [2 2 34 10], 'Color', 'w');
if any(GS.mask(:))
    chS = any(reshape(GS.mask, size(GS.mask,1), []), 2);
    ax = subplot(1,3,1);
    tm = squeeze(mean(GS.stat.stat(chS,:,:), 1));
    imagesc(ax, GS.stat.time, GS.stat.freq, tm); axis(ax, 'xy'); hold(ax, 'on');
    contour(ax, GS.stat.time, GS.stat.freq, double(squeeze(any(GS.mask(chS,:,:), 1))), [0.5 0.5], 'k', 'LineWidth', 1.2);
    colormap(ax, bwr(256)); c = max(abs(tm(:))); caxis(ax, [-c c]); colorbar(ax);
    xlabel(ax, 'Time (s)'); ylabel(ax, 'Frequency (Hz)');
    title(ax, sprintf('SME t (later old - similar), %d channels', sum(chS)));
    ax = subplot(1,3,2);
    topo(ax, GS.stat.label, sum(reshape(GS.mask, size(GS.mask,1), []), 2), C, haveFT, 'cluster points per channel');
    ax = subplot(1,3,3);
    group_dots(ax, {P.SME(ya), P.SME(~ya)}, [1 2], {col.YA, col.OA});
    set(ax, 'XTick', [1 2], 'XTickLabel', {'Younger', 'Older'}); ylabel(ax, 'SME (log10 power, old - similar)');
    a = R.age.SME; title(ax, sprintf('t(%d) = %.2f, p = %.3f, d = %.2f', a.df, a.t, a.p, a.d));
else
    annotation(f, 'textbox', [0.3 0.4 0.4 0.2], 'String', 'No significant SME cluster', 'EdgeColor', 'none', 'FontSize', 14);
end
save_fig(f, fullfile(C.dir.figs, 'Fig_SME'));

%% Fig_corr
f = figure('Units', 'centimeters', 'Position', [2 2 40 9], 'Color', 'w');
spec = {'Ispec', 'IM', R.Ispec_beh, 1; 'Ispec', 'LDI', R.Ispec_beh, 2; 'SME', 'IM', R.SME_beh, 1; ...
        'SME', 'LDI', R.SME_beh, 2; 'Ispec', 'SME', R.Ispec_SME, 1};
for k = 1:size(spec, 1)
    ax = subplot(1, 5, k);
    x = P.(spec{k,2}); y = P.(spec{k,1}); r = spec{k,3}; j = spec{k,4};
    scatter(ax, x(ya), y(ya), 20, col.YA, 'filled'); hold(ax, 'on');
    scatter(ax, x(~ya), y(~ya), 20, col.OA, 'filled');
    xlabel(ax, spec{k,2}); ylabel(ax, spec{k,1});
    title(ax, sprintf('partial r = %.2f, p_{FDR} = %.3f', r.r(j), r.pFDR(j)), 'FontWeight', 'normal');
end
save_fig(f, fullfile(C.dir.figs, 'Fig_corr'));

%% Fig_lags
f = figure('Units', 'centimeters', 'Position', [2 2 34 18], 'Color', 'w');
isI = strcmp(lagT.Type, 'Isim');
ax = subplot(2,2,1);
edges = 0:5:max(lagT.lag)+5;
histogram(ax, lagT.lag(isI), edges, 'FaceColor', [0.85 0.33 0.1]); hold(ax, 'on');
histogram(ax, lagT.lag(~isI), edges, 'FaceColor', [0.47 0.27 0.84]);
legend(ax, {'Isim pairs', 'Csim pairs'}); xlabel(ax, 'Lag (trial positions)'); ylabel(ax, 'Pairs (all persons)');
ax = subplot(2,2,2);
[g, lagBin] = findgroups(min(floor(lagT.lag/5)*5, 200));
mRef = splitapply(@mean, lagT.refSim, g);
plot(ax, lagBin, mRef, 'k.-'); xlabel(ax, 'Lag (5-trial bins, >= 200 pooled)'); ylabel(ax, 'Reference similarity (Fisher z)');
title(ax, 'Similarity of unrelated pairs by lag (= drift)');
ax = subplot(2,2,3);
group_dots(ax, {P.drift(ya), P.drift(~ya)}, [1 2], {col.YA, col.OA});
set(ax, 'XTick', [1 2], 'XTickLabel', {'Younger', 'Older'}); ylabel(ax, 'reference(Isim) - reference(Csim)');
title(ax, sprintf('Drift effect in the Ispec cluster (mean Ispec %.4f)', mean(P.Ispec)));
ax = subplot(2,2,4);
histogram(ax, lagT.refPosDist, 30, 'FaceColor', [0.5 0.5 0.5]);
xlabel(ax, 'Position distance of the farthest reference pair'); ylabel(ax, 'Pairs');
save_fig(f, fullfile(C.dir.figs, 'Fig_lags'));
fprintf('Figures saved to %s\n', C.dir.figs);
end

function group_dots(ax, vals, xs, cols)
hold(ax, 'on');
for k = 1:numel(vals)
    v = vals{k}; jit = (rand(size(v)) - 0.5) * 0.3;
    scatter(ax, xs(k) + jit, v, 18, cols{k}, 'filled', 'MarkerFaceAlpha', 0.6);
    plot(ax, xs(k) + [-0.25 0.25], [mean(v) mean(v)], 'k-', 'LineWidth', 2);
    se = std(v) / sqrt(numel(v));
    plot(ax, [xs(k) xs(k)], mean(v) + [-se se], 'k-', 'LineWidth', 1.5);
end
xlim(ax, [min(xs) - 0.7, max(xs) + 0.7]); box(ax, 'off');
end

function topo(ax, label, values, C, haveFT, ttl)
if ~haveFT
    axis(ax, 'off'); text(ax, 0.5, 0.5, 'layout file not found', 'HorizontalAlignment', 'center'); return
end
L = load(C.layoutFile); fn = fieldnames(L);
d = struct('label', {label(:)}, 'avg', double(values(:)), 'dimord', 'chan_time', 'time', 0);
cfg = []; cfg.layout = L.(fn{1}); cfg.comment = 'no'; cfg.marker = 'on'; cfg.interactive = 'no';
axes(ax); ft_topoplotER(cfg, d); title(ax, ttl);
end

function save_fig(f, name)
set(f, 'PaperPositionMode', 'auto');
print(f, [name '.png'], '-dpng', '-r300');
print(f, [name '.pdf'], '-dpdf', '-painters');
close(f);
end

function m = bwr(n)
% blue-white-red colour map
h = floor(n/2);
m = [[linspace(0,1,h)' linspace(0,1,h)' ones(h,1)]; [ones(n-h,1) linspace(1,0,n-h)' linspace(1,0,n-h)']];
end

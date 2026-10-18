function p1_figures()
% P1_FIGURES  Group-level figures (PLOT step only; no statistics computed).
% Reads the outputs of g1-g4 and draws:
%   Fig3_SimMatrices  group-average Isim and Csim matrices, YA and OA, as the
%                     original Fig3 (all channels, display smoothing 32 ms)
%   Fig_Ispec   t-map (mean over cluster channels) with cluster outline,
%               cluster topography, Isim/Csim and Ispec per age group,
%               split halves
%   Fig_SME     t-map (mean over SME cluster channels) with outline,
%               topography, SME per age group
%   Fig_GammaSync  gamma power change from baseline (mean over cluster
%               channels, YA and OA) with cluster outline, topography,
%               gSync per age group
%   Fig_GammaSME   the same for the gamma SME (log10 power, old - similar)
%   Fig_corr    between-person correlations (numbers from g4)
%   Fig_lags    lags of Isim / Csim pairs, baseline similarity by lag,
%               removed drift per person, pair positions in the session
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
title(ax, sprintf('Ispec group t, Isim vs Csim (mean over %d cluster channels)', sum(chI)));
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

%% Fig3_SimMatrices (group-average similarity) -- as the original Fig3
%  (rsa_plots.m): mean over ALL channels, display-only 2-D Gaussian smoothing
%  (sigma 32 ms, display_smooth_sigma.m), same axes, colour map and colour
%  limits. Statistics are unaffected by the smoothing.
GM = load(fullfile(C.dir.group, 'ispec_group_means.mat'));
sigmaSmp = 0.032 / mean(diff(GM.time));
disp3 = @(X) smooth2d(squeeze(mean(double(X), 1)), sigmaSmp);   % mean over channels
D = struct('YA_I', disp3(GM.GA.YA.Isim), 'YA_C', disp3(GM.GA.YA.Csim), ...
           'OA_I', disp3(GM.GA.OA.Isim), 'OA_C', disp3(GM.GA.OA.Csim));
allv = [D.YA_I(:); D.YA_C(:); D.OA_I(:); D.OA_C(:)];
clim3 = [max(prctile(allv, 5), 0) prctile(allv, 98)];
if clim3(2) <= clim3(1), clim3 = [prctile(allv, 5) prctile(allv, 98)]; end   % corrected values can be <= 0
f = figure('Units', 'centimeters', 'Position', [2 2 28 15], 'Color', 'w');
panels = {'YA_I', 'Younger adults - Item similarity'; 'YA_C', 'Younger adults - Category similarity'; ...
          'OA_I', 'Older adults - Item similarity';   'OA_C', 'Older adults - Category similarity'};
for k = 1:4
    ax = subplot(2, 2, k);
    imagesc(ax, GM.time, GM.time, D.(panels{k,1})); set(ax, 'YDir', 'normal');
    colormap(ax, flipud(autumn(256))); caxis(ax, clim3); colorbar(ax);
    xlabel(ax, 'Encoding 2nd (time, s)'); ylabel(ax, 'Encoding 1st (time, s)'); title(ax, panels{k,2});
end
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

%% Fig_GammaSync, Fig_GammaSME
gspec = {'gsync_cluster.mat', 'gSync', 'log10 power - baseline', 'Gamma sync (log10 power - baseline)', 'Fig_GammaSync'; ...
         'gsme_cluster.mat',  'gSME',  'old - similar (log10 power)', 'Gamma SME (log10 power, old - similar)', 'Fig_GammaSME'};
for g = 1:size(gspec, 1)
    gamma_fig(load(fullfile(C.dir.group, gspec{g,1})), P.(gspec{g,2}), R.age, gspec{g,2}, gspec(g,3:5), ya, col, C, haveFT);
end

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
plot(ax, lagBin, mRef, 'k.-'); xlabel(ax, 'Lag (5-trial bins, >= 200 pooled)'); ylabel(ax, 'Baseline similarity (Fisher z)');
title(ax, 'Baseline similarity by lag (= drift captured)');
ax = subplot(2,2,3);
group_dots(ax, {P.drift(ya), P.drift(~ya)}, [1 2], {col.YA, col.OA});
set(ax, 'XTick', [1 2], 'XTickLabel', {'Younger', 'Older'}); ylabel(ax, 'baseline(Isim) - baseline(Csim)');
title(ax, sprintf('Removed by the baseline in the Ispec cluster (mean Ispec %.4f)', mean(P.Ispec)));
ax = subplot(2,2,4);
histogram(ax, lagT.pos(strcmp(lagT.Type, 'Isim')), 30, 'FaceColor', [0.85 0.33 0.1]); hold(ax, 'on');
histogram(ax, lagT.pos(~strcmp(lagT.Type, 'Isim')), 30, 'FaceColor', [0.47 0.27 0.84]);
xlabel(ax, 'Pair position in the session (midpoint)'); ylabel(ax, 'Pairs (all persons)');
save_fig(f, fullfile(C.dir.figs, 'Fig_lags'));
fprintf('Figures saved to %s\n', C.dir.figs);
end

function gamma_fig(GG, v, age, vname, txt, ya, col, C, haveFT)
% YA and OA maps (mean over cluster channels, or all channels without a
% cluster) with the cluster outline, topography, per-person values by age
f = figure('Units', 'centimeters', 'Position', [2 2 40 10], 'Color', 'w');
anyG = any(GG.mask(:));
if anyG, chG = any(reshape(GG.mask, size(GG.mask,1), []), 2); else, chG = true(numel(GG.stat.label), 1); end
[~, iC] = ismember(GG.stat.label(chG), GG.label);
iF = arrayfun(@(x) find(abs(GG.freq - x) < 1e-6, 1), GG.stat.freq);
iT = arrayfun(@(x) find(abs(GG.time - x) < 1e-6, 1), GG.stat.time);
mG = {squeeze(mean(GG.GA_YA(iC, iF, iT), 1)), squeeze(mean(GG.GA_OA(iC, iF, iT), 1))};
c = max(cellfun(@(m) max(abs(m(:))), mG));
grp = {'Younger', 'Older'};
for k = 1:2
    ax = subplot(1,4,k);
    imagesc(ax, GG.stat.time, GG.stat.freq, mG{k}); axis(ax, 'xy'); hold(ax, 'on');
    if anyG
        contour(ax, GG.stat.time, GG.stat.freq, double(squeeze(any(GG.mask(chG,:,:), 1))), [0.5 0.5], 'k', 'LineWidth', 1.2);
    end
    colormap(ax, bwr(256)); caxis(ax, [-c c]); colorbar(ax);
    xlabel(ax, 'Time (s)'); ylabel(ax, 'Frequency (Hz)');
    title(ax, sprintf('%s: %s, %d channels', grp{k}, txt{1}, sum(chG)));
end
if anyG
    ax = subplot(1,4,3);
    topo(ax, GG.stat.label, sum(reshape(GG.mask, size(GG.mask,1), []), 2), C, haveFT, 'cluster points per channel');
    ax = subplot(1,4,4);
    group_dots(ax, {v(ya), v(~ya)}, [1 2], {col.YA, col.OA});
    set(ax, 'XTick', [1 2], 'XTickLabel', grp); ylabel(ax, txt{2});
    a = age.(vname); title(ax, sprintf('t(%d) = %.2f, p = %.3f, d = %.2f', a.df, a.t, a.p, a.d));
else
    annotation(f, 'textbox', [0.55 0.4 0.4 0.2], 'String', 'No significant positive cluster (all channels shown)', 'EdgeColor', 'none', 'FontSize', 12);
end
save_fig(f, fullfile(C.dir.figs, txt{3}));
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

function out = smooth2d(mat, sigmaSmp)
% 2-D Gaussian smoothing for DISPLAY only (unchanged from rsa_plots.m);
% normalised by the kernel mass inside the matrix; NaNs ignored.
if sigmaSmp <= 0, out = mat; return; end
h = ceil(3*sigmaSmp);
g = exp(-(-h:h).^2 / (2*sigmaSmp^2)); g = g / sum(g);
k = g' * g;
ok = ~isnan(mat); m = mat; m(~ok) = 0;
out = conv2(m, k, 'same') ./ conv2(double(ok), k, 'same');
out(~ok) = NaN;
end

function m = bwr(n)
% blue-white-red colour map
h = floor(n/2);
m = [[linspace(0,1,h)' linspace(0,1,h)' ones(h,1)]; [ones(n-h,1) linspace(1,0,n-h)' linspace(1,0,n-h)']];
end

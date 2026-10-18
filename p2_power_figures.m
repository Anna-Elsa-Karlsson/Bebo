function p2_power_figures()
% P2_POWER_FIGURES  Figures of the power analyses only (PLOT step; no
% statistics computed). Needs only the outputs of g2 and g5, so it can run
% right after `submit_pipeline.sh power`:
%   group/sme_cluster.mat, group/gsync_cluster.mat, group/gsme_cluster.mat,
%   tables/person_sme.csv, tables/person_gpower.csv
% If stats/person_stats.mat (g4) exists, the YA vs OA test is written in the
% titles; otherwise the per-person values are shown without a test.
%
%   Fig_SME        alpha SME: group t-map of the 1st-level t-maps (mean
%                  over cluster channels) with cluster outline, topography,
%                  SME per age group
%   Fig_GammaSync  gamma sync: group t-map, YA and OA mean (log10 power -
%                  baseline), topography, gSync per age group
%   Fig_GammaSME   gamma SME: the same (log10 power, old - similar)
% Without a significant cluster the maps show the mean over ALL channels.
% R2017b-compatible (subplot, print). Output: C.dir.figs (png + pdf).

C = bebo_config();
if ~exist(C.dir.figs, 'dir'), mkdir(C.dir.figs); end
haveFT = exist(C.layoutFile, 'file') == 2;
if haveFT, bebo_init_fieldtrip(C); end
col = struct('YA', [30 144 255]/255, 'OA', [80 80 80]/255);

PS = readtable(fullfile(C.dir.tables, 'person_sme.csv'));
PG = readtable(fullfile(C.dir.tables, 'person_gpower.csv'));
age = struct();
fs = fullfile(C.dir.stats, 'person_stats.mat');
if exist(fs, 'file'), S = load(fs, 'R'); age = S.R.age; end

spec = {'sme_cluster.mat',   PS.SME,   PS.YA == 1, 'SME',   'Alpha SME (log10 power, old - similar)',  '',                            'Fig_SME'; ...
        'gsync_cluster.mat', PG.gSync, PG.YA == 1, 'gSync', 'Gamma sync (log10 power - baseline)',      'log10 power - baseline',      'Fig_GammaSync'; ...
        'gsme_cluster.mat',  PG.gSME,  PG.YA == 1, 'gSME',  'Gamma SME (log10 power, old - similar)',   'old - similar (log10 power)', 'Fig_GammaSME'};
for k = 1:size(spec, 1)
    G = load(fullfile(C.dir.group, spec{k,1}));
    power_fig(G, spec{k,2}, spec{k,3}, age, spec(k,4:7), col, C, haveFT);
end
fprintf('Power figures saved to %s\n', C.dir.figs);
end

function power_fig(G, v, ya, age, txt, col, C, haveFT)
% txt = {value name, value label, map label ('' = no YA/OA maps), file name}
haveGA = ~isempty(txt{3}) && isfield(G, 'GA_YA');
nP = 3 + 2*haveGA;
f = figure('Units', 'centimeters', 'Position', [2 2 10*nP 10], 'Color', 'w');
anyC = any(G.mask(:));
if anyC, ch = any(reshape(G.mask, size(G.mask,1), []), 2); else, ch = true(numel(G.stat.label), 1); end
outline = double(squeeze(any(G.mask(ch,:,:), 1)));
chTxt = sprintf('%d channels', sum(ch)); if ~anyC, chTxt = 'all channels, no cluster'; end

% group t-map
ax = subplot(1, nP, 1);
tm = squeeze(mean(G.stat.stat(ch,:,:), 1));
imagesc(ax, G.stat.time, G.stat.freq, tm); axis(ax, 'xy'); hold(ax, 'on');
if anyC, contour(ax, G.stat.time, G.stat.freq, outline, [0.5 0.5], 'k', 'LineWidth', 1.2); end
colormap(ax, bwr(256)); c = max(abs(tm(:))); caxis(ax, [-c c]); colorbar(ax);
xlabel(ax, 'Time (s)'); ylabel(ax, 'Frequency (Hz)');
title(ax, sprintf('Group t (1st-level t-maps vs 0), %s', chTxt));

% YA and OA means in the original units (gamma only)
p = 2;
if haveGA
    [~, iC] = ismember(G.stat.label(ch), G.label);
    iF = arrayfun(@(x) find(abs(G.freq - x) < 1e-6, 1), G.stat.freq);
    iT = arrayfun(@(x) find(abs(G.time - x) < 1e-6, 1), G.stat.time);
    mG = {squeeze(mean(G.GA_YA(iC, iF, iT), 1)), squeeze(mean(G.GA_OA(iC, iF, iT), 1))};
    cl = max(cellfun(@(m) max(abs(m(:))), mG));
    grp = {'Younger', 'Older'};
    for g = 1:2
        ax = subplot(1, nP, p); p = p + 1;
        imagesc(ax, G.stat.time, G.stat.freq, mG{g}); axis(ax, 'xy'); hold(ax, 'on');
        if anyC, contour(ax, G.stat.time, G.stat.freq, outline, [0.5 0.5], 'k', 'LineWidth', 1.2); end
        colormap(ax, bwr(256)); caxis(ax, [-cl cl]); colorbar(ax);
        xlabel(ax, 'Time (s)'); ylabel(ax, 'Frequency (Hz)');
        title(ax, sprintf('%s: %s', grp{g}, txt{3}));
    end
end

if anyC
    ax = subplot(1, nP, p);
    topo(ax, G.stat.label, sum(reshape(G.mask, size(G.mask,1), []), 2), C, haveFT, 'cluster points per channel');
    ax = subplot(1, nP, p + 1);
    group_dots(ax, {v(ya), v(~ya)}, [1 2], {col.YA, col.OA});
    set(ax, 'XTick', [1 2], 'XTickLabel', {'Younger', 'Older'}); ylabel(ax, txt{2});
    if isfield(age, txt{1})
        a = age.(txt{1}); title(ax, sprintf('t(%d) = %.2f, p = %.3f, d = %.2f', a.df, a.t, a.p, a.d));
    else
        title(ax, 'per person (age test: run G4)');
    end
else
    ax = subplot(1, nP, p); axis(ax, 'off');
    text(ax, 0.5, 0.5, sprintf('No significant cluster\nin the expected direction'), 'HorizontalAlignment', 'center', 'FontSize', 12);
end
save_fig(f, fullfile(C.dir.figs, txt{4}));
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

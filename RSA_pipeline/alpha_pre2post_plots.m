function alpha_pre2post_plots(pdat, cfg_layout)
% Figure + stats for the alpha activation-vs-zero (pre-to-post) analysis.
% Port of Fig2_Alpha + Section 3 of BEBO_analyses.m, made R2017b-
% compatible for the grid (subplot instead of tiledlayout, caxis instead
% of clim, no mean(x,[1 2]), print instead of exportgraphics).
%
%   Row 1: baseline-corrected TFR (mean over subjects + channels), YA | OA
%   Row 2: activation-vs-zero t-values (mean over channels) with the
%          significant desync cluster outlined | topoplot of cluster-mean
%          t-values, black asterisks = channels in the significant cluster
%   Row 3: per-subject desync IN THE DESYNC CLUSTER, YA vs OA (raincloud) |
%          desync x IM scatter (partial r controlling age)
%   Stats: YA vs OA t-test on the cluster desync; desync x IM and
%          desync x LDI (partial r controlling age, within-group r, Fisher
%          z YA vs OA; BH-FDR within this 2-test family, separately for
%          partial r and Fisher z -- parallel to the Ispec correlations in
%          rsa_plots.m, which form their own family) -> text file
%
%   pdat       : struct with fields .alpha (Alpha_pre2post dir: GApow*.mat,
%                stat_alpha.mat), .perfFile (Performance_All.txt, matched
%                to ids_ALL by ID), .fig (output dir)
%
% (2026-10-05) Per-subject desync = mean baseline-corrected log power over
% ALL points (channel x frequency x time) of the significant desync
% cluster(s) -- where the desync is -- for the age comparison (YA vs OA)
% and desync x IM/LDI. NOT tested against 0 there (circular: the cluster
% was selected for desync vs 0). The cluster was chosen on the pooled
% sample (37 YA / 56 OA), which strictly leans the YA vs OA comparison
% slightly towards OA; the direct whole-head YA vs OA cluster test (run
% 2026-10-02, old_ageDiff/) showed no age difference either.
% Topoplot asterisks = the cluster's channels.
%   cfg_layout : path to the EEG layout .mat file

ft_available = exist('ft_topoplotER','file') == 2;

col_YA = [30 144 255]/255;
col_OA = [80 80 80]/255;
cmap_div = bwr_colormap();

ID_ya = ids_ALL('ya'); nYA = length(ID_ya);
ID_oa = ids_ALL('oa'); nOA = length(ID_oa);
nAll = nYA + nOA;
idx_YA = 1:nYA; idx_OA = nYA+1:nAll;
age_vec = [ones(nYA,1); 2*ones(nOA,1)];
allIDs = [ID_ya(:); ID_oa(:)];

% Behaviour (IM, LDI), matched to ids_ALL by ID
perf = readtable(pdat.perfFile);
[found, iPerf] = ismember(allIDs, perf.ID);
if ~all(found)
    error('alpha_pre2post_plots:perf', 'Performance file is missing subject(s): %s', mat2str(allIDs(~found)'));
end
IM = perf.IM(iPerf); LDI = perf.LDI(iPerf);

%% === Load ===
L = load(fullfile(pdat.alpha,'GApow.mat'));      GApow      = L.GApow;
L = load(fullfile(pdat.alpha,'GApow_YA.mat'));   GApow_YA   = L.GApow_YA;
L = load(fullfile(pdat.alpha,'GApow_OA.mat'));   GApow_OA   = L.GApow_OA;
L = load(fullfile(pdat.alpha,'stat_alpha.mat')); stat_alpha = L.stat_alpha;

if size(GApow.powspctrm,1) ~= nAll
    error('alpha_pre2post_plots:nSubj', 'GApow has %d subjects but ids_ALL has %d', size(GApow.powspctrm,1), nAll);
end

freq_vec = GApow.freq;
time_vec = GApow.time;

%% === Desync cluster + per-subject desync over the cluster's points ===
alpha_mask = alpha_cluster_mask(stat_alpha);            % chan x freq x time on stat_alpha's own (restricted) axes

% stat_alpha was tested on a restricted window (0.4-1.7 s, 6-30 Hz); GApow
% spans the full decomposed range -- restrict it to the test's axes, then
% average each subject's baseline-corrected log power (log ratio; NEGATIVE
% = desync) over the cluster points (plain mean).
[okC, iCh] = ismember(stat_alpha.label, GApow.label);
freqIdxStat = arrayfun(@(f) find(abs(freq_vec - f) < 1e-6, 1), stat_alpha.freq);
timeIdxStat = arrayfun(@(t) find(abs(time_vec - t) < 1e-6, 1), stat_alpha.time);
if ~all(okC), error('alpha_pre2post_plots:chan', 'stat_alpha channels missing in GApow'); end
pow_sub   = GApow.powspctrm(:, iCh, freqIdxStat, timeIdxStat);   % subj x chan x freq x time, aligned to alpha_mask
alpha_avg = reshape(pow_sub, nAll, []) * double(alpha_mask(:)) / sum(alpha_mask(:));

anyFreqSig     = squeeze(any(alpha_mask,2));                          % chan x time
timeTouchedAny = any(anyFreqSig,1);                                    % time points the cluster covers
clusterChans   = any(reshape(alpha_mask, size(alpha_mask,1), []), 2); % channels in the cluster at all

% Display matrices (R2017b: mean over one dim at a time)
GA_YA_display = squeeze(mean(mean(GApow_YA.powspctrm, 1), 2)); % freq x time
GA_OA_display = squeeze(mean(mean(GApow_OA.powspctrm, 1), 2));
stat_display  = squeeze(mean(stat_alpha.stat, 1));             % freq x time (restricted)
mask_display  = squeeze(any(alpha_mask, 1));

fprintf('Desync cluster: %.0f-%.0f Hz, %.2f-%.2f s, %d channels\n', ...
    min(stat_alpha.freq(any(squeeze(any(alpha_mask,1)),2))), max(stat_alpha.freq(any(squeeze(any(alpha_mask,1)),2))), ...
    min(stat_alpha.time(timeTouchedAny)), max(stat_alpha.time(timeTouchedAny)), sum(clusterChans));

%% === Figure ===
fig = figure('Units','centimeters','Position',[2 2 22 26],'Color','w','Visible','off');

ax = subplot(3,2,1);
imagesc(ax, time_vec, freq_vec, GA_YA_display); set(ax,'YDir','normal');
colormap(ax, cmap_div); caxis(ax, sym_clim(GA_YA_display));
xlabel(ax,'Time (s)'); ylabel(ax,'Frequency (Hz)'); title(ax,'Younger adults');
cb = colorbar(ax); ylabel(cb,'Power (log10, baseline-corrected)');

ax = subplot(3,2,2);
imagesc(ax, time_vec, freq_vec, GA_OA_display); set(ax,'YDir','normal');
colormap(ax, cmap_div); caxis(ax, sym_clim(GA_OA_display));
xlabel(ax,'Time (s)'); ylabel(ax,'Frequency (Hz)'); title(ax,'Older adults');
cb = colorbar(ax); ylabel(cb,'Power (log10, baseline-corrected)');

ax = subplot(3,2,3);
imagesc(ax, stat_alpha.time, stat_alpha.freq, stat_display); set(ax,'YDir','normal');
colormap(ax, cmap_div); caxis(ax, max(abs(stat_display(:)))*[-1 1]); hold(ax,'on');
contour(ax, stat_alpha.time, stat_alpha.freq, double(mask_display), [0.5 0.5], 'k', 'LineWidth', 1.5);
xlabel(ax,'Time (s)'); ylabel(ax,'Frequency (Hz)'); title(ax,'Activation vs. zero (desync)');
cb = colorbar(ax); ylabel(cb,'t-values');

ax = subplot(3,2,4);
if ft_available
    % Per channel: mean t over the cluster's time-frequency extent (cells
    % where the cluster is present on ANY channel), so every channel gets a
    % value. (The version in BEBO_analyses.m multiplied a chan x (freq*time)
    % matrix by the full chan*freq*time mask -> dimension error.)
    tf_mask   = mask_display(:);                                      % (freq*time) x 1
    stat_ct   = reshape(stat_alpha.stat, size(stat_alpha.stat,1), []); % chan x (freq*time)
    chan_tval = stat_ct(:, tf_mask) * ones(sum(tf_mask),1) / sum(tf_mask);
    clim_topo = max(abs(chan_tval)) * [-1 1];
    caxis(ax, clim_topo); colormap(ax, cmap_div);
    cb_topo = colorbar(ax); ylabel(cb_topo,'t-values');
    cb_topo.Position(3) = cb_topo.Position(3) * 0.5;
    dummy = struct(); dummy.label = stat_alpha.label;
    dummy.avg = chan_tval; dummy.time = 0; dummy.dimord = 'chan_time';
    cfg_topo = []; cfg_topo.layout = cfg_layout; cfg_topo.marker = 'off';
    cfg_topo.style = 'fill'; cfg_topo.colormap = cmap_div; cfg_topo.axes = ax;
    cfg_topo.comment = 'no'; cfg_topo.zlim = clim_topo;
    cfg_topo.highlight        = 'on';
    cfg_topo.highlightchannel = stat_alpha.label(clusterChans);   % the alpha cluster's own channels
    cfg_topo.highlightsymbol  = '*';
    cfg_topo.highlightsize    = 8;
    cfg_topo.highlightcolor   = [0 0 0];
    ft_topoplotER(cfg_topo, dummy);
    title(ax,{'Topographic distribution (mean t over cluster)', sprintf('* cluster channels (%d)', sum(clusterChans))},'FontSize',9);
else
    axis(ax,'off'); text(ax,0.5,0.5,'FieldTrip required','HorizontalAlignment','center');
end

% Desync x behaviour (partial r controlling age; within-group r; Fisher z)
behav = {IM, LDI}; blabel = {'IM', 'LDI'};
r_part = nan(1,2); p_part = nan(1,2); r_YA = nan(1,2); r_OA = nan(1,2); z_f = nan(1,2); p_f = nan(1,2);
for k = 1:2
    y = behav{k};
    [r_part(k), p_part(k)] = partialcorr(alpha_avg, y, age_vec);
    r_YA(k) = corr(alpha_avg(idx_YA), y(idx_YA));
    r_OA(k) = corr(alpha_avg(idx_OA), y(idx_OA));
    z_f(k) = (atanh(r_YA(k)) - atanh(r_OA(k))) / sqrt(1/(nYA-3) + 1/(nOA-3));
    p_f(k) = 2 * (1 - normcdf(abs(z_f(k))));
end
p_part_fdr = fdr_bh(p_part);
p_f_fdr    = fdr_bh(p_f);

ax = subplot(3,2,5);
raincloud_compare(alpha_avg(idx_YA), alpha_avg(idx_OA), col_YA, col_OA, ax);
title(ax,'Alpha desync (desync cluster)'); ylabel(ax,'Desync, cluster mean (log ratio)');
set(ax,'XTick',[1 2],'XTickLabel',{'Younger','Older'});

ax = subplot(3,2,6);
corr_scatter(alpha_avg, IM, age_vec, col_YA, col_OA, ax, 'Desync, cluster mean (log ratio)', 'Item memory (IM)', ...
    sprintf('Partial r=%.2f, p=%.3f (p_{FDR}=%.3f)', r_part(1), p_part(1), p_part_fdr(1)));

if ~exist(pdat.fig,'dir'), mkdir(pdat.fig); end
savefig_all(fig, fullfile(pdat.fig,'Fig_Alpha_pre2post'));


%% === Stats: YA vs OA on the cluster desync ===
[~, p_age, ~, st] = ttest2(alpha_avg(idx_YA), alpha_avg(idx_OA));
d_age = (mean(alpha_avg(idx_YA)) - mean(alpha_avg(idx_OA))) / ...
    sqrt(((nYA-1)*var(alpha_avg(idx_YA)) + (nOA-1)*var(alpha_avg(idx_OA))) / (nAll-2));

fid = fopen(fullfile(pdat.fig,'Alpha_pre2post_stats.txt'),'w');
fprintf(fid,'Alpha activation vs. zero -- %s\n\n', datestr(now));
p_crit = stat_alpha.cfg.alpha;
negP = []; if isfield(stat_alpha,'negclusters') && ~isempty(stat_alpha.negclusters), negP = [stat_alpha.negclusters.prob]; end % absent when no cluster formed
fprintf(fid,'Significant negative (desync) cluster(s): %s (p_crit=%.3f)\n', mat2str(negP(negP < p_crit),3), p_crit);
fprintf(fid,'Cluster extent: %.0f-%.0f Hz, %.2f-%.2f s, %d channels\n\n', ...
    min(stat_alpha.freq(any(squeeze(any(alpha_mask,1)),2))), max(stat_alpha.freq(any(squeeze(any(alpha_mask,1)),2))), ...
    min(stat_alpha.time(timeTouchedAny)), max(stat_alpha.time(timeTouchedAny)), sum(clusterChans));
fprintf(fid,'Cluster channels: %s\n\n', strjoin(stat_alpha.label(clusterChans)', ' '));
fprintf(fid,'Per-subject desync below = mean over ALL points of the significant desync cluster(s) (log ratio; NEGATIVE = desync).\n');
fprintf(fid,'Not tested against 0 here (circular: the cluster was selected for desync vs 0).\n\n');
fprintf(fid,'--- Age difference on desync (desync cluster) ---\n');
fprintf(fid,'Younger: M=%.3f SD=%.3f  Older: M=%.3f SD=%.3f\n', ...
    mean(alpha_avg(idx_YA)), std(alpha_avg(idx_YA)), mean(alpha_avg(idx_OA)), std(alpha_avg(idx_OA)));
fprintf(fid,'t(%d)=%.3f, p=%.4f, Cohen d=%.3f\n\n', st.df, st.tstat, p_age, d_age);
fprintf(fid,'--- Desync x behaviour (N=%d: %d YA, %d OA) ---\n', nAll, nYA, nOA);
fprintf(fid,'Desync = per-subject power change vs. baseline over the desync cluster (NEGATIVE = stronger desync),\n');
fprintf(fid,'so a NEGATIVE r means stronger desync goes with a higher score.\n');
fprintf(fid,'Partial r controlling age (full sample); within-group Pearson r; Fisher z-test YA vs OA.\n');
fprintf(fid,'FDR (Benjamini-Hochberg) within this family (2 tests), separately for partial r and Fisher z.\n\n');
for k = 1:2
    fprintf(fid,'Desync ~ %s\n', blabel{k});
    fprintf(fid,'   partial r=%.3f, p=%.4f, p_FDR=%.4f\n', r_part(k), p_part(k), p_part_fdr(k));
    fprintf(fid,'   YA r=%.3f (n=%d), OA r=%.3f (n=%d)\n', r_YA(k), nYA, r_OA(k), nOA);
    fprintf(fid,'   Fisher z=%.2f, p=%.4f, p_FDR=%.4f\n\n', z_f(k), p_f(k), p_f_fdr(k));
end
fclose(fid);
writetable(table(stat_alpha.label(clusterChans), 'VariableNames', {'channel'}), fullfile(pdat.fig, 'Alpha_cluster_channels.csv'));
fprintf('Alpha pre2post figure + stats saved -> %s\n', pdat.fig);

end


%% =========================================================
%% LOCAL HELPER FUNCTIONS (same as rsa_plots.m)
%% =========================================================

function corr_scatter(x, y, age_vec, col_YA, col_OA, ax, xlab, ylab, title_str)
% As in rsa_plots.m; regression line spans each group's own x-range
iY = age_vec == 1; iO = age_vec == 2;
hold(ax,'on');
scatter(ax, x(iO), y(iO), 30, col_OA, 'x', 'LineWidth', 1.2);
scatter(ax, x(iY), y(iY), 30, col_YA, 'o', 'filled', 'MarkerFaceAlpha', 0.6);
for g = 1:2
    if g == 1, sel = iY; c = col_YA; else, sel = iO; c = col_OA; end
    xs = x(sel); ys = y(sel); ok = ~isnan(xs) & ~isnan(ys);
    pf = polyfit(xs(ok), ys(ok), 1); xl = [min(xs(ok)) max(xs(ok))];
    plot(ax, xl, polyval(pf, xl), '-', 'Color', c, 'LineWidth', 1.2);
end
xlabel(ax, xlab); ylabel(ax, ylab); title(ax, title_str, 'FontSize', 9); box(ax,'off');
legend(ax, {'Older','Younger'}, 'Location','best', 'Box','off', 'FontSize',8);
end

function p_adj = fdr_bh(p)
% Benjamini-Hochberg adjusted p-values (step-up), same order as input
m = numel(p); [ps, order] = sort(p(:)');
adj = ps .* m ./ (1:m);
adj = min(1, fliplr(cummin(fliplr(adj))));
p_adj = nan(size(p)); p_adj(order) = adj;
end

function cl = sym_clim(data)
v = max(abs(prctile(data(:), [2 98])));
cl = [-v v];
end

function raincloud_compare(data_YA, data_OA, col_YA, col_OA, ax)
hold(ax,'on');
raincloud_single(data_YA(:), 1, col_YA, ax);
raincloud_single(data_OA(:), 2, col_OA, ax);
xlim(ax,[0.3 2.7]); box(ax,'off');
end

function raincloud_single(data, xpos, col, ax)
data = data(~isnan(data));
[f, xi] = ksdensity(data, 'NumPoints', 200);
f = f / max(f) * 0.35;
patch(ax, [xpos + f, xpos*ones(1,length(f))], [xi, fliplr(xi)], col, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
q1 = prctile(data,25); q3 = prctile(data,75); med = median(data);
wlo = max(min(data), q1-1.5*(q3-q1)); whi = min(max(data), q3+1.5*(q3-q1));
bw = 0.07;
patch(ax, xpos+[-bw bw bw -bw], [q1 q1 q3 q3], col, 'FaceAlpha', 0.85, 'EdgeColor', col*0.6, 'LineWidth', 0.8);
line(ax, xpos+[-bw bw], [med med], 'Color','w', 'LineWidth', 2);
line(ax, [xpos xpos], [wlo q1], 'Color', col*0.6, 'LineWidth', 1);
line(ax, [xpos xpos], [q3 whi], 'Color', col*0.6, 'LineWidth', 1);
jitter = (rand(length(data),1)-0.5)*0.15;
scatter(ax, xpos-0.22+jitter, data, 18, col, 'filled', 'MarkerFaceAlpha', 0.45, 'MarkerEdgeColor','none');
end

function savefig_all(fig, filepath)
% exportgraphics (R2020a+) is not available on the grid's R2017b
set(fig, 'Units','centimeters'); pos = get(fig,'Position');
set(fig, 'PaperUnits','centimeters', 'PaperSize', pos(3:4), 'PaperPosition', [0 0 pos(3:4)]); % PDF page = figure size (otherwise cut off)
print(fig, [filepath '.pdf'], '-dpdf');
print(fig, [filepath '.png'], '-dpng', '-r300');
print(fig, [filepath '.tif'], '-dtiff', '-r300');
fprintf('Saved: %s (.pdf/.png/.tif)\n', filepath);
end

function cmap = bwr_colormap(n)
if nargin < 1, n = 256; end
half = floor(n/2);
blue_to_white = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
white_to_red = [ones(n-half,1), linspace(1,0,n-half)', linspace(1,0,n-half)'];
cmap = [blue_to_white; white_to_red];
end

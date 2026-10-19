function plot_singleTrial_distributions(freqLo, freqHi)
% Per-subject distributions of the single-trial values that go into the
% glmer: Isim (T_item.Sim, similarity within the Ispec cluster) and alpha
% power (T_item.pow1st, log10 alpha in the Ispec ROI at presentation 1).
% Diagnostic only -- reads the per-subject output of
% rsa_singleTrial_merged.m, writes nothing the pipeline uses.
%
% Output (Results/Plots/SingleTrialDistributions/):
%   Dist_Isim_page<k>.png/.pdf    histograms, one panel per subject
%   Dist_Alpha_page<k>.png/.pdf
%   singleTrial_distribution_summary.csv   per subject: n, mean, SD,
%                                          median, min, max, skewness,
%                                          n NaN, for both measures
%
% Usage on the grid (from RSA_pipeline/):
%   matlab -r "plot_singleTrial_distributions(30,100); exit"
% R2017b-compatible (subplot, print; no tiledlayout/exportgraphics).

if nargin < 2, freqLo = 30; freqHi = 100; end

pipelineRoot = fileparts(mfilename('fullpath'));
addpath(pipelineRoot); % ids_ALL.m
variantLabel = sprintf('Ispec_%dto%dHz', freqLo, freqHi);
inDir  = fullfile(pipelineRoot, '..', 'Results', 'SingleTrialMerged', variantLabel);
outDir = fullfile(pipelineRoot, '..', 'Results', 'Plots', 'SingleTrialDistributions');
if ~exist(outDir,'dir'), mkdir(outDir); end

ids = ids_ALL('all'); ids = ids(:);
nSubj = numel(ids);
col_YA = [30 144 255]/255;
col_OA = [80 80 80]/255;

%% === Load ===
isim  = cell(nSubj,1);
alpha = cell(nSubj,1);
for s = 1:nSubj
    f = fullfile(inDir, sprintf('%d_singleTrial_within-item.mat', ids(s)));
    if ~exist(f,'file')
        warning('Missing %s -- subject left empty', f);
        continue
    end
    L = load(f, 'T_item');
    isim{s}  = L.T_item.Sim;
    alpha{s} = L.T_item.pow1st;
end

%% === Summary table ===
stats = @(x) [sum(~isnan(x)), mean(x,'omitnan'), std(x,'omitnan'), median(x,'omitnan'), ...
              min(x), max(x), skewness(x(~isnan(x))), sum(isnan(x))];
S_i = nan(nSubj,8); S_a = nan(nSubj,8);
for s = 1:nSubj
    if ~isempty(isim{s}),  S_i(s,:) = stats(isim{s});  end
    if ~isempty(alpha{s}), S_a(s,:) = stats(alpha{s}); end
end
names = {'n','mean','SD','median','min','max','skew','nNaN'};
T = array2table([ids, double(ids < 4000), S_i, S_a], 'VariableNames', ...
    [{'ID','YA'}, strcat('Isim_', names), strcat('Alpha_', names)]);
writetable(T, fullfile(outDir, 'singleTrial_distribution_summary.csv'));

%% === Plots: one panel per subject, 6 x 8 = 48 per page ===
plot_pages(isim,  ids, 'Isim (Fisher z, within cluster)', 'Dist_Isim',  outDir, col_YA, col_OA);
plot_pages(alpha, ids, 'Alpha power (log10, pres. 1)',    'Dist_Alpha', outDir, col_YA, col_OA);

fprintf('Done: %d subjects -> %s\n', nSubj, outDir);

end


function plot_pages(vals, ids, xlab, fileStem, outDir, col_YA, col_OA)
nRow = 6; nCol = 8; perPage = nRow*nCol;
nSubj = numel(ids);
nPage = ceil(nSubj/perPage);

% Common x-range across subjects (1st-99th percentile of all values), so
% panels are comparable; values outside are clipped into the edge bins.
allv = vertcat(vals{:}); allv = allv(~isnan(allv));
if isempty(allv), warning('No data for %s', fileStem); return; end
lims = prctile(allv, [1 99]);
edges = linspace(lims(1), lims(2), 26);

for pg = 1:nPage
    fig = figure('Units','centimeters','Position',[1 1 42 30],'Color','w','Visible','off');
    for k = 1:perPage
        s = (pg-1)*perPage + k;
        if s > nSubj, break; end
        ax = subplot(nRow, nCol, k);
        x = vals{s};
        if isempty(x)
            axis(ax,'off'); text(ax,0.5,0.5,sprintf('%d\nmissing',ids(s)),'HorizontalAlignment','center');
            continue
        end
        x = x(~isnan(x));
        x = min(max(x, edges(1)), edges(end)); % clip into edge bins
        if ids(s) < 4000, c = col_YA; else, c = col_OA; end
        histogram(ax, x, edges, 'FaceColor', c, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
        hold(ax,'on');
        yl = ylim(ax);
        plot(ax, median(x)*[1 1], yl, 'k-', 'LineWidth', 1); % median
        title(ax, sprintf('%d (n=%d)', ids(s), numel(x)), 'FontSize', 7);
        set(ax, 'FontSize', 6, 'XLim', edges([1 end])); box(ax,'off');
    end
    annotation(fig, 'textbox', [0 0.965 1 0.03], 'String', ...
        sprintf('%s -- page %d/%d (blue = YA, grey = OA, line = median, x clipped to 1st-99th pct of all subjects)', xlab, pg, nPage), ...
        'EdgeColor','none','HorizontalAlignment','center','FontSize',10);
    fp = fullfile(outDir, sprintf('%s_page%d', fileStem, pg));
    print(fig, [fp '.png'], '-dpng', '-r200');
    print(fig, [fp '.pdf'], '-dpdf', '-fillpage');
    close(fig);
    fprintf('Saved %s\n', fp);
end
end

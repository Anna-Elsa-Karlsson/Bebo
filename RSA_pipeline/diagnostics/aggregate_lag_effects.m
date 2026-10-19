function aggregate_lag_effects()
% DIAGNOSTIC ONLY -- does not touch the production pipeline. Aggregates
% test_lag_effects_onesubj.m's per-subject output across the full sample
% (2026-08-27), after that array job has completed for all subjects.
%
% Two outputs:
%
% (1) Grouped bar plot: 3 condition-groups (within-item, within-cat,
%     between-cat) x 3 bars each (real, actVsBL, blVsBl -- all RAW,
%     uncorrected values, per user request). Error bars = SEM ACROSS
%     SUBJECTS (per-subject mean computed first, then group mean +/- SEM
%     -- not a pooled distribution across all trial-pairs, which would
%     understate true between-subject variability).
%
% (2) Group-level lag decay curve: each subject's binned
%     (lag -> mean raw similarity) curve from test_lag_effects_onesubj.m,
%     averaged bin-by-bin across subjects (mean +/- SEM across subjects
%     per bin). Bins aligned by INDEX (1-9, matching the fixed edges
%     array used in test_lag_effects_onesubj.m: [0 5 10 20 40 80 160 320
%     640 max(lag)+1]) -- subjects with fewer populated bins are
%     NaN-padded rather than misaligned.
%
% Also reports, per condition, a one-sample t-test of each subject's own
% corr(lag, real) against zero -- does trial lag reliably predict raw
% similarity across the population, not just in one subject's data.
%
% Reads from Results/BaselineDiagnostic/<subj>_lag_effects.mat (table T)
% and Results/BaselineDiagnostic/<subj>_lagDecayCurve.mat, written by
% test_lag_effects_onesubj.m. Saves aggregate outputs to the same folder.
%
% Usage: matlab -r "aggregate_lag_effects; exit"

pipelineRoot = fileparts(mfilename('fullpath'));
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/');
ft_defaults

diagDir = fullfile(pipelineRoot, 'Results', 'BaselineDiagnostic');
subjects = ids_ALL('all');
conds = {'within-item','within-cat','between-cat'};
measures = {'realVal','actVsBL','blVsBl'};
measureLabels = {'real (activ-activ)','actVsBL (activ-baseline)','blVsBl (baseline-baseline)'};

%% (1) Per-subject condition x measure means, for the bar plot

nSubj = length(subjects);
subjMeans = nan(nSubj, length(conds), length(measures)); % [subj x cond x measure]
subjLagCorr = nan(nSubj, length(conds)); % corr(lag, realVal) per subject per condition
included = false(nSubj,1);

for s = 1:nSubj
    curr_subj = num2str(subjects(s));
    fpath = fullfile(diagDir, [curr_subj '_lag_effects.mat']);
    if ~exist(fpath,'file')
        warning('Missing lag_effects file for subject %s -- skipping', curr_subj);
        continue
    end
    L = load(fpath, 'T');
    T = L.T;
    included(s) = true;

    for c = 1:length(conds)
        rows = strcmp(T.condition, conds{c});
        for m = 1:length(measures)
            subjMeans(s,c,m) = mean(T.(measures{m})(rows), 'omitnan');
        end
        subjLagCorr(s,c) = corr(T.lag(rows), T.realVal(rows), 'rows','pairwise');
    end
end

fprintf('Aggregating over %d/%d subjects with available data\n', sum(included), nSubj);
subjMeans = subjMeans(included,:,:);
subjLagCorr = subjLagCorr(included,:);
nIncluded = sum(included);

groupMean = squeeze(mean(subjMeans,1,'omitnan'));          % [cond x measure]
groupSEM  = squeeze(std(subjMeans,0,1,'omitnan')) / sqrt(nIncluded); % [cond x measure]

%% Bar plot: 3 conditions x 3 measures, error bars = SEM across subjects

fig = figure('Visible','off','Units','centimeters','Position',[2 2 22 16]);
b = bar(groupMean, 'grouped');
hold on;
[nGroups, nBars] = size(groupMean);
groupWidth = min(0.8, nBars/(nBars+1.5));
for m = 1:nBars
    xpos = (1:nGroups) - groupWidth/2 + (2*m-1) * groupWidth / (2*nBars);
    errorbar(xpos, groupMean(:,m), groupSEM(:,m), 'k', 'linestyle','none', 'LineWidth',1);
end
set(gca,'XTickLabel', conds);
ylabel('Mean raw correlation (Fisher-z)');
legend(measureLabels, 'Location','best','Box','off');
title(sprintf('Raw similarity by condition and comparison type (n=%d subjects)', nIncluded));
box off;

outDir = diagDir;
print(fig, fullfile(outDir, 'groupLevel_rawSimilarity_barplot.png'), '-dpng', '-r200');
fprintf('Saved bar plot -> %s\n', fullfile(outDir, 'groupLevel_rawSimilarity_barplot.png'));

fprintf('\n=== Group-level raw similarity means (n=%d subjects) ===\n', nIncluded);
fprintf('%-14s', ''); fprintf('%22s', measureLabels{:}); fprintf('\n');
for c = 1:length(conds)
    fprintf('%-14s', conds{c});
    for m = 1:length(measures)
        fprintf('%13.4f (%.4f)', groupMean(c,m), groupSEM(c,m));
    end
    fprintf('\n');
end

%% Group-level test: does lag predict raw similarity, per condition, across subjects?

fprintf('\n=== Group-level lag-predicts-similarity test (one-sample t-test of per-subject r vs 0) ===\n');
for c = 1:length(conds)
    r_c = subjLagCorr(:,c);
    [~,p,~,stats] = ttest(r_c);
    fprintf('%-14s mean r=%.3f  t(%d)=%.3f  p=%.4f\n', conds{c}, mean(r_c,'omitnan'), stats.df, stats.tstat, p);
end

%% (2) Group-level lag decay curve, aligned by bin index

nBinsMax = 9; % matches the fixed edges array in test_lag_effects_onesubj.m: numel(edges)-1
allBinMean = nan(nSubj, nBinsMax);
allBinCenter = nan(nSubj, nBinsMax);

for s = 1:nSubj
    if ~included(s), continue; end
    curr_subj = num2str(subjects(s));
    fpath = fullfile(diagDir, [curr_subj '_lagDecayCurve.mat']);
    if ~exist(fpath,'file')
        warning('Missing lagDecayCurve file for subject %s -- skipping from decay curve aggregation', curr_subj);
        continue
    end
    L = load(fpath, 'binMean','binCenter');
    nB = length(L.binMean);
    allBinMean(s,1:nB) = L.binMean(:)';
    allBinCenter(s,1:nB) = L.binCenter(:)';
end

groupBinMean = mean(allBinMean,1,'omitnan');
groupBinSEM  = std(allBinMean,0,1,'omitnan') ./ sqrt(sum(~isnan(allBinMean),1));
groupBinCenter = mean(allBinCenter,1,'omitnan'); % average bin center across subjects, for x-position

validBins = ~isnan(groupBinMean);
fig2 = figure('Visible','off','Units','centimeters','Position',[2 2 20 14]);
errorbar(groupBinCenter(validBins), groupBinMean(validBins), groupBinSEM(validBins), '-o','LineWidth',1.2);
xlabel('Trial lag'); ylabel('Mean raw similarity (\pm SEM across subjects)');
title(sprintf('Group-level raw similarity vs. trial lag (n=%d subjects)', nIncluded));
box off;
print(fig2, fullfile(outDir, 'groupLevel_lagDecayCurve.png'), '-dpng', '-r200');
fprintf('\nSaved group-level lag decay curve -> %s\n', fullfile(outDir, 'groupLevel_lagDecayCurve.png'));

save(fullfile(outDir, 'groupLevel_lag_effects_aggregate.mat'), ...
    'subjMeans','subjLagCorr','groupMean','groupSEM','allBinMean','allBinCenter', ...
    'groupBinMean','groupBinSEM','groupBinCenter','included','subjects','-v7.3');
fprintf('Saved aggregate data -> %s\n', fullfile(outDir, 'groupLevel_lag_effects_aggregate.mat'));

end

function b1_behavior()
% B1_BEHAVIOR  Behaviour (analysis plan section 4): taken over UNCHANGED from
% RSA_pipeline/behavior_analysis.m (validated; output identical), only paths
% (bebo_config) and the ID function (bebo_ids) adapted. Independent of the
% EEG steps. Output: C.dir.figs/Fig_Behavior.*, C.dir.stats/Behavior_stats.txt
C = bebo_config();
pdat.perfFile = C.perfFile; pdat.fig = C.dir.figs; pdat.stats = C.dir.stats;
if ~exist(pdat.stats, 'dir'), mkdir(pdat.stats); end

% Behavioural figure + stats (MATLAB/R2017b port of section 1 of
% BEBO_behavior_glmer.R, i.e. Fig1_Behavior):
%
%   A: response proportions per item type (confusion matrix), YA | OA,
%      mean +/- SE (bars: x = item type, colour = response)
%
% COLUMN CONVENTION in Performance_All.txt: <response>2<itemtype>, e.g.
% old2sim = P(respond "old" | similar item), new2old = P(respond "new" |
% old item); item 'foil' = new item. Confirmed 2026-08-31 and by the data
% (proportions over responses sum to ~1 per item type only this way).
% NOTE: the original Fig1 (BEBO_behavior_glmer.R) read the first token as
% item type -- same bars, but its x-axis ("Response") and legend ("Item
% type") labels were swapped; corrected here.
%   B: IM and LDI per age group (raincloud)
%   C: category-specific false recognition (port of
%      BEBO_oldnessBias_categoryVsNew.m, 2026-08-31): P("old" | similar)
%      = old2sim vs P("old" | new) = old2foil per age group, and the
%      difference old2sim - old2foil (= age x item-type interaction test)
%   Stats: IM and LDI (Welch t + d); category-specific false recognition:
%      (1) age x item-type interaction = Welch t-test YA vs OA on the
%          difference old2sim - old2foil (Welch rather than the equal-
%          variance ANOVA interaction: variances of the difference differ
%          between age groups), then (2) simple effects: YA vs OA on
%          old2sim and on old2foil separately (Welch), Holm-corrected over
%          the two tests -> text file
%
% Only subjects in ids_ALL('all') are used (Performance_All.txt may
% contain more rows, e.g. excluded IDs); every ids_ALL subject must be
% present, matched by ID (not row order).
%
%   pdat : struct with fields .perfFile (path to Performance_All.txt),
%          .fig (output dir)

col_YA = [30 144 255]/255;
col_OA = [80 80 80]/255;
col_resp = [149 117 205; 76 175 80; 240 98 146]/255; % response Old, Similar, New (same colours as the R figure)

ID_ya = bebo_ids('ya'); ID_oa = bebo_ids('oa');
allIDs = [ID_ya(:); ID_oa(:)];
nYA = numel(ID_ya); nOA = numel(ID_oa);
idx_YA = 1:nYA; idx_OA = nYA+1:nYA+nOA;

%% === Load + restrict to ids_ALL ===
perfAll = readtable(pdat.perfFile);
[found, iRow] = ismember(allIDs, perfAll.ID);
if ~all(found)
    error('behavior_analysis:missing', 'Performance file is missing subject(s): %s', mat2str(allIDs(~found)'));
end
perf = perfAll(iRow,:);  % ids_ALL order, extra rows dropped
fprintf('Behaviour: %d of %d rows in %s used (ids_ALL: %d YA, %d OA)\n', height(perf), height(perfAll), pdat.perfFile, nYA, nOA);

% Confusion cells: <response>2<itemtype> (see header); item 'foil' = New
resps = {'old','sim','new'}; items = {'old','sim','foil'};
respLab = {'Old','Similar','New'}; itemLab = {'Old','Similar','New'};
P = nan(nYA+nOA, 3, 3); % subj x response x item  = P(response | item)
for r = 1:3
    for i = 1:3
        P(:,r,i) = perf.([resps{r} '2' items{i}]);
    end
end

%% === Figure ===
fig = figure('Units','centimeters','Position',[2 2 30 36],'Color','w','Visible','off');

% --- A: confusion bars, one panel per age group ---
grpIdx = {idx_YA, idx_OA}; grpName = {'Younger','Older'};
yMax = 0;
axA = gobjects(1,2);
for g = 1:2
    axA(g) = subplot(3,2,g); ax = axA(g); hold(ax,'on');
    M  = squeeze(mean(P(grpIdx{g},:,:),1));                              % response x item
    SE = squeeze(std(P(grpIdx{g},:,:),0,1)) / sqrt(numel(grpIdx{g}));
    bw = 0.25; offs = [-bw 0 bw];
    h = gobjects(1,3);
    for r = 1:3
        x = (1:3) + offs(r);                                              % x = item type
        h(r) = bar(ax, x, M(r,:), bw*0.9, 'FaceColor', col_resp(r,:), 'EdgeColor','none');
        errorbar(ax, x, M(r,:), SE(r,:), 'k', 'LineStyle','none', 'CapSize', 4);
    end
    yMax = max(yMax, max(M(:)+SE(:)));
    set(ax, 'XTick', 1:3, 'XTickLabel', itemLab, 'XLim', [0.5 3.5]); box(ax,'off');
    xlabel(ax,'Item type'); title(ax, grpName{g});
    if g == 1, ylabel(ax,'Proportion of responses'); end
    if g == 2, lg = legend(ax, h, respLab, 'Location','northeastoutside','Box','off'); title(lg, 'Response'); end
end
for g = 1:2, ylim(axA(g), [0 yMax*1.08]); end
text(axA(1), -0.12, 1.08, 'A', 'Units','normalized', 'FontSize',16, 'FontWeight','bold');

% --- B: IM and LDI rainclouds ---
meas = {'IM','LDI'};
for m = 1:2
    ax = subplot(3,2,2+m);
    v = perf.(meas{m});
    raincloud_compare(v(idx_YA), v(idx_OA), col_YA, col_OA, ax);
    set(ax,'XTick',[1 2],'XTickLabel',{'Younger','Older'});
    title(ax, meas{m});
    if m == 1
        ylabel(ax,'Score');
        text(ax, -0.12, 1.08, 'B', 'Units','normalized', 'FontSize',16, 'FontWeight','bold');
    end
end

% --- C: category-specific false recognition ---
fa_sim  = perf.old2sim;   % P("old" | similar item) -- category information available
fa_new  = perf.old2foil;  % P("old" | new item)     -- no category information
fa_diff = fa_sim - fa_new;

ax = subplot(3,2,5); hold(ax,'on');
raincloud_single(fa_sim(idx_YA), 1.0, col_YA, ax);           raincloud_single(fa_new(idx_YA), 1.7, col_YA*0.4+0.6, ax);
raincloud_single(fa_sim(idx_OA), 3.0, col_OA, ax);           raincloud_single(fa_new(idx_OA), 3.7, [170 170 170]/255, ax);
xlim(ax,[0.3 4.4]); box(ax,'off');
plot(ax, [2.35 2.35], ylim(ax), '--', 'Color',[0.75 0.75 0.75], 'LineWidth', 0.8);
set(ax,'XTick',[1.35 3.35],'XTickLabel',{'Younger','Older'});
h1 = patch(ax,NaN,NaN,col_YA,'EdgeColor','none'); h2 = patch(ax,NaN,NaN,col_YA*0.4+0.6,'EdgeColor','none');
legend(ax,[h1 h2], {'"old" | similar','"old" | new'},'Location','best','Box','off','FontSize',8);
ylabel(ax, 'P(respond "old")'); title(ax, 'False recognition by item type');
text(ax, -0.12, 1.08, 'C', 'Units','normalized', 'FontSize',16, 'FontWeight','bold');

ax = subplot(3,2,6);
raincloud_compare(fa_diff(idx_YA), fa_diff(idx_OA), col_YA, col_OA, ax);
set(ax,'XTick',[1 2],'XTickLabel',{'Younger','Older'});
ylabel(ax, '"old"|similar - "old"|new'); title(ax, 'Category-specific false recognition');

if ~exist(pdat.fig,'dir'), mkdir(pdat.fig); end
fp = fullfile(pdat.fig, 'Fig_Behavior');
print(fig, [fp '.pdf'], '-dpdf', '-fillpage');
print(fig, [fp '.png'], '-dpng', '-r300');
print(fig, [fp '.tif'], '-dtiff', '-r300');
fprintf('Saved: %s (.pdf/.png/.tif)\n', fp);

%% === Stats: IM and LDI, YA vs OA (Welch t-test, as in the R script) ===
fid = fopen(fullfile(pdat.stats, 'Behavior_stats.txt'), 'w');
fprintf(fid, 'BEBO behaviour -- %s\nN = %d (%d YA, %d OA), IDs from ids_ALL.m\n\n', datestr(now), nYA+nOA, nYA, nOA);
for m = 1:2
    v = perf.(meas{m}); ya = v(idx_YA); oa = v(idx_OA);
    [~, p, ~, st] = ttest2(ya, oa, 'Vartype', 'unequal');
    sp = sqrt(((nYA-1)*var(ya) + (nOA-1)*var(oa)) / (nYA+nOA-2));
    d  = (mean(ya) - mean(oa)) / sp;
    fprintf(fid, '--- %s: YA vs OA ---\nYA: M=%.3f SD=%.3f  OA: M=%.3f SD=%.3f\n', meas{m}, mean(ya), std(ya), mean(oa), std(oa));
    fprintf(fid, 'Welch t(%.1f)=%.3f, p=%.4f, Cohen d=%.3f\n\n', st.df, st.tstat, p, d);
end
fprintf(fid, '--- Response proportions P(response | item type) (M, SE) ---\n');
for g = 1:2
    fprintf(fid, '%s:\n', grpName{g});
    for i = 1:3
        for r = 1:3
            x = P(grpIdx{g},r,i);
            fprintf(fid, '  item %-7s -> response %-7s: M=%.3f SE=%.3f\n', itemLab{i}, respLab{r}, mean(x), std(x)/sqrt(numel(x)));
        end
    end
end
fprintf(fid, '\n=== CATEGORY-SPECIFIC FALSE RECOGNITION (old2sim vs old2foil) ===\n');
fprintf(fid, 'If OA''s "old" bias is driven by category familiarity, it should appear for\n');
fprintf(fid, 'similar items (same category) but not for new items (unstudied categories).\n\n');
% (1) Interaction: age x item type (similar vs new) = Welch t-test YA vs OA
%     on the difference score (main test -- unequal variances, see header)
fprintf(fid, '(1) AGE x ITEM-TYPE INTERACTION: YA vs OA on old2sim - old2foil (Welch t-test = main test)\n');
xYA = fa_diff(idx_YA); xOA = fa_diff(idx_OA);
[~, p_int, ~, st_int] = ttest2(xYA, xOA, 'Vartype','unequal');
d_int = (mean(xYA) - mean(xOA)) / sqrt(((numel(xYA)-1)*var(xYA) + (numel(xOA)-1)*var(xOA)) / (numel(xYA)+numel(xOA)-2));
p_lev = vartestn([xYA(:); xOA(:)], [zeros(numel(xYA),1); ones(numel(xOA),1)], 'TestType','LeveneAbsolute','Display','off');
fprintf(fid, 'YA: M=%.4f SD=%.4f  OA: M=%.4f SD=%.4f\n', mean(xYA), std(xYA), mean(xOA), std(xOA));
fprintf(fid, 'Welch t(%.1f)=%.3f, p=%.4f, Cohen d=%.3f  (Levene p=%.4f)\n\n', st_int.df, st_int.tstat, p_int, d_int, p_lev);

% (2) Simple effects: age difference per item type, Welch, Holm over the 2 tests
fprintf(fid, '(2) SIMPLE EFFECTS: YA vs OA per item type (Welch t-tests, Holm-corrected over 2 tests)\n');
lbl = {'old2sim  = P("old" | similar item)', 'old2foil = P("old" | new item)    '};
dat = {fa_sim, fa_new};
p_raw = nan(1,2); st_se = cell(1,2); d_se = nan(1,2);
for k = 1:2
    xYA = dat{k}(idx_YA); xOA = dat{k}(idx_OA);
    [~, p_raw(k), ~, st_se{k}] = ttest2(xYA, xOA, 'Vartype','unequal');
    d_se(k) = (mean(xYA) - mean(xOA)) / sqrt(((numel(xYA)-1)*var(xYA) + (numel(xOA)-1)*var(xOA)) / (numel(xYA)+numel(xOA)-2));
end
p_holm = holm(p_raw);
for k = 1:2
    xYA = dat{k}(idx_YA); xOA = dat{k}(idx_OA);
    fprintf(fid, '%s  YA: M=%.4f SD=%.4f  OA: M=%.4f SD=%.4f  Welch t(%.1f)=%.3f, p=%.4f, p_Holm=%.4f, d=%.3f\n', ...
        lbl{k}, mean(xYA), std(xYA), mean(xOA), std(xOA), st_se{k}.df, st_se{k}.tstat, p_raw(k), p_holm(k), d_se(k));
end
fclose(fid);
fprintf('Behaviour stats saved -> %s\n', fullfile(pdat.fig, 'Behavior_stats.txt'));

end


%% ===== helpers =====

function p_adj = holm(p)
% Holm-Bonferroni adjusted p-values (step-down), same order as input
m = numel(p); [ps, order] = sort(p(:)');
adj = min(1, cummax((m - (1:m) + 1) .* ps));
p_adj = nan(size(p)); p_adj(order) = adj;
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

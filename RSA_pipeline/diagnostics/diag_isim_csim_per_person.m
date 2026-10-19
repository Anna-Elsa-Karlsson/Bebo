function diag_isim_csim_per_person(resultsDir, variant)
% DIAGNOSTIC ONLY -- does not touch the production pipeline (2026-10-19).
% Per participant: mean item similarity (Isim) and category similarity
% (Csim) in the Ispec ROI, with SE over pairs and the number of pairs, from
% the single-trial tables (one value per pair = mean over the Ispec ROI,
% baseline-corrected Fisher z). Ispec = mean Isim - mean Csim.
%
% Consistency check: the per-person Isim / Csim from the single-trial
% tables must equal the per-person values of the group-level analysis
% (GA/<variant>/<variant>_mean_specA/B_in_ROI.mat, Ispec_derive_alphaROI.m),
% because both are means over the same pairs and the same ROI. Before the
% 2026-10-19 rerun Csim did not match (r = .93) -- after the rerun the
% differences should be ~0.
%
%   resultsDir : folder with GA/, SingleTrialMerged/ (default: the Results
%                folder next to RSA_pipeline/, i.e. ~/BEBO2/Results)
%   variant    : default 'Ispec_30to100Hz'
%
% Reads   SingleTrialMerged/<variant>/<variant>_singleTrial_within-item_fullSample.txt
%         SingleTrialMerged/<variant>/<variant>_singleTrial_within-cat_fullSample.txt
%         GA/<variant>/<variant>_mean_specA_in_ROI.mat, _specB_ (if present)
% Writes  Plots/<variant>/Diag_IsimCsim_perPerson.png/.pdf   (every participant)
%         Plots/<variant>/Diag_IsimCsim_summary.png/.pdf     (4 summary panels)
%         Stats/<variant>/Diag_IsimCsim_perPerson.csv         (all values)
% R2017b-compatible (subplot, errorbar, print).

pipelineRoot = fileparts(fileparts(mfilename('fullpath')));   % RSA_pipeline/
addpath(pipelineRoot);                                          % ids_ALL.m
if nargin < 1 || isempty(resultsDir), resultsDir = fullfile(pipelineRoot, '..', 'Results'); end
if nargin < 2 || isempty(variant), variant = 'Ispec_30to100Hz'; end
stmDir  = fullfile(resultsDir, 'SingleTrialMerged', variant);
plotDir = fullfile(resultsDir, 'Plots', variant);
statDir = fullfile(resultsDir, 'Stats', variant);
for d = {plotDir, statDir}, if ~exist(d{1}, 'dir'), mkdir(d{1}); end, end
col_YA = [30 144 255]/255; col_OA = [80 80 80]/255;

%% Per-pair values -> per person
I = readtable(fullfile(stmDir, [variant '_singleTrial_within-item_fullSample.txt']));
C = readtable(fullfile(stmDir, [variant '_singleTrial_within-cat_fullSample.txt']));
idsAll = [ids_ALL('ya'), ids_ALL('oa')]';
ids = idsAll(ismember(idsAll, unique(I.ID)));
if numel(ids) < numel(idsAll)
    warning('diag:ids', '%d participant(s) of ids_ALL are not in the tables: %s', ...
        numel(idsAll) - numel(ids), mat2str(setdiff(idsAll, ids)'));
end
n = numel(ids);
P = table(ids, ids < 4000, nan(n,1), nan(n,1), nan(n,1), nan(n,1), nan(n,1), nan(n,1), ...
    'VariableNames', {'ID','YA','nIsim','nCsim','Isim','Isim_SE','Csim','Csim_SE'});
for s = 1:n
    a = I.Sim(I.ID == ids(s)); b = C.Sim(C.ID == ids(s));
    P.nIsim(s) = numel(a); P.nCsim(s) = numel(b);
    P.Isim(s) = mean(a); P.Isim_SE(s) = std(a) / sqrt(numel(a));
    P.Csim(s) = mean(b); P.Csim_SE(s) = std(b) / sqrt(numel(b));
end
P.Ispec = P.Isim - P.Csim;
P.z_Ispec = nan(n,1);   % z within age group
for g = [true false]
    k = P.YA == g; P.z_Ispec(k) = (P.Ispec(k) - mean(P.Ispec(k))) / std(P.Ispec(k));
end

%% Consistency with the group-level per-person values
P.Isim_GA = nan(n,1); P.Csim_GA = nan(n,1);
fA = fullfile(resultsDir, 'GA', variant, [variant '_mean_specA_in_ROI.mat']);
fB = fullfile(resultsDir, 'GA', variant, [variant '_mean_specB_in_ROI.mat']);
haveGA = exist(fA, 'file') == 2 && exist(fB, 'file') == 2;
consTxt = 'group-level values not found';
if haveGA
    LA = load(fA, 'specA'); LB = load(fB, 'specB');
    gaIDs = [ids_ALL('ya'), ids_ALL('oa')]';
    vA = double([LA.specA.YA(:); LA.specA.OA(:)]); vB = double([LB.specB.YA(:); LB.specB.OA(:)]);
    if numel(vA) == numel(gaIDs)
        [~, ig] = ismember(ids, gaIDs);
        P.Isim_GA = vA(ig); P.Csim_GA = vB(ig);
        consTxt = sprintf('tables vs group level: Isim r = %.4f (max |diff| %.2g), Csim r = %.4f (max |diff| %.2g)', ...
            corr(P.Isim, P.Isim_GA), max(abs(P.Isim - P.Isim_GA)), corr(P.Csim, P.Csim_GA), max(abs(P.Csim - P.Csim_GA)));
    else
        consTxt = sprintf('group-level files hold %d participants, ids_ALL %d -- not compared', numel(vA), numel(gaIDs));
    end
end
writetable(P, fullfile(statDir, 'Diag_IsimCsim_perPerson.csv'));

%% Report
fprintf('\n=== Isim / Csim per participant (%s, n = %d: %d YA, %d OA) ===\n', variant, n, sum(P.YA), sum(~P.YA));
fprintf('Consistency: %s\n', consTxt);
for g = [true false]
    k = P.YA == g; lab = 'OA'; if g, lab = 'YA'; end
    fprintf('%s: Isim M=%.4f SD=%.4f | Csim M=%.4f SD=%.4f | Ispec M=%.4f SD=%.4f | Isim pairs %d-%d, Csim pairs %d-%d\n', lab, ...
        mean(P.Isim(k)), std(P.Isim(k)), mean(P.Csim(k)), std(P.Csim(k)), mean(P.Ispec(k)), std(P.Ispec(k)), ...
        min(P.nIsim(k)), max(P.nIsim(k)), min(P.nCsim(k)), max(P.nCsim(k)));
end
fprintf('Ispec < 0 (Csim > Isim): %d participant(s): %s\n', sum(P.Ispec < 0), mat2str(P.ID(P.Ispec < 0)'));
fprintf('|z(Ispec)| > 3 within age group: %s\n', mat2str(P.ID(abs(P.z_Ispec) > 3)'));
fprintf('Fewer than 10 Csim pairs: %s\n', mat2str(P.ID(P.nCsim < 10)'));
fprintf('r(number of Csim pairs, |z(Ispec)|) = %.2f\n', corr(P.nCsim, abs(P.z_Ispec)));

%% Figure 1: every participant (sorted by Ispec within age group)
f = figure('Units', 'centimeters', 'Position', [2 2 42 24], 'Color', 'w', 'Visible', 'off');
grp = {'Younger adults', 'Older adults'}; isYA = [true false]; cols = {col_YA, col_OA};
for g = 1:2
    Q = sortrows(P(P.YA == isYA(g), :), 'Ispec');
    x = (1:height(Q))';
    ax = subplot(2, 1, g); hold(ax, 'on');
    plot(ax, [0 height(Q)+1], [0 0], 'k:');
    for s = 1:height(Q), plot(ax, [x(s) x(s)], [Q.Csim(s) Q.Isim(s)], '-', 'Color', [.75 .75 .75]); end
    hI = errorbar(ax, x - 0.15, Q.Isim, Q.Isim_SE, 'o', 'Color', cols{g}, 'MarkerFaceColor', cols{g}, 'MarkerSize', 4, 'CapSize', 0);
    hC = errorbar(ax, x + 0.15, Q.Csim, Q.Csim_SE, 's', 'Color', cols{g}, 'MarkerFaceColor', 'w', 'MarkerSize', 4, 'CapSize', 0);
    yl = ylim(ax);
    for s = 1:height(Q)
        text(ax, x(s), yl(1), sprintf('%d\n%d/%d', Q.ID(s), Q.nIsim(s), Q.nCsim(s)), 'FontSize', 5, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', 90);
    end
    xlim(ax, [0 height(Q)+1]); set(ax, 'XTick', []);
    ylabel(ax, 'Similarity (ROI mean, Fisher z)');
    title(ax, sprintf('%s (n = %d): mean \\pm SE over pairs, sorted by Ispec; labels: ID, Isim/Csim pairs', grp{g}, height(Q)));
    legend(ax, [hI hC], {'Item similarity (Isim)', 'Category similarity (Csim)'}, 'Location', 'northwest', 'Box', 'off');
end
save_fig(f, fullfile(plotDir, 'Diag_IsimCsim_perPerson'));

%% Figure 2: summary
f = figure('Units', 'centimeters', 'Position', [2 2 34 26], 'Color', 'w', 'Visible', 'off');
ax = subplot(2,2,1); hold(ax, 'on');
lims = [min([P.Isim; P.Csim]) max([P.Isim; P.Csim])]; lims = lims + [-1 1] * 0.05 * diff(lims);
plot(ax, lims, lims, 'k:');
scatter(ax, P.Csim(P.YA), P.Isim(P.YA), 25, col_YA, 'filled');
scatter(ax, P.Csim(~P.YA), P.Isim(~P.YA), 25, col_OA, 'filled');
axis(ax, [lims lims]); axis(ax, 'square');
xlabel(ax, 'Csim'); ylabel(ax, 'Isim'); title(ax, 'Isim vs Csim (above the line: Ispec > 0)');
legend(ax, {'', 'Younger', 'Older'}, 'Location', 'southeast', 'Box', 'off');

ax = subplot(2,2,2); hold(ax, 'on');
scatter(ax, P.nCsim(P.YA), P.Ispec(P.YA), 25, col_YA, 'filled');
scatter(ax, P.nCsim(~P.YA), P.Ispec(~P.YA), 25, col_OA, 'filled');
plot(ax, xlim(ax), [0 0], 'k:');
xlabel(ax, 'Number of Csim pairs'); ylabel(ax, 'Ispec (Isim - Csim)');
title(ax, sprintf('Ispec vs number of pairs (r(n, |z|) = %.2f)', corr(P.nCsim, abs(P.z_Ispec))));

ax = subplot(2,2,3); hold(ax, 'on');
for g = 1:2
    k = P.YA == isYA(g); off = 3 * (g - 1);
    for s = find(k)', plot(ax, off + [1 2], [P.Isim(s) P.Csim(s)], '-', 'Color', cols{g}*0.35 + 0.65); end
    plot(ax, off + [1 2], [mean(P.Isim(k)) mean(P.Csim(k))], '-o', 'Color', cols{g}, 'LineWidth', 2.5, 'MarkerFaceColor', cols{g});
end
set(ax, 'XTick', [1 2 4 5], 'XTickLabel', {'Isim', 'Csim', 'Isim', 'Csim'}); xlim(ax, [0.5 5.5]);
ylabel(ax, 'Similarity (Fisher z)'); title(ax, 'Per participant: Isim \rightarrow Csim (left YA, right OA)');

ax = subplot(2,2,4); hold(ax, 'on');
if haveGA && ~all(isnan(P.Csim_GA))
    scatter(ax, P.Csim_GA, P.Csim, 25, [0.85 0.33 0.1], 'filled');
    scatter(ax, P.Isim_GA, P.Isim, 25, [0.2 0.6 0.2], 'filled');
    l = [min([P.Csim_GA; P.Isim_GA; P.Csim; P.Isim]) max([P.Csim_GA; P.Isim_GA; P.Csim; P.Isim])];
    plot(ax, l, l, 'k:'); axis(ax, 'square');
    legend(ax, {'Csim', 'Isim'}, 'Location', 'southeast', 'Box', 'off');
    xlabel(ax, 'Group-level per-person value (GA, ROI)'); ylabel(ax, 'Single-trial tables (mean over pairs)');
else
    axis(ax, 'off');
end
title(ax, {'Consistency check (should lie on the line)', consTxt}, 'FontSize', 7);
save_fig(f, fullfile(plotDir, 'Diag_IsimCsim_summary'));
fprintf('Saved Diag_IsimCsim_perPerson / _summary (png, pdf) -> %s\n      Diag_IsimCsim_perPerson.csv -> %s\n', plotDir, statDir);
end

function save_fig(f, name)
set(f, 'PaperUnits', 'centimeters'); pos = get(f, 'Position');
set(f, 'PaperSize', pos(3:4), 'PaperPosition', [0 0 pos(3:4)]);
print(f, [name '.png'], '-dpng', '-r200');
print(f, [name '.pdf'], '-dpdf', '-painters');
close(f);
end

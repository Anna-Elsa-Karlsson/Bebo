function rsa_plots(contrastType, freqLo, freqHi, pdat, baseDir, cfg_layout)
% Reproduces the Fig3 (sim matrices) and Fig4 (spec t-val matrix + topo +
% raincloud + correlations) structure from BEBO_analyses.m, parameterized
% by contrast type and frequency band so it runs once per variant instead
% of being hand-copied per analysis.
%
%   contrastType  : 'Ispec' or 'Catspec'
%   freqLo,freqHi : analysis sub-band in Hz
%   pdat          : struct with fields .ga, .stat, .fig, .stm (paths;
%                   .stm = SingleTrialMerged output, for the alpha values)
%   baseDir       : BEBO project root (for Performance_All.txt, Apow/)
%   cfg_layout    : path to the EEG layout .mat file (for topoplots)
%
% Run only after rsa_group_level_fullSample.m has completed for this
% variant.

switch contrastType
    case 'Ispec'
        nameA = 'within-item'; nameB = 'within-cat';
        labA = 'Item similarity'; labB = 'Category similarity';
        specLabel = 'Item specificity'; specYlabel = 'Item - category similarity (r)';
    case 'Catspec'
        nameA = 'within-cat'; nameB = 'between-cat';
        labA = 'Within-category similarity'; labB = 'Between-category similarity';
        specLabel = 'Category specificity'; specYlabel = 'Within - between category similarity (r)';
    otherwise
        error('rsa_plots:contrastType', 'Unknown contrastType "%s"', contrastType)
end

bandLabel = sprintf('%dto%dHz', freqLo, freqHi);
variantLabel = sprintf('%s_%s', contrastType, bandLabel);
gaDir = fullfile(pdat.ga, variantLabel);
statDir = fullfile(pdat.stat, variantLabel);
plotDir = fullfile(pdat.fig, variantLabel);
if ~exist(plotDir,'dir'), mkdir(plotDir); end

ft_available = exist('ft_topoplotER','file') == 2;

%% === Setup: colors, colormaps (matches BEBO_analyses.m exactly) ===

col_YA = [30 144 255]/255;
col_OA = [80 80 80]/255;
cmap_div = bwr_colormap();
cmap_sim = flipud(autumn(256));

ID_ya = ids_ALL('ya'); nYA = length(ID_ya);
ID_oa = ids_ALL('oa'); nOA = length(ID_oa);
nAll = nYA + nOA;
idx_YA = 1:nYA; idx_OA = nYA+1:nAll;
age_vec = [ones(nYA,1); 2*ones(nOA,1)];

%% === Load variant data ===

L = load(fullfile(gaDir, ['YA_' nameA '_sim_GA.mat']),'sim'); sim_YA_A = L.sim;
L = load(fullfile(gaDir, ['YA_' nameB '_sim_GA.mat']),'sim'); sim_YA_B = L.sim;
L = load(fullfile(gaDir, ['OA_' nameA '_sim_GA.mat']),'sim'); sim_OA_A = L.sim;
L = load(fullfile(gaDir, ['OA_' nameB '_sim_GA.mat']),'sim'); sim_OA_B = L.sim;

L = load(fullfile(statDir, [variantLabel '_2nd_level_cluster_stat_vsBL.mat']),'stat2'); stat2 = L.stat2;
% per-subject values IN THE ROI (whole Ispec cluster; Ispec_derive_alphaROI.m)
L = load(fullfile(gaDir, [variantLabel '_mean_specA_in_ROI.mat']),'specA'); specA = L.specA;
L = load(fullfile(gaDir, [variantLabel '_mean_specB_in_ROI.mat']),'specB'); specB = L.specB;
L = load(fullfile(gaDir, [variantLabel '_mean_spec_in_ROI.mat']),'spec'); spec = L.spec;

spec_z = double([spec.YA; spec.OA]);
specA_all = double([specA.YA; specA.OA]);

sim_time = sim_YA_A.time;
% Significant POSITIVE cluster(s) only (as in fullSample / single-trial /
% alpha ROI); stat2.mask would also include negative clusters.
ispec_mask = false(size(stat2.stat));
if isfield(stat2,'posclusters') && ~isempty(stat2.posclusters) % fields are absent when no cluster formed
    for sc = find([stat2.posclusters.prob] < stat2.cfg.alpha)
        ispec_mask = ispec_mask | (stat2.posclusterslabelmat == sc);
    end
end
ispec_chan_active = squeeze(any(reshape(ispec_mask, size(ispec_mask,1), []), 2));

%% === Behavioral + alpha data (for correlation panels) ===
% Both matched to ids_ALL by ID (not row position), and error out if any
% subject is missing, so a stale/short file can't silently misalign.

allIDs = [ID_ya(:); ID_oa(:)];

perf = readtable(fullfile(baseDir,'Performance_All.txt')); % BEBO2/Performance_All.txt (same file as behavior_analysis.m)
[found, iPerf] = ismember(allIDs, perf.ID);
if ~all(found)
    error('rsa_plots:perf', 'Performance_All.txt is missing subject(s): %s', mat2str(allIDs(~found)'));
end
IM = perf.IM(iPerf);
LDI = perf.LDI(iPerf);

% Alpha (2026-09-26, replaces the circular SME/Apow mask): per-subject
% log10 alpha power (alpha_band.m, 8-25 Hz) in the Ispec-derived ROI (cluster channels
% x 0.4-1.7 s, alpha_time_window.m), averaged over exactly the trials that enter
% the Ispec analysis -- computed by rsa_singleTrial_merged.m, collected
% by aggregate_singleTrial_merged.m. Always the Ispec variant's alpha
% (the ROI is Ispec-derived), also when plotting Catspec.
La = load(fullfile(pdat.stm, sprintf('Ispec_%dto%dHz', freqLo, freqHi), ...
    sprintf('Ispec_%dto%dHz_alphaIspec_fullSample.mat', freqLo, freqHi)), 'ID', 'alphaIspec');
[found, iAlpha] = ismember(allIDs, La.ID);
if ~all(found)
    error('rsa_plots:alpha', 'alphaIspec file is missing subject(s): %s', mat2str(allIDs(~found)'));
end
alpha_avg = La.alphaIspec(iAlpha);

% BASELINE-CORRECTED alpha in the same Ispec ROI (2026-09-30, added for
% comparison -- the raw log-power alpha above is kept unchanged). Source:
% the per-subject alpha_pre2post files (log10 wavelet power, per-trial
% baseline correction 'absolute' on log power, -0.6 to -0.2 s = log ratio
% to baseline; averaged over presentation 1-2 'old'/'similar' trials of
% all exemplars, i.e. ~the Ispec trial set), cut to the Ispec ROI:
% cluster channels (Ispec_alphaROI_config.mat) x 0.4-1.7 s (alpha_time_window.m) x
% alphaBandBL. Negative = desync.
alphaBandBL = alpha_band(); % same band as everywhere (alpha_band.m)
Lroi = load(fullfile(pdat.stm, 'Ispec_alphaROI_config.mat'), 'roiChannels');
Lroi.timeWindow = alpha_time_window(); % 0.4-1.7 s (alpha_time_window.m, 2026-10-06), as the single-trial alpha
alpha_bl = nan(nAll,1);
for s = 1:nAll
    Lp = load(fullfile(pdat.alpha, sprintf('%d_alphaPow.mat', allIDs(s))), 'pow');
    P  = Lp.pow;
    ch = ismember(P.label, Lroi.roiChannels);
    fr = P.freq >= alphaBandBL(1) & P.freq <= alphaBandBL(2);
    tm = P.time >= Lroi.timeWindow(1) & P.time <= Lroi.timeWindow(2);
    if ~any(ch) || ~any(fr) || ~any(tm)
        error('rsa_plots:alphaBL', 'Subject %d: empty ROI selection (chan %d, freq %d, time %d)', allIDs(s), sum(ch), sum(fr), sum(tm));
    end
    tmp = P.powspctrm(ch, fr, tm);   % ROI = channels of the Ispec cluster
    alpha_bl(s) = mean(tmp(:), 'omitnan');
end

%% === Display matrices ===

% DISPLAY-ONLY smoothing (2026-10-02): the time x time matrices show
% row/column stripes (each row = one time point of trial A correlated with
% all of trial B; 15 frequencies, 8 ms steps). For the figures only, they
% are smoothed with a 2-D Gaussian (sigma = smoothSigma s, edge-corrected).
% Statistics, cluster outlines (contours from the UNSMOOTHED mask) and all
% values are unaffected. State in the figure captions.
smoothSigma = display_smooth_sigma();                  % s (display_smooth_sigma.m; 32 ms = 4 samples at 8 ms)
sigmaSmp = smoothSigma / mean(diff(sim_YA_A.time));    % in samples

sim_display = @(s) smooth2d(squeeze(mean(mean(double(s.powspctrm), 1), 2)), sigmaSmp); % R2017b: mean() only takes one dim at a time, not [1 2]
GA_YA_A_disp = sim_display(sim_YA_A); GA_YA_B_disp = sim_display(sim_YA_B);
GA_OA_A_disp = sim_display(sim_OA_A); GA_OA_B_disp = sim_display(sim_OA_B);

clim_sim3 = sym_clim_pos([GA_YA_A_disp(:); GA_YA_B_disp(:); GA_OA_A_disp(:); GA_OA_B_disp(:)]);

%% === "Fig3"-equivalent: sim matrices ===

fig3 = figure('Units','centimeters','Position',[2 2 28 22],'Color','w');
% R2017b: tiledlayout/nexttile (R2019b+) not available -- using subplot.
% Layout (2026-10-02, diagonal panels dropped -- not used in the
% manuscript): rows 1-2 = similarity matrices (2 columns), row 3 =
% raincloud (2/3 width) + interaction (1/3), on a 3x3 grid.

plot_simmat(subplot(3,2,1), GA_YA_A_disp, sim_time, clim_sim3, cmap_sim, ['Younger adults — ' labA]);
plot_simmat(subplot(3,2,2), GA_YA_B_disp, sim_time, clim_sim3, cmap_sim, ['Younger adults — ' labB]);
plot_simmat(subplot(3,2,3), GA_OA_A_disp, sim_time, clim_sim3, cmap_sim, ['Older adults — ' labA]);
plot_simmat(subplot(3,2,4), GA_OA_B_disp, sim_time, clim_sim3, cmap_sim, ['Older adults — ' labB]);

ax = subplot(3,3,[7 8]);
raincloud_4group(specA.YA, specA.OA, specB.YA, specB.OA, col_YA, col_OA, col_YA*0.4+0.6, [170 170 170]/255, ax, {labA, labB});
ylabel(ax,'Similarity (r)');

% Age x similarity-type interaction, next to the raincloud (as in the
% original Fig3 in BEBO_analyses.m -- was missing from this port)
ax = subplot(3,3,9);
plot_interaction_sim(ax, specA.YA, specA.OA, specB.YA, specB.OA, [0 0 0], [120 70 215]/255, labA, labB, 'Similarity (r)');
title(ax, 'Age x similarity-type interaction');

savefig_all(fig3, fullfile(plotDir, ['Fig3_' variantLabel '_SimMatrices']));

%% === "Fig4"-equivalent: spec t-val matrix + topo + raincloud + correlations ===

spec_tval_display = smooth2d(squeeze(mean(double(stat2.stat), 1)), sigmaSmp); % display-only smoothing (see above)
clim_tval = sym_clim(spec_tval_display);

fig4 = figure('Units','centimeters','Position',[2 2 28 24],'Color','w');
% R2017b: subplot instead of tiledlayout/nexttile, same mapping approach as Fig3 above.

ax = subplot(3,3,[1 2 4 5]);
imagesc(ax, sim_time, sim_time, spec_tval_display); set(ax,'YDir','normal');
colormap(ax,cmap_div); caxis(ax,clim_tval); axis(ax,'square'); hold(ax,'on');
mask2d = squeeze(any(ispec_mask,1));
contour(ax, sim_time, sim_time, double(mask2d), [0.5 0.5], 'k', 'LineWidth', 1.5);
xlabel(ax,'Encoding 2nd (time, s)'); ylabel(ax,'Encoding 1st (time, s)');
title(ax, specLabel); cb = colorbar(ax); ylabel(cb,'t-values');

ax = subplot(3,3,3);
if ft_available
    % Colour: mean t per channel over the cluster's extent only (time x
    % time points where the cluster exists on any channel) -- same rule as
    % the alpha topoplot.
    stat_ct = reshape(double(stat2.stat), size(stat2.stat,1), []);
    spec_tval_perchan = stat_ct(:, mask2d(:)) * ones(sum(mask2d(:)),1) / sum(mask2d(:));
    % Asterisks: the ROI channels = all channels of the significant
    % cluster (Ispec_derive_alphaROI.m)
    coreCh = ismember(sim_YA_A.label, Lroi.roiChannels);
    clim_topo = max(abs(spec_tval_perchan)) * [-1 1];
    caxis(ax, clim_topo); colormap(ax, cmap_div);
    cb_topo = colorbar(ax,'Location','eastoutside'); cb_topo.Position(3) = cb_topo.Position(3)*0.5;
    ylabel(cb_topo,'t-values');
    dummy = struct(); dummy.label = sim_YA_A.label; dummy.avg = spec_tval_perchan; dummy.time = 0; dummy.dimord = 'chan_time';
    cfg_topo = []; cfg_topo.layout = cfg_layout; cfg_topo.marker = 'off'; cfg_topo.style = 'fill';
    cfg_topo.colormap = cmap_div; cfg_topo.axes = ax; cfg_topo.comment = 'no'; cfg_topo.zlim = clim_topo;
    cfg_topo.highlight = 'on'; cfg_topo.highlightchannel = sim_YA_A.label(coreCh);
    cfg_topo.highlightsymbol = '*'; cfg_topo.highlightsize = 8; cfg_topo.highlightcolor = [0 0 0];
    ft_topoplotER(cfg_topo, dummy);
    title(ax,{'Topographic distribution (mean t over cluster)', sprintf('* cluster channels = ROI (%d channels)', sum(coreCh))},'FontSize',9);
else
    axis(ax,'off'); text(ax,0.5,0.5,'FieldTrip required','HorizontalAlignment','center');
end

ax = subplot(3,3,[6 9]);
raincloud_compare(spec.YA, spec.OA, col_YA, col_OA, ax);
title(ax, specLabel); ylabel(ax, specYlabel);
set(ax,'XTick',[1 2],'XTickLabel',{'Younger','Older'});

[r_avg,p_avg] = partialcorr(alpha_avg, spec_z, age_vec);
ax = subplot(3,3,7);
corr_scatter(alpha_avg, spec_z, age_vec, col_YA, col_OA, ax, 'Alpha power (log10, Ispec ROI)', specYlabel, sprintf('Partial r=%.2f, p=%.3f', r_avg, p_avg));

[r_im,p_im] = partialcorr(IM, spec_z, age_vec);
ax = subplot(3,3,8);
corr_scatter(IM, spec_z, age_vec, col_YA, col_OA, ax, 'Item memory performance', specYlabel, sprintf('Partial r=%.2f, p=%.3f', r_im, p_im));

savefig_all(fig4, fullfile(plotDir, ['Fig4_' variantLabel '_' strrep(specLabel,' ','')]));

%% === Fig4b: Ispec x alpha, raw log power vs baseline-corrected (comparison) ===
[r_bl, p_bl] = partialcorr(alpha_bl, spec_z, age_vec);
fig4b = figure('Units','centimeters','Position',[2 2 24 11],'Color','w');
ax = subplot(1,2,1);
corr_scatter(alpha_avg, spec_z, age_vec, col_YA, col_OA, ax, 'Alpha power (log10, Ispec ROI)', specYlabel, sprintf('Raw log power: partial r=%.2f, p=%.3f', r_avg, p_avg));
ax = subplot(1,2,2);
corr_scatter(alpha_bl, spec_z, age_vec, col_YA, col_OA, ax, 'Alpha power change vs. baseline (log ratio, Ispec ROI)', specYlabel, sprintf('Baseline-corrected: partial r=%.2f, p=%.3f', r_bl, p_bl));
savefig_all(fig4b, fullfile(plotDir, ['Fig4b_' variantLabel '_alpha_comparison']));


%% === Stats text file ===

L = load(fullfile(statDir, [variantLabel '_spec_in_ROI_3rd_level_stats_Pairwise.mat']),'YA_vs_OA','STATS');
YA_vs_OA = L.YA_vs_OA; STATS = L.STATS;

fid = fopen(fullfile(plotDir, [variantLabel '_stats.txt']), 'w');
fprintf(fid, '%s Statistical Results (%s, %d-%dHz)\nGenerated: %s\n\n', contrastType, contrastType, freqLo, freqHi, datestr(now));
fprintf(fid, '--- Younger vs Older %s (per-subject values in the significant cluster) ---\n', lower(specLabel));
fprintf(fid, 'Younger: M=%.4f, SD=%.4f\n', YA_vs_OA.M_YA, YA_vs_OA.SD_YA);
fprintf(fid, 'Older:   M=%.4f, SD=%.4f\n', YA_vs_OA.M_OA, YA_vs_OA.SD_OA);
fprintf(fid, 't(%.0f)=%.3f, p=%.4f, Cohen d=%.3f\n\n', STATS.df, STATS.tstat, YA_vs_OA.P, YA_vs_OA.cohensD);
% Correlation set (as section 8 of BEBO_analyses.m): partial r controlling
% age (full sample), within-group r (YA, OA), Fisher z-test YA vs OA; FDR
% (Benjamini-Hochberg) separately within the partial-r family and the
% Fisher-z family (3 tests each). Only alpha and IM are plotted (Fig4).
preds  = {alpha_avg, IM, LDI};
plabel = {'alpha (log10, Ispec ROI, all Ispec trials)', 'IM', 'LDI'};
nP = numel(preds);
r_part = nan(1,nP); p_part = nan(1,nP); r_YA = nan(1,nP); r_OA = nan(1,nP); z_f = nan(1,nP); p_f = nan(1,nP);
for k = 1:nP
    x = preds{k};
    [r_part(k), p_part(k)] = partialcorr(x, spec_z, age_vec);
    r_YA(k) = corr(x(idx_YA), spec_z(idx_YA));
    r_OA(k) = corr(x(idx_OA), spec_z(idx_OA));
    z_f(k) = (atanh(r_YA(k)) - atanh(r_OA(k))) / sqrt(1/(nYA-3) + 1/(nOA-3));
    p_f(k) = 2 * (1 - normcdf(abs(z_f(k))));
end
p_part_fdr = fdr_bh(p_part);
p_f_fdr    = fdr_bh(p_f);

fprintf(fid, '--- Correlations with %s (N=%d: %d YA, %d OA) ---\n', specLabel, nAll, nYA, nOA);
fprintf(fid, 'Partial r controlling age (full sample); within-group Pearson r; Fisher z-test YA vs OA.\n');
fprintf(fid, 'FDR (Benjamini-Hochberg) applied separately within the partial-r family and the Fisher-z family (%d tests each).\n\n', nP);
for k = 1:nP
    fprintf(fid, '%s ~ %s\n', specLabel, plabel{k});
    fprintf(fid, '   partial r=%.3f, p=%.4f, p_FDR=%.4f\n', r_part(k), p_part(k), p_part_fdr(k));
    fprintf(fid, '   YA r=%.3f (n=%d), OA r=%.3f (n=%d)\n', r_YA(k), nYA, r_OA(k), nOA);
    fprintf(fid, '   Fisher z=%.2f, p=%.4f, p_FDR=%.4f\n\n', z_f(k), p_f(k), p_f_fdr(k));
end
% --- Comparison: baseline-corrected alpha (same Ispec ROI) ---
% FDR computed as if BL alpha REPLACED the raw log-power alpha in the
% family above (BL alpha, IM, LDI), i.e. the numbers to report if we switch.
r_blYA = corr(alpha_bl(idx_YA), spec_z(idx_YA));
r_blOA = corr(alpha_bl(idx_OA), spec_z(idx_OA));
z_bl = (atanh(r_blYA) - atanh(r_blOA)) / sqrt(1/(nYA-3) + 1/(nOA-3));
p_zbl = 2 * (1 - normcdf(abs(z_bl)));
p_part_alt = fdr_bh([p_bl, p_part(2:3)]);
p_f_alt    = fdr_bh([p_zbl, p_f(2:3)]);
fprintf(fid, '--- COMPARISON: %s ~ alpha, BASELINE-CORRECTED (log ratio to -0.6..-0.2 s; same Ispec ROI, %g-%g Hz) ---\n', specLabel, alphaBandBL(1), alphaBandBL(2));
fprintf(fid, 'Source: alpha_pre2post per-subject files (presentation 1-2, old/similar, all exemplars). Negative = desync.\n');
fprintf(fid, 'p_FDR computed with BL alpha replacing raw log alpha in the family (BL alpha, IM, LDI).\n');
fprintf(fid, '   partial r=%.3f, p=%.4f, p_FDR=%.4f\n', r_bl, p_bl, p_part_alt(1));
fprintf(fid, '   YA r=%.3f (n=%d), OA r=%.3f (n=%d)\n', r_blYA, nYA, r_blOA, nOA);
fprintf(fid, '   Fisher z=%.2f, p=%.4f, p_FDR=%.4f\n\n', z_bl, p_zbl, p_f_alt(1));

% --- CONTROL: does itemRep explain the age difference in Ispec? (2026-09-30) ---
% Per subject, Ispec split by itemRep condition, from the aggregated
% single-trial tables (per-pair cluster-mean similarity 'Sim', itemRep
% -1 = low, 1 = high): Ispec_r = mean(Isim | r) - mean(Csim | r).
% Averaging over pairs and over the cluster commute, so this is the same
% quantity as the fullSample Ispec, split by condition.
if strcmp(contrastType, 'Ispec')
    stDir = fullfile(pdat.stm, variantLabel);
    Li = load(fullfile(stDir, [variantLabel '_singleTrial_within-item_fullSample.mat']), 'T_item');
    Lc = load(fullfile(stDir, [variantLabel '_singleTrial_within-cat_fullSample.mat']),  'T_cat');
    Ti = Li.T_item; Tc = Lc.T_cat;
    ispec_rep = nan(nAll, 2); % columns: low (-1), high (1) itemRep
    repLev = [-1 1];
    for s = 1:nAll
        for r = 1:2
            mi = mean(Ti.Sim(Ti.ID == allIDs(s) & Ti.itemRep == repLev(r)), 'omitnan');
            mc = mean(Tc.Sim(Tc.ID == allIDs(s) & Tc.itemRep == repLev(r)), 'omitnan');
            ispec_rep(s,r) = mi - mc;
        end
    end
    dRep = ispec_rep(:,2) - ispec_rep(:,1); % high - low

    fprintf(fid, '--- CONTROL: itemRep and the age difference in %s ---\n', specLabel);
    fprintf(fid, 'Per-subject %s split by itemRep (from single-trial tables): low (normal, high-exemplar) vs high (high-repetition, high-repetition+exemplar).\n', specLabel);
    fprintf(fid, 'YA low: M=%.4f SD=%.4f | YA high: M=%.4f SD=%.4f\n', mean(ispec_rep(idx_YA,1),'omitnan'), std(ispec_rep(idx_YA,1),'omitnan'), mean(ispec_rep(idx_YA,2),'omitnan'), std(ispec_rep(idx_YA,2),'omitnan'));
    fprintf(fid, 'OA low: M=%.4f SD=%.4f | OA high: M=%.4f SD=%.4f\n\n', mean(ispec_rep(idx_OA,1),'omitnan'), std(ispec_rep(idx_OA,1),'omitnan'), mean(ispec_rep(idx_OA,2),'omitnan'), std(ispec_rep(idx_OA,2),'omitnan'));

    % (1) Age x itemRep: Welch t YA vs OA on (high - low)
    [~, p_int, ~, st_int] = ttest2(dRep(idx_YA), dRep(idx_OA), 'Vartype','unequal');
    fprintf(fid, '(1) Age x itemRep interaction (Welch t, YA vs OA on high - low): t(%.1f)=%.3f, p=%.4f\n', st_int.df, st_int.tstat, p_int);

    % (2) Age difference within each itemRep level (Welch), Holm over 2
    p_lev = nan(1,2); st_lev = cell(1,2); d_lev = nan(1,2); lvName = {'low itemRep ', 'high itemRep'};
    for r = 1:2
        a = ispec_rep(idx_YA,r); b = ispec_rep(idx_OA,r); a = a(~isnan(a)); b = b(~isnan(b));
        [~, p_lev(r), ~, st_lev{r}] = ttest2(a, b, 'Vartype','unequal');
        d_lev(r) = (mean(a) - mean(b)) / sqrt(((numel(a)-1)*var(a) + (numel(b)-1)*var(b)) / (numel(a)+numel(b)-2));
    end
    p_lev_holm = holm(p_lev);
    for r = 1:2
        fprintf(fid, '(2) YA vs OA, %s: Welch t(%.1f)=%.3f, p=%.4f, p_Holm=%.4f, d=%.3f\n', lvName{r}, st_lev{r}.df, st_lev{r}.tstat, p_lev(r), p_lev_holm(r), d_lev(r));
    end

    % (3) Main effect of itemRep (paired, all subjects)
    ok = ~isnan(dRep);
    [~, p_main, ~, st_main] = ttest(ispec_rep(ok,2), ispec_rep(ok,1));
    fprintf(fid, '(3) itemRep main effect (paired t, high vs low, n=%d): t(%d)=%.3f, p=%.4f, dz=%.3f\n\n', sum(ok), st_main.df, st_main.tstat, p_main, mean(dRep(ok))/std(dRep(ok)));
end
fclose(fid);

fprintf('Plots + stats saved for %s -> %s\n', variantLabel, plotDir);

end


%% =========================================================
%% LOCAL HELPER FUNCTIONS (matching BEBO_analyses.m conventions exactly)
%% =========================================================

function cl = sym_clim(data)
v = max(abs(prctile(data(:), [2 98])));
cl = [-v v];
end

function cl = sym_clim_pos(data)
% positive-only similarity data: clip low end near typical similarity floor
lo = prctile(data(:), 5); hi = prctile(data(:), 98);
cl = [max(lo,0) hi];
end

function raincloud_4group(ya_h, oa_h, ya_m, oa_m, c_yah, c_oah, c_yam, c_oam, ax, labs)
hold(ax,'on');
raincloud_single(ya_h, 1.0, c_yah, ax); raincloud_single(ya_m, 1.7, c_yam, ax);
raincloud_single(oa_h, 3.0, c_oah, ax); raincloud_single(oa_m, 3.7, c_oam, ax);
xlim(ax,[0.3 4.4]); box(ax,'off');
plot(ax, [2.35 2.35], ylim(ax), '--', 'Color',[0.75 0.75 0.75], 'LineWidth', 0.8); % R2017b: xline (R2018b+) not available
set(ax,'XTick',[1.35 3.35],'XTickLabel',{'Younger','Older'});
h_hit = patch(ax,NaN,NaN,c_yah,'EdgeColor','none'); h_miss = patch(ax,NaN,NaN,c_yam,'EdgeColor','none');
legend(ax,[h_hit h_miss], labs,'Location','best','Box','off','FontSize',8);
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

function out = smooth2d(mat, sigmaSmp)
% 2-D Gaussian smoothing for DISPLAY only; normalised by the kernel mass
% inside the matrix (no darkening at the edges); NaNs ignored.
if sigmaSmp <= 0, out = mat; return; end
h = ceil(3*sigmaSmp);
g = exp(-(-h:h).^2 / (2*sigmaSmp^2)); g = g / sum(g);
k = g' * g;
ok = ~isnan(mat); m = mat; m(~ok) = 0;
out = conv2(m, k, 'same') ./ conv2(double(ok), k, 'same');
out(~ok) = NaN;
end

function plot_simmat(ax, mat, t, clim_sim, cmap, ttl)
imagesc(ax, t, t, mat); set(ax,'YDir','normal'); colormap(ax, cmap); caxis(ax, clim_sim);
xlabel(ax,'Encoding 2nd (time, s)'); ylabel(ax,'Encoding 1st (time, s)'); title(ax, ttl); colorbar(ax);
end

function plot_interaction_sim(ax, ya1, oa1, ya2, oa2, col1, col2, lab1, lab2, ylab)
% Interaction plot: mean +/- SEM per age group, one line per similarity
% type (ported unchanged from BEBO_analyses.m)
ya1 = ya1(~isnan(ya1)); oa1 = oa1(~isnan(oa1)); ya2 = ya2(~isnan(ya2)); oa2 = oa2(~isnan(oa2));
m1 = [mean(ya1) mean(oa1)]; se1 = [std(ya1)/sqrt(length(ya1)) std(oa1)/sqrt(length(oa1))];
m2 = [mean(ya2) mean(oa2)]; se2 = [std(ya2)/sqrt(length(ya2)) std(oa2)/sqrt(length(oa2))];
hold(ax,'on');
errorbar(ax, [1 2], m1, se1, '-o', 'Color', col1, 'LineWidth', 1.5, ...
    'MarkerFaceColor', col1, 'MarkerSize', 5, 'CapSize', 6);
errorbar(ax, [1 2], m2, se2, '-o', 'Color', col2, 'LineWidth', 1.5, ...
    'MarkerFaceColor', col2, 'MarkerSize', 5, 'CapSize', 6);
xlim(ax,[0.7 2.3]); box(ax,'off');
set(ax,'XTick',[1 2],'XTickLabel',{'Younger','Older'});
ylabel(ax, ylab);
legend(ax,{lab1,lab2},'Location','best','Box','off','FontSize',8);
end

function corr_scatter(x, y, age_vec, col_YA, col_OA, ax, xlab, ylab, title_str)
idx_YA = age_vec == 1; idx_OA = age_vec == 2;
hold(ax,'on');
scatter(ax, x(idx_OA), y(idx_OA), 30, col_OA, 'x', 'LineWidth', 1.2);
scatter(ax, x(idx_YA), y(idx_YA), 30, col_YA, 'o', 'filled', 'MarkerFaceAlpha', 0.6);
add_regline(ax, x(idx_YA), y(idx_YA), col_YA); add_regline(ax, x(idx_OA), y(idx_OA), col_OA);
xlabel(ax, xlab); ylabel(ax, ylab); title(ax, title_str, 'FontSize', 9); box(ax,'off');
legend(ax, {'Older','Younger'}, 'Location','best', 'Box','off', 'FontSize',8);
end

function add_regline(ax, x, y, col)
p = polyfit(x(~isnan(x)&~isnan(y)), y(~isnan(x)&~isnan(y)), 1);
xl = [min(x)*1.05, max(x)*1.05];
plot(ax, xl, polyval(p,xl), '-', 'Color', col, 'LineWidth', 1.2);
end

function savefig_all(fig, filepath)
% exportgraphics (R2020a+) is not available on the grid's MATLAB R2017b --
% use print() instead (2026-08-27 fix, same root cause as the earlier
% sum(...,'all',...) R2017b-compatibility bug).
set(fig, 'Units','centimeters'); pos = get(fig,'Position');
set(fig, 'PaperUnits','centimeters', 'PaperSize', pos(3:4), 'PaperPosition', [0 0 pos(3:4)]); % PDF page = figure size (otherwise cut off)
print(fig, [filepath '.pdf'], '-dpdf');
print(fig, [filepath '.png'], '-dpng', '-r300');
print(fig, [filepath '.tif'], '-dtiff', '-r300');
fprintf('Saved: %s (.pdf/.png/.tif)\n', filepath);
end

function p_adj = holm(p)
% Holm-Bonferroni adjusted p-values (step-down), same order as input
m = numel(p); [ps, order] = sort(p(:)');
adj = min(1, cummax((m - (1:m) + 1) .* ps));
p_adj = nan(size(p)); p_adj(order) = adj;
end

function p_adj = fdr_bh(p)
% Benjamini-Hochberg adjusted p-values (step-up), same order as input
m = numel(p); [ps, order] = sort(p(:)');
adj = ps .* m ./ (1:m);
adj = min(1, fliplr(cummin(fliplr(adj))));
p_adj = nan(size(p)); p_adj(order) = adj;
end

function cmap = bwr_colormap(n)
if nargin < 1, n = 256; end
half = floor(n/2);
blue_to_white = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
white_to_red = [ones(n-half,1), linspace(1,0,n-half)', linspace(1,0,n-half)'];
cmap = [blue_to_white; white_to_red];
end

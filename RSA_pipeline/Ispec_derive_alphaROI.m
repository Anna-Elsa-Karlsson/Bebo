function Ispec_derive_alphaROI(freqLo, freqHi, pdat)
% Derives the alpha-power extraction ROI (channels + time window) from
% the Ispec 2nd-level cluster stat -- replaces the old, circular SME
% (Item-vs-Cat)-derived alpha ROI. Uses the SAME cluster (and the same
% POSITIVE-cluster-only convention) that rsa_singleTrial_merged.m already
% uses to build its RSA similarity mask, so the alpha-power ROI is tied
% to the exact same non-circular, already-validated effect.
%
% The Ispec cluster's two axes are BOTH activation time (time1 x time2,
% an RSA correlation matrix -- sim_mat.freq = sim_mat.time convention, so
% stat2.freq holds time1's values and stat2.time holds time2's values).
% There is no frequency axis to inherit (Ispec is a 30-100Hz gamma
% effect; alpha power is a different construct) -- this derives ONLY the
% spatial (channels) and temporal (single time window, from the union of
% time1/time2) extent. The actual alpha BAND (alpha_band.m) is fixed
% separately, from the alpha_pre2post activation-vs-zero cluster, and is
% applied directly in rsa_singleTrial_merged.m.
%
%   freqLo,freqHi : Ispec variant band, e.g. 30,100
%   pdat          : struct with fields .stat (Ispec Stats root), .ga (GA root) and .out
%                   (where to save Ispec_alphaROI_config.mat)

variantLabel = sprintf('Ispec_%dto%dHz', freqLo, freqHi);
L = load(fullfile(pdat.stat, variantLabel, [variantLabel '_2nd_level_cluster_stat_vsBL.mat']), 'stat2');
stat2 = L.stat2;

p_crit = stat2.cfg.alpha;
pos_signif_clust = [];
if isfield(stat2,'posclusters') && ~isempty(stat2.posclusters) % fields are absent when no cluster formed
    pos_cluster_pvals = [stat2.posclusters(:).prob];
    pos_signif_clust = find(pos_cluster_pvals < p_crit);
end
if isempty(pos_signif_clust)
    error('Ispec_derive_alphaROI:noCluster', 'No significant positive cluster found in %s''s 2nd-level stat.', variantLabel);
end

mask = false(size(stat2.stat));
for sc = pos_signif_clust
    mask = mask | (stat2.posclusterslabelmat == sc);
end

% ROI = the WHOLE significant cluster (2026-10-02; the core-channel
% selection -- above-average cluster mass -- was DROPPED: it was an extra,
% post-hoc selection that further inflates the Ispec values and is a
% researcher degree of freedom). roiChannels = ALL channels that are part
% of the significant cluster(s); no ranking or per-channel profile.
touchedChan = any(reshape(mask, size(mask,1), []), 2);
roiChannels = stat2.label(touchedChan);

% Time window spans the cluster's extent on BOTH axes (kept deliberately,
% 2026-09-26): time1 (stat2.freq) = presentation 1 / exemplar 1 trial,
% time2 (stat2.time) = presentation 2 / exemplar 2 trial (see
% correlate_pairs_and_package: corr(tf1,tf2) -> rows = pair member 1).
% The per-subject alpha measure averages over the trials from BOTH pair
% members, so the window covers both.
anyMask = squeeze(any(mask(touchedChan,:,:),1)); % time1 x time2
[t1idx, t2idx] = find(anyMask);
allTouchedTimeVals = [stat2.freq(t1idx), stat2.time(t2idx)]; % freq field actually holds time1 values
timeWindow = [min(allTouchedTimeVals), max(allTouchedTimeVals)];

fprintf('Derived Ispec ROI (%s): whole cluster, %d channels, time window=[%.2f %.2f]s\n', ...
    variantLabel, length(roiChannels), timeWindow(1), timeWindow(2));
fprintf('ROI channels: %s\n', strjoin(roiChannels(:)', ' '));

if ~exist(pdat.out,'dir'), mkdir(pdat.out); end
save(fullfile(pdat.out, 'Ispec_alphaROI_config.mat'), 'roiChannels', 'timeWindow');
writetable(table(roiChannels(:), 'VariableNames', {'channel'}), fullfile(pdat.out, 'Ispec_cluster_channels.csv'));

%% Per-subject Ispec values IN THE ROI (2026-10-02): every analysis after
%  the group-level stats uses the ROI = all points of the significant
%  positive cluster(s) (identical to the whole-cluster means of
%  rsa_group_level_fullSample.m).
%  Replaces the whole-cluster means from rsa_group_level_fullSample.m (those
%  files are left untouched; plots/eBOSC now read the *_in_ROI files).
%  Per subject: mean over all ROI points (channel x time1 x time2) of the
%  grand-average similarity, within-item (A) and within-cat (B); spec = A-B.
roiMask = mask & repmat(touchedChan, [1 size(mask,2) size(mask,3)]);
gaDir   = fullfile(pdat.ga, variantLabel);
statDir = fullfile(pdat.stat, variantLabel);
specA = struct(); specB = struct(); spec = struct();
for g = {'YA','OA'}
    La = load(fullfile(gaDir, [g{1} '_within-item_sim_GA.mat']), 'sim');
    Lb = load(fullfile(gaDir, [g{1} '_within-cat_sim_GA.mat']), 'sim');
    Pa = reshape(double(La.sim.powspctrm), size(La.sim.powspctrm,1), []);
    Pb = reshape(double(Lb.sim.powspctrm), size(Lb.sim.powspctrm,1), []);
    specA.(g{1}) = mean(Pa(:, roiMask(:)), 2, 'omitnan');
    specB.(g{1}) = mean(Pb(:, roiMask(:)), 2, 'omitnan');
    spec.(g{1})  = specA.(g{1}) - specB.(g{1});
    clear La Lb Pa Pb
end
[H,P,CI,STATS] = ttest2(spec.YA, spec.OA); %#ok<ASGLU>
n1 = numel(spec.YA); n2 = numel(spec.OA);
YA_vs_OA = struct('H',H,'P',P,'CI',CI,'STATS',STATS, ...
    'cohensD', (mean(spec.YA)-mean(spec.OA)) / sqrt(((n1-1)*var(spec.YA)+(n2-1)*var(spec.OA))/(n1+n2-2)), ...
    'M_YA',mean(spec.YA),'SD_YA',std(spec.YA),'M_OA',mean(spec.OA),'SD_OA',std(spec.OA));
save(fullfile(gaDir, [variantLabel '_mean_specA_in_ROI.mat']), 'specA');
save(fullfile(gaDir, [variantLabel '_mean_specB_in_ROI.mat']), 'specB');
save(fullfile(gaDir, [variantLabel '_mean_spec_in_ROI.mat']), 'spec');
save(fullfile(statDir, [variantLabel '_spec_in_ROI_3rd_level_stats_Pairwise.mat']), 'YA_vs_OA', 'STATS');
fprintf('ROI Ispec (whole cluster, %d channels): YA M=%.4f SD=%.4f, OA M=%.4f SD=%.4f, t(%d)=%.3f, p=%.4f, d=%.3f\n', ...
    sum(touchedChan), YA_vs_OA.M_YA, YA_vs_OA.SD_YA, YA_vs_OA.M_OA, YA_vs_OA.SD_OA, STATS.df, STATS.tstat, P, YA_vs_OA.cohensD);
fprintf('Saved Ispec alpha ROI config -> %s\n', fullfile(pdat.out, 'Ispec_alphaROI_config.mat'));

end

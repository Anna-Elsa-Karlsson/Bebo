function PF_derive_roi(pdat)
% Derives PF's channel ROI, time window, and frequency search band
% DIRECTLY from the alpha_pre2post activation-vs-zero cluster stat
% (stat_alpha.mat) -- replaces the old SME (Item-vs-Cat)-cluster-derived
% values, which were circular (SME's Item/Category conditions ARE the
% old/similar response later predicted downstream). Desync = NEGATIVE
% cluster (power < baseline) -- see alpha_pre2post design notes.
%
% Channel selection (2026-09-26): ALL channels that are part of the
% significant desync cluster at any time/frequency -- the full cluster,
% no further thresholding. (Replaces the union of top-quartile coverage
% and top-quartile peak |t|, which kept only the posterior core of the
% effect.) coverage/peakT are still computed and saved for reference.
%
% Time window = min/max extent (across all channels) of the significant
% cluster. Frequency SEARCH BAND = the fixed alpha band (alpha_band.m,
% 8-16 Hz, 2026-09-30) -- not the cluster's frequency extent.
%
%   pdat : struct with fields .alphaDir (Alpha_pre2post output dir,
%          containing stat_alpha.mat) and .out (where to save
%          PF_roi_config.mat)

L = load(fullfile(pdat.alphaDir, 'stat_alpha.mat'), 'stat_alpha');
stat = L.stat_alpha;

p_crit = stat.cfg.alpha;
neg_signif_clust = [];
if ~isempty(stat.negclusters)
    neg_pvals = [stat.negclusters(:).prob];
    neg_signif_clust = find(neg_pvals < p_crit);
end
if isempty(neg_signif_clust)
    error('PF_derive_roi:noCluster', 'No significant negative (desync) cluster found in stat_alpha.mat.');
end

mask = false(size(stat.negclusterslabelmat));
for sc = neg_signif_clust
    mask = mask | (stat.negclusterslabelmat == sc);
end

% dimord chan_freq_time (ft_freqstatistics default)
nChan = size(mask,1);
coverage = nan(nChan,1);
peakT    = nan(nChan,1);
for c = 1:nChan
    m = squeeze(mask(c,:,:));
    coverage(c) = mean(m(:));
    if any(m(:))
        t = squeeze(stat.stat(c,:,:));
        peakT(c) = max(abs(t(m)));
    else
        peakT(c) = 0;
    end
end

candidates = find(coverage > 0);
if isempty(candidates)
    error('PF_derive_roi:noChannels', 'No channels participate in the significant cluster.');
end

selChan = candidates; % full cluster

roiChannels = stat.label(selChan);

% Window + frequency band: min/max extent of the cluster across ALL
% channels (not just the selected ROI subset) -- "where the effect
% exists at all", same philosophy as the original design.
anyChanMask = squeeze(any(mask,1)); % freq x time
[freqIdx, timeIdx] = find(anyChanMask);
clusterFreq = [stat.freq(min(freqIdx)), stat.freq(max(freqIdx))]; % reported only
% PF SEARCH BAND = the fixed alpha band used by all alpha analyses
% (alpha_band.m, 8-16 Hz) -- NOT the cluster's frequency extent, which can
% reach into beta (e.g. 8-30 Hz) when alpha and beta desync form one
% cluster. Channels and time window still come from the cluster.
freqBand = alpha_band();
winLat   = [stat.time(min(timeIdx)), stat.time(max(timeIdx))];

fprintf('Derived PF ROI: %d channels, window=[%.2f %.2f]s, freq=[%.1f %.1f]Hz\n', ...
    length(roiChannels), winLat(1), winLat(2), freqBand(1), freqBand(2));
fprintf('(cluster frequency extent, for reference: %.1f-%.1f Hz)\n', clusterFreq(1), clusterFreq(2));
disp(roiChannels');

if ~exist(pdat.out,'dir'), mkdir(pdat.out); end
save(fullfile(pdat.out, 'PF_roi_config.mat'), 'roiChannels', 'winLat', 'freqBand', 'clusterFreq', 'coverage', 'peakT');
fprintf('Saved PF ROI config -> %s\n', fullfile(pdat.out, 'PF_roi_config.mat'));

end

function [mask, prof] = alpha_cluster_core(stat_alpha)
% The alpha_pre2post desync cluster and its per-channel profile, defined
% ONCE for every analysis that uses the alpha cluster as ROI
% (alpha_pre2post_plots.m: per-subject desync; ebosc_results.m: eBOSC
% measures).
%
%   mask : chan x freq x time, significant NEGATIVE (desync) cluster(s),
%          p < stat_alpha.cfg.alpha
%   prof : table from cluster_channel_profile.m (duration, strength, mass,
%          core); prof.core = the ROI channels (above-average cluster mass)

mask = false(size(stat_alpha.stat));
if isfield(stat_alpha,'negclusters') && ~isempty(stat_alpha.negclusters) % fields are absent when no cluster formed
    for sc = find([stat_alpha.negclusters.prob] < stat_alpha.cfg.alpha)
        mask = mask | (stat_alpha.negclusterslabelmat == sc);
    end
end
if ~any(mask(:))
    error('alpha_cluster_core:noCluster', 'No significant negative (desync) cluster in stat_alpha');
end
prof = cluster_channel_profile(mask, stat_alpha.stat, stat_alpha.label);
end

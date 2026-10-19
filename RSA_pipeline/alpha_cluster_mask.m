function mask = alpha_cluster_mask(stat_alpha)
% The alpha_pre2post desync cluster: union of all significant NEGATIVE
% clusters (p < stat_alpha.cfg.alpha), chan x freq x time. No channel
% ranking or selection (2026-10-02; replaces alpha_cluster_core.m).

mask = false(size(stat_alpha.stat));
if isfield(stat_alpha,'negclusters') && ~isempty(stat_alpha.negclusters) % fields are absent when no cluster formed
    for sc = find([stat_alpha.negclusters.prob] < stat_alpha.cfg.alpha)
        mask = mask | (stat_alpha.negclusterslabelmat == sc);
    end
end
if ~any(mask(:))
    error('alpha_cluster_mask:noCluster', 'No significant negative (desync) cluster in stat_alpha');
end
end

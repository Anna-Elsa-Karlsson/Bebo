function P = cluster_channel_profile(mask, stat, label)
% Per-channel profile of a significant cluster, and its CORE channels --
% the channels with the strongest / most sustained contribution. Used for
% BOTH the Ispec cluster (chan x time1 x time2) and the alpha_pre2post
% cluster (chan x freq x time), so the definition is identical.
%
%   mask  : chan x A x B logical, the significant cluster(s)
%   stat  : chan x A x B t-values (same size)
%   label : chan x 1 channel labels
%
% Returns table P (one row per channel):
%   inCluster  channel is part of the cluster at all
%   duration   share of the cluster's extent (A x B points where the
%              cluster exists on ANY channel) covered by this channel
%   strength   mean |t| over this channel's cluster points
%   mass       sum of |t| over this channel's cluster points
%              (= duration x extent size x strength: strength AND duration)
%   core       inCluster AND mass above the mean mass of the cluster
%              channels (= more than its 'fair share' of the cluster mass);
%              the ROI channels
% The ranking is DESCRIPTIVE: the cluster test shows that an effect exists
% in the cluster, not that it is strongest on these channels.

nCh  = size(mask,1);
M    = reshape(mask, nCh, []);
T    = reshape(abs(double(stat)), nCh, []);
ext  = any(M, 1);                                 % the cluster's extent (union over channels)
inCluster = any(M, 2);
nPts = sum(M, 2);
mass = sum(T .* M, 2);
strength = mass ./ max(nPts, 1); strength(~inCluster) = NaN;
duration = nPts / sum(ext);
core = inCluster & mass > mean(mass(inCluster));
P = table(label(:), inCluster, duration, strength, mass, core, ...
    'VariableNames', {'label','inCluster','duration','strength','mass','core'});
end

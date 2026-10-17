function stat = bebo_cluster_test(A, B, label, ax2, ax3, C, seed, sel)
% BEBO_CLUSTER_TEST  Group-level cluster permutation test, paired
% (depsamplesT: each person's A vs the same person's B), as in the plan:
% clusteralpha, minnbchan, maxsum, two-tailed, alpha .025, 1000
% permutations, neighbours by distance (55 mm) -- all from bebo_config.
%
%   A, B   [subject x chan x ax2 x ax3]  (B = zeros for a test against 0)
%   label  channel labels; ax2, ax3 = values of the 2nd and 3rd axis
%          (Ispec: time1 and time2, passed as FieldTrip 'freq' and 'time';
%           SME: frequency and time)
%   seed   random seed of the permutation draw (reproducible)
%   sel    optional struct with .frequency / .latency ranges to test
%
% Returns the FieldTrip stat structure (posclusters, negclusters, *labelmat).

if nargin < 8, sel = struct(); end
n = size(A, 1);
mk = @(X) struct('dimord', 'subj_chan_freq_time', 'label', {label(:)}, ...
    'freq', ax2(:)', 'time', ax3(:)', 'powspctrm', double(X));
dA = mk(A); dB = mk(B);

L = load(C.elecFile, 'elec');
elec = L.elec; elec.elecpos = elec.pnt; elec.chanpos = elec.pnt;
cfgn = []; cfgn.method = 'distance'; cfgn.neighbourdist = C.clus.neighbourdist;
neighbours = ft_prepare_neighbours(cfgn, elec);

cfg = [];
cfg.channel          = 'all';
cfg.parameter        = 'powspctrm';
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = C.clus.clusteralpha;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = C.clus.minnbchan;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = C.clus.alpha;
cfg.numrandomization = C.clus.numrandomization;
cfg.neighbours       = neighbours;
cfg.randomseed       = seed;
cfg.design           = [ones(1,n) 2*ones(1,n); 1:n 1:n];
cfg.ivar             = 1;
cfg.uvar             = 2;
if isfield(sel, 'frequency'), cfg.frequency = sel.frequency; end
if isfield(sel, 'latency'),   cfg.latency   = sel.latency;   end
ft_warning off 'FieldTrip:getdimord:warning_dimord_could_not_be_determined'
stat = ft_freqstatistics(cfg, dA, dB);
end

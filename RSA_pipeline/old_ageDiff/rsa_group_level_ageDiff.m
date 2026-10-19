function rsa_group_level_ageDiff(contrastType, freqLo, freqHi, pdat)
% Age difference in Ispec WITHOUT any ROI (2026-10-02): cluster-based
% permutation test YA vs OA on the per-subject 1st-level maps (the same
% Isim-vs-Csim t-maps that enter the fullSample test against 0), over the
% whole channel x time1 x time2 space. No selection is involved, so this
% is the unbiased test of the age difference (an ROI chosen on the pooled
% sample is not orthogonal to the group contrast with unequal group sizes,
% 37 YA vs 56 OA).
%
% Settings identical to rsa_group_level_fullSample.m except the design:
% indepsamplesT (YA = 1, OA = 2), montecarlo, clusteralpha 0.01, maxsum,
% minnbchan 2, two-tailed (alpha 0.025 per tail), 1000 permutations,
% neighbours 55 mm (Elec/elec.mat). Positive t = YA > OA.
%
% Requires rsa_group_level_perAge.m for both groups (1st-level stat files).
% Output: <variant>_ageDiff_cluster_stat.mat ('stat_age') in pdat.stat/<variant>/
%
%   pdat : .stat, .elec

switch contrastType
    case 'Ispec', nameA = 'within-item'; nameB = 'within-cat';
    otherwise, error('rsa_group_level_ageDiff:contrastType', 'Only Ispec is supported');
end
variantLabel = sprintf('%s_%dto%dHz', contrastType, freqLo, freqHi);
statDir = fullfile(pdat.stat, variantLabel);

Lya = load(fullfile(statDir, ['YA_' nameA '_vs_' nameB '_1st_level_stat_vsBL.mat']), 'stat');
Loa = load(fullfile(statDir, ['OA_' nameA '_vs_' nameB '_1st_level_stat_vsBL.mat']), 'stat');
ya = Lya.stat; oa = Loa.stat; clear Lya Loa
nYA = size(ya.stat,1); nOA = size(oa.stat,1);
ya.stat = ft_struct2single(ya.stat);
oa.stat = ft_struct2single(oa.stat);

cfg = load(fullfile(pdat.elec, 'elec.mat')); % as rsa_group_level_fullSample.m
cfg.elec.elecpos = cfg.elec.pnt;
cfg.elec.chanpos = cfg.elec.pnt;
cfg.layout = ft_prepare_layout(cfg);

cfg.channel = 'all'; cfg.latency = 'all'; cfg.frequency = 'all';
cfg.avgoverchan = 'no'; cfg.avgovertime = 'no'; cfg.avgoverfreq = 'no';
cfg.method = 'montecarlo'; cfg.statistic = 'indepsamplesT';
cfg.correctm = 'cluster'; cfg.clusteralpha = 0.01; cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2; cfg.tail = 0; cfg.clustertail = 0; cfg.alpha = 0.025;
cfg.numrandomization = 1000;
cfg.parameter = 'stat';

cfg_neighb = []; cfg_neighb.method = 'distance'; cfg_neighb.neighbourdist = 55;
ft_warning off 'FieldTrip:getdimord:warning_dimord_could_not_be_determined'
cfg.neighbours = ft_prepare_neighbours(cfg_neighb, cfg.elec);

cfg.design = [ones(1,nYA), 2*ones(1,nOA)];
cfg.ivar = 1;

stat_age = ft_freqstatistics(cfg, ya, oa); %#ok<NASGU>
save(fullfile(statDir, [variantLabel '_ageDiff_cluster_stat.mat']), 'stat_age', '-v7.3')
fprintf('Saved age-difference cluster stat (%d YA vs %d OA) for %s -> %s\n', nYA, nOA, variantLabel, statDir);
end

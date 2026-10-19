function rsa_group_level_fullSample(contrastType, freqLo, freqHi, pdat, simField)
% Consolidated fullSample step: 2nd-level cluster stat (step4b, YA+OA
% pooled vs a zero comparator) + 3rd-level cluster extraction (step4c,
% per-subject cluster-mean values + YA vs OA t-test).
%
% Requires rsa_group_level_perAge.m to have already completed for BOTH
% 'YA' and 'OA' for this exact contrastType/freqLo/freqHi -- run this
% only after both of those jobs finish (explicit grid dependency in
% submit_rsa_all_variants.sh, not an implicit assumption about run order
% like the original pipeline had).
%
%   contrastType  : 'Ispec' or 'Catspec'
%   freqLo,freqHi : analysis sub-band in Hz
%   pdat          : struct with fields .ga, .stat, .elec (paths; .elec = dir with elec.mat)
%   simField      : (optional, default 'powspctrm') 'powspctrm'
%                   (baseline-corrected, primary) or 'powspctrmOrig'
%                   (uncorrected) -- 2026-09-XX addition, must match
%                   whatever rsa_group_level_perAge.m was run with for
%                   these subjects. Loads the corresponding GA files
%                   (sim_GA.mat vs sim_GA_orig.mat, both already saved
%                   regardless by that step) and 1st-level stat file
%                   (plain vs '_uncorrected' suffix), and tags all of
%                   this step's own outputs with the same suffix so the
%                   two parallel runs never overwrite each other.
%
% Output naming: uses a generic "specA_minus_specB" framing so it works
% for both contrast types without contrast-specific hardcoding --
% Ispec:   spec = within-item minus within-cat  (matches original "ispec")
% Catspec: spec = within-cat  minus between-cat (matches original "catspec")

if nargin < 5 || isempty(simField), simField = 'powspctrm'; end
gaSuffix = ''; statSuffix = '';
if ~strcmp(simField, 'powspctrm')
    gaSuffix = '_orig';      % load sim_GA_orig.mat instead of sim_GA.mat
    statSuffix = '_uncorrected';
end
% Note: regardless of gaSuffix, the loaded GA struct's own field is always
% named 'powspctrm' (that's just what grand_average_one_condition calls
% it internally in both files) -- the suffix only picks which FILE to load.

switch contrastType
    case {'Ispec', 'IspecAll'} % IspecAll (2026-10-18): Ispec without the response filter, own folders
        nameA = 'within-item'; nameB = 'within-cat';
    case 'IspecPaired' % 2026-10-19: item-specific Csim, row-aligned with Isim -> paired 1st level
        nameA = 'within-item'; nameB = 'within-catItem';
    case 'Catspec'
        nameA = 'within-cat';  nameB = 'between-cat';
    otherwise
        error('rsa_group_level_fullSample:contrastType', 'Unknown contrastType "%s"', contrastType)
end

bandLabel = sprintf('%dto%dHz', freqLo, freqHi);
variantLabel = sprintf('%s_%s', contrastType, bandLabel);
variantOutLabel = [variantLabel statSuffix];
gaDir = fullfile(pdat.ga, variantLabel);
statDir = fullfile(pdat.stat, variantLabel);

%% === Step 4b equivalent: 2nd-level cluster stat, YA+OA pooled vs zero ===

statIya = load(fullfile(statDir, ['YA_' nameA '_vs_' nameB '_1st_level_stat_vsBL' statSuffix '.mat']), 'stat');
statIoa = load(fullfile(statDir, ['OA_' nameA '_vs_' nameB '_1st_level_stat_vsBL' statSuffix '.mat']), 'stat');

tmp = [];
for y = 1:size(statIya.stat.stat,1), tmp(end+1,:,:,:) = statIya.stat.stat(y,:,:,:); end %#ok<AGROW>
for y = 1:size(statIoa.stat.stat,1), tmp(end+1,:,:,:) = statIoa.stat.stat(y,:,:,:); end %#ok<AGROW>
stat_I = statIya.stat;
stat_I.stat = tmp;
stat_I.df = [statIya.stat.df, statIoa.stat.df];
clear statIya statIoa

cfg = load(fullfile(pdat.elec, 'elec.mat')); % Elec/elec.mat -- same file as alpha_pre2post_group_level.m (2026-09-26)
cfg.elec.elecpos = cfg.elec.pnt;
cfg.elec.chanpos = cfg.elec.pnt;
cfg.layout = ft_prepare_layout(cfg);

stat = stat_I;
null = stat_I; null.stat = zeros(size(stat.stat));

cfg.channel = 'all'; cfg.latency = 'all'; cfg.frequency = 'all';
cfg.avgoverchan = 'no'; cfg.avgovertime = 'no'; cfg.avgoverfreq = 'no';
% ONE-SAMPLE TEST AGAINST 0 via depsamplesT (2026-09-26 fix): each
% subject is paired with their own zero (cfg.uvar), which gives the
% one-sample t (df = n-1) and a sign-flip permutation null. The previous
% indepsamplesT treated the zeros as a separate group of n subjects --
% same t-value but df = 2n-2, and permutations mixed real values and
% zeros across the two "groups", which is not a valid null. Same fix as
% alpha_pre2post_group_level.m.
cfg.method = 'montecarlo'; cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster'; cfg.clusteralpha = 0.01; cfg.clusterstatistic = 'maxsum'; % 2026-10-01: 0.01 (was 0.05), as the alpha_pre2post test
cfg.minnbchan = 1; cfg.tail = 0; cfg.clustertail = 0; cfg.alpha = 0.025; % minnbchan 1 = setting of the final Ispec cluster stat (2026-10-03)
cfg.numrandomization = 1000; % 2026-09-26: 1000 (MC SE of p near 0.025 ~0.005; with 500 it was ~0.007)
cfg.parameter = 'stat';

cfg_neighb = []; cfg_neighb.method = 'distance';
% neighbourdist set EXPLICITLY (2026-09-26): the original step4b scripts
% left it unset, i.e. FieldTrip's default 40mm (~3.4 neighbours/channel,
% min 1). 55mm (~6 neighbours, min 2) = same as alpha_pre2post_group_level.m.
cfg_neighb.neighbourdist = 55;
ft_warning off 'FieldTrip:getdimord:warning_dimord_could_not_be_determined'
cfg.neighbours = ft_prepare_neighbours(cfg_neighb, cfg.elec);

subj = size(stat.stat,1);
design = zeros(2,subj*2);
design(1,:) = [ones(1,subj), 2*ones(1,subj)]; % condition: real vs zero
design(2,:) = [1:subj, 1:subj];               % subject: pairs each subject with their own zero
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design = design;

stat.stat = ft_struct2single(stat.stat);
null.stat = ft_struct2single(null.stat);

stat2 = ft_freqstatistics(cfg, stat, null); %#ok<NASGU>
save(fullfile(statDir, [variantOutLabel '_2nd_level_cluster_stat_vsBL.mat']), 'stat2', '-v7.3')
fprintf('Saved 2nd-level cluster stat (simField=%s) for %s -> %s\n', simField, variantOutLabel, statDir);

%% === Step 4c equivalent: extract cluster means, per-subject spec = A-B, YA vs OA ===

spec = struct(); specA = struct(); specB = struct();

for groups = {'YA','OA'}
    group = cell2mat(groups);

    L = load(fullfile(gaDir, [group '_' nameA '_sim_GA' gaSuffix '.mat']), 'sim'); wi = L.sim;
    L = load(fullfile(gaDir, [group '_' nameB '_sim_GA' gaSuffix '.mat']), 'sim'); wc = L.sim;

    p_crit = stat2.cfg.alpha;

    pos_signif_clust = [];
    if isfield(stat2,'posclusters') && ~isempty(stat2.posclusters) % fields are absent when no cluster formed
        pos_cluster_pvals = [stat2.posclusters(:).prob];
        pos_signif_clust = find(pos_cluster_pvals < p_crit);
        disp([num2str(length(pos_signif_clust)) ' positive cluster(s) with p < ' num2str(p_crit)])
    end
    % POSITIVE CLUSTERS ONLY (2026-09-26): same convention as
    % rsa_singleTrial_merged.m and Ispec_derive_alphaROI.m, so the
    % cluster-mean values, the single-trial similarity mask and the alpha
    % ROI all come from the same cluster(s). The effect of interest is
    % A > B (item > cat); a negative cluster would be a different effect
    % and averaging it in could cancel part of the positive one.
    if isfield(stat2,'negclusters') && ~isempty(stat2.negclusters) && any([stat2.negclusters(:).prob] < p_crit)
        warning('%d significant NEGATIVE cluster(s) found -- NOT included in the cluster means', sum([stat2.negclusters(:).prob] < p_crit));
    end
    %
    % BUGFIX (2026-09-05): must iterate over the actual cluster index
    % VALUES in pos_signif_clust/neg_signif_clust, not over 1:count --
    % the original `for sc = 1:length(pos_signif_clust)` silently only
    % worked when the significant clusters happened to be exactly
    % {1,...,N} (FieldTrip's top-ranked clusters); with any other
    % significant-cluster pattern (e.g. clusters {2,3} significant but
    % not 1) it would grab the wrong cluster(s) and/or miss a real one.
    mask = NaN(size(stat2.stat));
    for sc = pos_signif_clust, mask(stat2.posclusterslabelmat == sc) = 1; end

    for s = 1:size(wi.powspctrm,1)
        wi.powspctrm2(s,:,:,:) = single(squeeze(wi.powspctrm(s,:,:,:)) .* mask);
        wc.powspctrm2(s,:,:,:) = single(squeeze(wc.powspctrm(s,:,:,:)) .* mask);
    end

    cfg2 = []; cfg2.parameter = 'powspctrm2'; cfg2.operation = 'subtract';
    ft_warning off 'FieldTrip:getdimord:warning_dimord_could_not_be_determined'
    specStruct = ft_math(cfg2, wi, wc);

    cfg3 = []; cfg3.avgoverchan = 'yes'; cfg3.avgovertime = 'yes'; cfg3.avgoverfreq = 'yes'; cfg3.nanmean = 'yes';
    wi_means = ft_selectdata(cfg3, wi); specA.(group) = wi_means.powspctrm2;
    wc_means = ft_selectdata(cfg3, wc); specB.(group) = wc_means.powspctrm2;
    spec_means = ft_selectdata(cfg3, specStruct); spec.(group) = spec_means.powspctrm2;

    clear wi wc wi_means wc_means specStruct spec_means
end

[H,P,CI,STATS] = ttest2(spec.YA, spec.OA); %#ok<ASGLU>
YA_vs_OA.H = H; YA_vs_OA.P = P; YA_vs_OA.CI = CI; YA_vs_OA.STATS = STATS;

n1 = length(spec.YA); n2 = length(spec.OA);
pooled_sd = sqrt(((n1-1)*var(spec.YA) + (n2-1)*var(spec.OA)) / (n1+n2-2));
YA_vs_OA.cohensD = (mean(spec.YA) - mean(spec.OA)) / pooled_sd;
YA_vs_OA.M_YA = mean(spec.YA); YA_vs_OA.SD_YA = std(spec.YA);
YA_vs_OA.M_OA = mean(spec.OA); YA_vs_OA.SD_OA = std(spec.OA);

save(fullfile(statDir, [variantOutLabel '_spec_in_cluster_3rd_level_stats_Pairwise.mat']), 'YA_vs_OA', 'STATS')
save(fullfile(gaDir, [variantOutLabel '_mean_specA_in_cluster.mat']), 'specA') % A = within-item (Ispec) / within-cat (Catspec)
save(fullfile(gaDir, [variantOutLabel '_mean_specB_in_cluster.mat']), 'specB') % B = within-cat (Ispec) / between-cat (Catspec)
save(fullfile(gaDir, [variantOutLabel '_mean_spec_in_cluster.mat']), 'spec')   % spec = A - B

fprintf('%s (simField=%s): YA M=%.4f SD=%.4f, OA M=%.4f SD=%.4f, t(%.0f)=%.3f, p=%.4f, d=%.3f\n', ...
    variantOutLabel, simField, YA_vs_OA.M_YA, YA_vs_OA.SD_YA, YA_vs_OA.M_OA, YA_vs_OA.SD_OA, ...
    STATS.df, STATS.tstat, YA_vs_OA.P, YA_vs_OA.cohensD);

end

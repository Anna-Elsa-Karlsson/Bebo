function alpha_pre2post_group_level(pdat)
% Group-level step for the alpha activation-vs-zero (pre-to-post)
% analysis: loads all subjects' per-subject descriptives (from
% alpha_pre2post_subject_level.m), builds grand averages (all/YA/OA),
% and runs the depsamplesT-vs-zero cluster stat (paired design -- each
% subject's real value against their own constant zero, NOT an
% indepsamplesT against a fake duplicated-N zero "group", which inflates
% df and was the root cause of the earlier "everything is significant"
% result).
%
%   pdat : struct with fields .out (per-subject files + this step's
%          output), .elec (dir containing elec.mat)
%
% Output: saves GApow.mat, GApow_YA.mat, GApow_OA.mat, stat_alpha.mat
% into pdat.out.

ids = ids_ALL('all');

AllSubjPow = {};
AllSubjPow_YA = {};
AllSubjPow_OA = {};

for j = 1:length(ids)
    subj = num2str(ids(j));
    L = load(fullfile(pdat.out, [subj '_alphaPow.mat']), 'pow');
    AllSubjPow{end+1,1} = L.pow; %#ok<AGROW>
    if ids(j) < 4000
        AllSubjPow_YA{end+1,1} = L.pow; %#ok<AGROW>
    else
        AllSubjPow_OA{end+1,1} = L.pow; %#ok<AGROW>
    end
end

%% Grand averages
cfgPow = [];
cfgPow.keepindividual = 'yes';
cfgPow.parameter = 'powspctrm';

GApow    = ft_freqgrandaverage(cfgPow, AllSubjPow{:});
GApow_YA = ft_freqgrandaverage(cfgPow, AllSubjPow_YA{:});
GApow_OA = ft_freqgrandaverage(cfgPow, AllSubjPow_OA{:});

%% Neighbours -- exact settings from the validated local script
load(fullfile(pdat.elec, 'elec.mat'), 'elec');
cfg_neigh = [];
cfg_neigh.elec = elec;
cfg_neigh.method = 'distance';
cfg_neigh.neighbourdist = 55;
neighbours = ft_prepare_neighbours(cfg_neigh);

%% Test against 0 -- paired (depsamplesT), each subject vs their own zero
stat_real = GApow;
stat_null = GApow;
stat_null.powspctrm = zeros(size(GApow.powspctrm));

cfg = [];
cfg.channel          = 'all';
cfg.latency          = alpha_time_window();   % 0.4-1.7 s (alpha_time_window.m): after the evoked response
cfg.frequency        = [6 30];
cfg.avgoverchan      = 'no'; cfg.avgovertime = 'no'; cfg.avgoverfreq = 'no';
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.01;                 % final desync cluster stat (2026-10-02)
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 4;                 % 4 (2026-10-02: briefly 2, back to 4 -- user)
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000; % 2026-09-26: 1000 (was 500)
cfg.parameter        = 'powspctrm';
cfg.neighbours       = neighbours;

nSubj = size(stat_real.powspctrm, 1);
cfg.design      = zeros(2, nSubj*2);
cfg.design(1,:) = [1:nSubj, 1:nSubj];
cfg.design(2,:) = [ones(1,nSubj), 2*ones(1,nSubj)];
cfg.uvar = 1;
cfg.ivar = 2;

stat_alpha = ft_freqstatistics(cfg, stat_real, stat_null); %#ok<NASGU>

%% Save
if ~exist(pdat.out, 'dir'), mkdir(pdat.out); end
save(fullfile(pdat.out, 'stat_alpha.mat'), 'stat_alpha', '-v7.3')
save(fullfile(pdat.out, 'GApow.mat'), 'GApow', '-v7.3')
save(fullfile(pdat.out, 'GApow_YA.mat'), 'GApow_YA', '-v7.3')
save(fullfile(pdat.out, 'GApow_OA.mat'), 'GApow_OA', '-v7.3')
fprintf('Alpha pre2post group-level done -> %s\n', pdat.out);

end

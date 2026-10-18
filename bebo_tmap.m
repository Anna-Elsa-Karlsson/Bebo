function t = bebo_tmap(X1, X2, label, ax2, ax3, type)
% BEBO_TMAP  First-level (within-person) t-map: FieldTrip ft_freqstatistics,
% analytic (as the original SME, TFR_GA_old.m, and the Ispec t-map). The
% t-maps are then tested against 0 at the group level.
%
%   X1, X2  [trial x chan x ax2 x ax3]
%   type    'indep'  indepsamplesT: two different sets of trials
%                    (later 'old' vs later 'similar')
%           'dep'    depsamplesT: the same trials in two states (row k of X1
%                    and X2 = trial k; post-stimulus vs the trial's baseline)
% Returns t [chan x ax2 x ax3].

mk = @(X) struct('dimord', 'rpt_chan_freq_time', 'label', {label(:)}, ...
    'freq', ax2(:)', 'time', ax3(:)', 'powspctrm', double(X));
n1 = size(X1, 1); n2 = size(X2, 1);
cfg = [];
cfg.channel = 'all'; cfg.latency = 'all'; cfg.frequency = 'all';
cfg.avgoverchan = 'no'; cfg.avgovertime = 'no'; cfg.avgoverfreq = 'no';
cfg.parameter = 'powspctrm';
cfg.method    = 'analytic';
cfg.tail      = 0;
switch type
    case 'indep'
        cfg.statistic = 'indepsamplesT';
        cfg.design    = [ones(1, n1) 2*ones(1, n2)];
        cfg.ivar      = 1;
    case 'dep'
        if n1 ~= n2, error('bebo_tmap:dep', 'depsamplesT needs the same trials in X1 and X2'); end
        cfg.statistic = 'depsamplesT';
        cfg.design    = [ones(1, n1) 2*ones(1, n1); 1:n1 1:n1];
        cfg.ivar      = 1;
        cfg.uvar      = 2;
    otherwise
        error('bebo_tmap:type', 'type must be ''indep'' or ''dep''');
end
ft_warning off 'FieldTrip:getdimord:warning_dimord_could_not_be_determined'
st = ft_freqstatistics(cfg, mk(X1), mk(X2));
t = st.stat;
end

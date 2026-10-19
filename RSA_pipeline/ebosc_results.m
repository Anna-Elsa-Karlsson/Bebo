function ebosc_results(pdat)
% eBOSC results (rewritten 2026-10-05): per-subject rhythm measures PER
% BAND (alpha, beta; ebosc_bands.m) in the ISPEC CLUSTER channels, computed
% from the episodes saved by ebosc_subject.m, + all statistics.
% Computation only -- ebosc_plots.m only loads what is saved here.
%
% Channels: all channels of the Ispec cluster (Ispec_alphaROI_config.mat),
% i.e. where the main analyses measure alpha. Time: the eBOSC analysis
% window 0-1.5 s.
%
% Per subject and band (an episode belongs to a band if lo <= its mean
% frequency < hi), first per channel, then the plain mean over the cluster
% channels:
%   abundance  fraction of the analysis window covered by episodes of the
%              band (mean over trials; 0 if a channel has none)
%   freq       MAIN frequency measure: duration-weighted mean episode
%              frequency (channels without episodes skipped)
%   amp        amplitude above background: mean log10(power/background)
%              over all episode time points (channels without episodes skipped)
%   PF         secondary: highest local maximum of the band's episode
%              spectrum within the band (NaN if none); PF_nChan = channels
%              contributing
%   nChanEp    number of cluster channels with at least one episode
% Also per subject: spectrum above background (log10(power/background),
% all time points, mean over the cluster channels) -> individual spectra;
% and PF per band on EVERY channel -> group-mean PF topographies (figure).
% Stats per band: YA vs OA (Welch t, Cohen d; Holm over freq, abundance,
% amp within the band), PF and nChanEp YA vs OA (reported separately);
% freq x Ispec (partial r controlling age, within-group r, Fisher z).
%
%   pdat : .subj (per-subject *_ebosc.mat), .stm (Ispec_alphaROI_config.mat),
%          .ga (Results/GA, per-subject Ispec), .out

ID_ya = ids_ALL('ya'); ID_oa = ids_ALL('oa');
ids = [ID_ya(:); ID_oa(:)];
nYA = numel(ID_ya); nOA = numel(ID_oa); nAll = nYA + nOA;
idx_YA = 1:nYA; idx_OA = nYA+1:nAll;
age_vec = [ones(nYA,1); 2*ones(nOA,1)];

L = load(fullfile(pdat.stm, 'Ispec_alphaROI_config.mat'), 'roiChannels');
clusterChans = L.roiChannels(:);

M = struct(); specAbove = []; F = []; bands = []; topoLabel = {}; pfAll = struct();
for s = 1:nAll
    E = load(fullfile(pdat.subj, sprintf('%d_ebosc.mat', ids(s))));
    if s == 1
        F = E.F; bands = E.bands; nB = numel(bands.names);
        specAbove = nan(nAll, numel(F));
        topoLabel = E.label(:);
        for b = 1:nB, pfAll.(bands.names{b}) = nan(nAll, numel(topoLabel)); end
        for b = 1:nB
            for f = {'abundance','freq','amp','PF','PF_nChan','nChanEp'}
                M.(bands.names{b}).(f{1}) = nan(nAll,1);
            end
            M.(bands.names{b}).specEp = nan(nAll, numel(F));
        end
    end
    chIdx = find(ismember(E.label, clusterChans));
    specAbove(s,:) = mean(E.specAll_log10(chIdx,:) - E.bg_log10(chIdx,:), 1, 'omitnan');
    ep = E.ep;
    for b = 1:nB
        rng_ = bands.ranges(b,:);
        inB = ep.freq >= rng_(1) & ep.freq < rng_(2);
        fi = find(F >= rng_(1) & F < rng_(2));
        ab = zeros(numel(chIdx),1); fr = nan(numel(chIdx),1); am = nan(numel(chIdx),1); pf = nan(numel(chIdx),1);
        for k = 1:numel(chIdx)
            r = find(inB & ep.chan == chIdx(k));
            if isempty(r), continue; end
            cover = false(E.nTrials, E.nWin);
            for q = r'
                cover(ep.trial(q), ep.on(q):ep.off(q)) = true;
            end
            ab(k) = mean(mean(cover, 2));
            fr(k) = sum(ep.freq(r) .* ep.nPts(r)) / sum(ep.nPts(r));
            am(k) = sum(ep.snrSum(r)) / sum(ep.nPts(r));
            sp = E.specEp_log10(chIdx(k), fi, b);
            if ~all(isnan(sp))
                mx = THG_find_local_max_1D_20150117(sp, 1);
                if ~isempty(mx), [~, w] = max(sp(mx)); pf(k) = F(fi(mx(w))); end
            end
        end
        nm = bands.names{b};
        M.(nm).abundance(s) = mean(ab);
        M.(nm).freq(s)      = mean(fr, 'omitnan');
        M.(nm).amp(s)       = mean(am, 'omitnan');
        M.(nm).PF(s)        = mean(pf, 'omitnan');
        M.(nm).PF_nChan(s)  = sum(~isnan(pf));
        M.(nm).nChanEp(s)   = sum(~isnan(fr));
        M.(nm).specEp(s,:)  = mean(squeeze(E.specEp_log10(chIdx,:,b)), 1, 'omitnan');
        % PF on every channel (for the topographies)
        [okT, iT] = ismember(topoLabel, E.label);
        for k = find(okT)'
            sp = E.specEp_log10(iT(k), fi, b);
            if all(isnan(sp)), continue; end
            mx = THG_find_local_max_1D_20150117(sp, 1);
            if ~isempty(mx), [~, w] = max(sp(mx)); pfAll.(nm)(s,k) = F(fi(mx(w))); end
        end
    end
end

%% Per-subject Ispec (cluster means, ids_ALL order)
Ls = load(fullfile(pdat.ga, 'Ispec_30to100Hz', 'Ispec_30to100Hz_mean_spec_in_ROI.mat'), 'spec');
ispec = double([Ls.spec.YA(:); Ls.spec.OA(:)]);

%% Stats per band
grp = struct(); corrF = struct();
for b = 1:nB
    nm = bands.names{b};
    meas = {'freq','abundance','amp','PF','nChanEp'};
    p_w = nan(1, numel(meas));
    for m = 1:numel(meas)
        v = M.(nm).(meas{m});
        a = v(idx_YA); c = v(idx_OA); a = a(~isnan(a)); c = c(~isnan(c));
        [~, p_w(m), ~, stw] = ttest2(a, c, 'Vartype','unequal');
        grp.(nm).(meas{m}) = struct('M_YA',mean(a),'SD_YA',std(a),'n_YA',numel(a),'M_OA',mean(c),'SD_OA',std(c),'n_OA',numel(c), ...
            't',stw.tstat,'df',stw.df,'p',p_w(m), ...
            'd',(mean(a)-mean(c))/sqrt(((numel(a)-1)*var(a)+(numel(c)-1)*var(c))/(numel(a)+numel(c)-2)));
    end
    ph = holm(p_w(1:3));
    for m = 1:3, grp.(nm).(meas{m}).p_Holm = ph(m); end
    corrF.(nm) = corr_with_ispec(M.(nm).freq, ispec, age_vec);
end

% group spectra (mean +/- SE): above background (all time points) and per-band episode spectra
spec = struct('F', F, 'bands', bands);
% group-mean PF per channel (subjects without a peak on that channel skipped)
topo = struct('label', {topoLabel});
for b = 1:nB
    nm = bands.names{b};
    topo.(nm).PF_YA = mean(pfAll.(nm)(idx_YA,:), 1, 'omitnan')';
    topo.(nm).PF_OA = mean(pfAll.(nm)(idx_OA,:), 1, 'omitnan')';
    topo.(nm).n_YA  = sum(~isnan(pfAll.(nm)(idx_YA,:)), 1)';
    topo.(nm).n_OA  = sum(~isnan(pfAll.(nm)(idx_OA,:)), 1)';
end
spec.above_YA = specAbove(idx_YA,:); spec.above_OA = specAbove(idx_OA,:);
for b = 1:nB
    nm = bands.names{b};
    spec.(nm).ep_YA_M = mean(M.(nm).specEp(idx_YA,:),1,'omitnan'); spec.(nm).ep_OA_M = mean(M.(nm).specEp(idx_OA,:),1,'omitnan');
end

T = table(ids, double(age_vec == 1), ispec, 'VariableNames', {'ID','YA','Ispec'});
for b = 1:nB
    nm = bands.names{b};
    for f = {'freq','abundance','amp','PF','PF_nChan','nChanEp'}
        T.([nm '_' f{1}]) = M.(nm).(f{1});
    end
end

%% Save
if ~exist(pdat.out,'dir'), mkdir(pdat.out); end
save(fullfile(pdat.out, 'eBOSC_results.mat'), 'T', 'grp', 'corrF', 'spec', 'topo', 'clusterChans', 'idx_YA', 'idx_OA', 'age_vec', 'bands');
writetable(T, fullfile(pdat.out, 'eBOSC_per_subject.csv'));

fid = fopen(fullfile(pdat.out, 'eBOSC_stats.txt'), 'w');
fprintf(fid, 'eBOSC rhythms (Kosciessa et al., 2020) -- %s\n', datestr(now));
fprintf(fid, 'Bands (ebosc_bands.m): %s; analysis window 0-1.5 s; channels: %d Ispec cluster channels.\n', ...
    strjoin(arrayfun(@(b) sprintf('%s %g-%g Hz', bands.names{b}, bands.ranges(b,1), bands.ranges(b,2)), 1:nB, 'uni', 0), ', '), numel(clusterChans));
fprintf(fid, 'N = %d (%d YA, %d OA). freq = duration-weighted mean episode frequency (main); PF = local maximum of the episode spectrum (secondary).\n\n', nAll, nYA, nOA);
lab = struct('freq','Frequency (Hz; mean episode frequency)','abundance','Abundance (fraction of time with the rhythm)', ...
    'amp','Amplitude above background (log10 ratio)','PF','Peak frequency of the episode spectrum (Hz; secondary)', ...
    'nChanEp','Cluster channels with at least one episode');
for b = 1:nB
    nm = bands.names{b};
    fprintf(fid, '==================== %s (%g-%g Hz) ====================\n', upper(nm), bands.ranges(b,1), bands.ranges(b,2));
    fprintf(fid, '--- YA vs OA (Welch t; Holm over frequency, abundance, amplitude) ---\n');
    for m = {'freq','abundance','amp','PF','nChanEp'}
        g = grp.(nm).(m{1});
        fprintf(fid, '%s\n   YA: M=%.3f SD=%.3f (n=%d)  OA: M=%.3f SD=%.3f (n=%d)\n   Welch t(%.1f)=%.3f, p=%.4f', ...
            lab.(m{1}), g.M_YA, g.SD_YA, g.n_YA, g.M_OA, g.SD_OA, g.n_OA, g.df, g.t, g.p);
        if isfield(g, 'p_Holm'), fprintf(fid, ', p_Holm=%.4f', g.p_Holm); end
        fprintf(fid, ', d=%.3f\n\n', g.d);
    end
    c = corrF.(nm);
    fprintf(fid, '--- %s frequency x Ispec ---\npartial r (age)=%.3f, p=%.4f | YA r=%.3f (n=%d), OA r=%.3f (n=%d) | Fisher z=%.2f, p=%.4f\n\n', ...
        nm, c.r_part, c.p_part, c.r_YA, c.n_YA, c.r_OA, c.n_OA, c.z, c.p_z);
end
fclose(fid);
fprintf('eBOSC results saved -> %s\n', pdat.out);
end


function C = corr_with_ispec(x, ispec, age_vec)
% partial r controlling age, within-group r, Fisher z YA vs OA
ok = ~isnan(x) & ~isnan(ispec);
[r_part, p_part] = partialcorr(x(ok), ispec(ok), age_vec(ok));
oY = ok & age_vec == 1; oO = ok & age_vec == 2;
r_YA = corr(x(oY), ispec(oY)); r_OA = corr(x(oO), ispec(oO));
z_f = (atanh(r_YA) - atanh(r_OA)) / sqrt(1/(sum(oY)-3) + 1/(sum(oO)-3));
p_f = 2 * (1 - normcdf(abs(z_f)));
C = struct('r_part',r_part,'p_part',p_part,'r_YA',r_YA,'n_YA',sum(oY),'r_OA',r_OA,'n_OA',sum(oO),'z',z_f,'p_z',p_f);
end

function p_adj = holm(p)
m = numel(p); [ps, order] = sort(p(:)');
adj = min(1, cummax((m - (1:m) + 1) .* ps));
p_adj = nan(size(p)); p_adj(order) = adj;
end

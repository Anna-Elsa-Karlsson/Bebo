function ebosc_subject(subj, pdat)
% eBOSC rhythm detection for one subject (Kosciessa et al., 2020), adapted
% from bosc_main.m (CFC_CMY_bosc project). Replaces the FFT-based PF
% estimation.
%
% Per channel (all EEG channels; channel selection happens in
% ebosc_results.m): wavelet transform per trial -> robust 1/f background
% fit over all trials (theta / alpha / beta peaks left out, ebosc_bands.m)
% -> power + duration thresholds -> rhythmic episodes.
%
% 2026-10-05: EVERY detected episode is saved (not only alpha summaries),
% so ebosc_results.m can apply any band without rerunning this step.
%
% Trials: presentations 1-2, exemplars 1-2, 'old'/'similar' (same rule
% as every other analysis). Segment -2..3.5 s; after trialPad (2 s each
% side) the analysis window is 0-1.5 s (encoding).
%
% Saves (pdat.out/<ID>_ebosc.mat):
%   label, F, tWin, nWin (samples in the analysis window), nTrials, bands
%   bg_log10, pt_log10   chan x freq: background fit and power threshold
%   specAll_log10        chan x freq: mean log10 power over all time points
%   specEp_log10         chan x freq x band: mean log10 power over the time
%                        points inside episodes of that band (NaN if none)
%   ep                   one row per detected episode (all frequencies):
%                        chan, trial, freq (mean episode frequency), on, off
%                        (samples in the analysis window), nPts, snrSum (sum
%                        over the episode's points of log10(power/background)
%                        at the episode's frequency)
%   pdat : .rawEnc (PreprocData root), .fun (eBOSC function folder), .out

subj = num2str(subj);
addpath(pdat.fun);

%% eBOSC parameters (as bosc_main.m unless noted)
cfg = [];
cfg.eBOSC.F                 = 2.^(0:.125:6);                    % 1-64 Hz, log-spaced
cfg.eBOSC.wavenumber        = 6;
cfg.eBOSC.ncyc              = repmat(3, 1, numel(cfg.eBOSC.F)); % duration threshold (cycles)
cfg.eBOSC.percentile        = .95;                              % power threshold (percentile of background)
cfg.eBOSC.fsample           = 250;
cfg.eBOSC.WLpadding         = 250;                              % wavelet edge padding [samples]
cfg.eBOSC.detectedPad       = 250;                              % shoulder for duration threshold
cfg.eBOSC.trialPad          = 500;                              % WLpadding + detectedPad
cfg.eBOSC.BGpad             = 250;                              % CHANGED (template: 750): BEBO segments are 5.5 s, 750 would leave no data; 250 = WLpadding removes wavelet edge artefacts
cfg.eBOSC.fstp              = 1;
cfg.eBOSC.BiasCorrection    = 'yes';                            % CHANGED (template: 'no'): temporal wavelet-bias correction (Kosciessa et al., 2020); without it abundance is inflated by the wavelet's temporal spread
cfg.eBOSC.method            = 'MaxBias';
cfg.eBOSC.edgeOnly          = 'yes';
cfg.eBOSC.effSignal         = 'PT';
cfg.eBOSC.waveseg           = [-2 3.5];                         % full BEBO epoch; analysis window after trialPad = 0-1.5 s
bands = ebosc_bands();
cfg.eBOSC.exclFranges       = bands.exclBG;                     % peaks removed before the background fit (theta, alpha, beta)

%% Load + select trials + channels + segment
L = load(fullfile(pdat.rawEnc, [subj '_Enc_allcond.mat']), 'data');
dat = L.data; clear L
if dat.fsample ~= cfg.eBOSC.fsample
    rc = []; rc.resamplefs = cfg.eBOSC.fsample; rc.detrend = 'no'; rc.demean = 'no';
    ti2 = dat.trialinfo2; dat = ft_resampledata(rc, dat); dat.trialinfo2 = ti2;
end
ti = dat.trialinfo2;
trl = find(cell2mat(ti(:,9)) <= 2 & cell2mat(ti(:,8)) < 3 & (strcmp(ti(:,11),'old') | strcmp(ti(:,11),'similar')));

sc = []; sc.trials = trl; sc.channel = {'all','-A1','-IOL','-LHEOG','-RHEOG'};
dat = ft_selectdata(sc, dat);
rt = [];
rt.begsample = find(dat.time{1} >  cfg.eBOSC.waveseg(1), 1, 'first');
rt.endsample = find(dat.time{1} <= cfg.eBOSC.waveseg(2), 1, 'last');
if any(cellfun(@numel, dat.time) ~= numel(dat.time{1}))
    error('ebosc_subject:triallength', 'Subject %s: trials differ in length', subj);
end
dat = ft_redefinetrial(rt, dat);
if numel(dat.time{1}) <= 2*cfg.eBOSC.trialPad
    error('ebosc_subject:short', 'Subject %s: segment (%d samples) shorter than 2 x trialPad', subj, numel(dat.time{1}));
end
tWin = dat.time{1}(cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad);   % analysis window (~0-1.5 s)

nCh = numel(dat.label); nF = numel(cfg.eBOSC.F); nT = numel(dat.trial); nB = numel(bands.names);
F = cfg.eBOSC.F;

out = struct();
out.label = dat.label; out.F = F; out.bands = bands; out.tWin = tWin([1 end]); out.nWin = numel(tWin);
out.nTrials = nT; out.cfg = cfg;
out.specAll_log10 = nan(nCh, nF); out.specEp_log10 = nan(nCh, nF, nB);
out.bg_log10 = nan(nCh, nF); out.pt_log10 = nan(nCh, nF);
EP = cell(nCh, 1);   % per channel: [chan trial freq on off nPts snrSum]

%% eBOSC per channel
for e = 1:nCh
    % wavelet transform per trial
    B = cell(1,nT);
    for t = 1:nT
        B{t} = BOSC_tf(dat.trial{t}(e,:), F, cfg.eBOSC.fsample, cfg.eBOSC.wavenumber);
    end

    % robust background fit over all trials (peaks in exclFranges removed, FWHM around the empirical peak)
    BG = cell2mat(cellfun(@(b) b(:, cfg.eBOSC.BGpad+1:end-cfg.eBOSC.BGpad), B, 'UniformOutput', false));
    freqKeep = true(size(F));
    for k = 1:size(cfg.eBOSC.exclFranges,1)
        i1 = find(F >= cfg.eBOSC.exclFranges(k,1), 1, 'first');
        i2 = find(F <= cfg.eBOSC.exclFranges(k,2), 1, 'last');
        [~, indPos] = max(mean(BG(i1:i2,:),2));
        indPos = i1 + indPos - 1;   % index of the empirical peak (template had +indPos without -1, i.e. one bin too high)
        LowFreq = F(indPos) - ((2/cfg.eBOSC.wavenumber)*F(indPos))/2;
        UpFreq  = F(indPos) + ((2/cfg.eBOSC.wavenumber)*F(indPos))/2;
        freqKeep(F >= LowFreq & F <= UpFreq) = false;
    end
    pv = eBOSC_bgfit_robust(F(freqKeep), BG(freqKeep,:));
    mp = 10.^(polyval(pv, log10(F)));
    [pt, dt] = BOSC_thresholds(cfg.eBOSC.fsample, cfg.eBOSC.percentile, cfg.eBOSC.ncyc, F, mp);
    clear BG

    sumEp = zeros(nB, nF); nEp = zeros(nB,1); sumAll = zeros(1,nF); nAll = 0;
    rows = zeros(0,7);
    for t = 1:nT
        TF_all = B{t};
        TF_WL  = TF_all(:, cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding);
        TF_Trl = TF_all(:, cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad);

        % detection (power + duration threshold)
        detected = zeros(size(TF_WL));
        for f = 1:nF
            detected(f,:) = BOSC_detect(TF_WL(f,:), pt(f), dt(f), cfg.eBOSC.fsample);
        end

        % episodes (frequency/temporal leakage control), then remove detectedPad
        cfg.eBOSC.npnts = size(detected,2); cfg.eBOSC.pt = pt;
        [detected_abn, episodes] = eBOSC_createEpisodes(TF_WL, detected, cfg); %#ok<ASGLU>
        cfg.eBOSC = rmfield(cfg.eBOSC, {'npnts','pt'});
        ind1 = cfg.eBOSC.detectedPad + 1;
        ind2 = size(detected,2) - cfg.eBOSC.detectedPad;
        rmv = false(size(episodes,1),1);
        for j = 1:size(episodes,1)
            ex = episodes{j,1}(:,2) < ind1 | episodes{j,1}(:,2) > ind2;
            episodes{j,1}(ex,:) = [];
            episodes{j,1}(:,2) = episodes{j,1}(:,2) - ind1 + 1;   % -> column index in TF_Trl
            episodes{j,2}(ex,:) = [];
            if isempty(episodes{j,1}), rmv(j) = true; else, episodes{j,3} = mean(episodes{j,2}(:,1)); end
        end
        episodes(rmv,:) = [];

        % every episode -> one row; band spectra from the episode time points
        idxB = false(nB, size(TF_Trl,2));
        for j = 1:size(episodes,1)
            pts = double(episodes{j,1});
            fr  = episodes{j,3};
            snr = sum(log10(double(episodes{j,2}(:,2))) - log10(mp(pts(:,1)))');
            rows(end+1,:) = [e, t, fr, min(pts(:,2)), max(pts(:,2)), size(pts,1), snr]; %#ok<AGROW>
            b = find(fr >= bands.ranges(:,1) & fr < bands.ranges(:,2), 1);
            if ~isempty(b), idxB(b, pts(:,2)) = true; end
        end
        for b = 1:nB
            if any(idxB(b,:))
                sumEp(b,:) = sumEp(b,:) + sum(log10(TF_Trl(:,idxB(b,:))), 2)';
                nEp(b) = nEp(b) + sum(idxB(b,:));
            end
        end
        sumAll = sumAll + sum(log10(TF_Trl), 2)';
        nAll = nAll + size(TF_Trl,2);
    end

    out.bg_log10(e,:) = log10(mp); out.pt_log10(e,:) = log10(pt);
    out.specAll_log10(e,:) = sumAll / nAll;
    for b = 1:nB
        if nEp(b) > 0, out.specEp_log10(e,:,b) = sumEp(b,:) / nEp(b); end
    end
    EP{e} = rows;
    fprintf('%s ch %d/%d (%s): %d episodes (%s)\n', subj, e, nCh, dat.label{e}, size(rows,1), ...
        strjoin(arrayfun(@(b) sprintf('%s %d', bands.names{b}, sum(rows(:,3) >= bands.ranges(b,1) & rows(:,3) < bands.ranges(b,2))), 1:nB, 'uni', 0), ', '));
end
E = vertcat(EP{:});
out.ep = struct('chan', E(:,1), 'trial', E(:,2), 'freq', E(:,3), 'on', E(:,4), 'off', E(:,5), 'nPts', E(:,6), 'snrSum', E(:,7));

if ~exist(pdat.out,'dir'), mkdir(pdat.out); end
save(fullfile(pdat.out, [subj '_ebosc.mat']), '-struct', 'out');
fprintf('Saved eBOSC episodes for %s (%d trials, %d channels, %d episodes) -> %s\n', subj, nT, nCh, size(E,1), pdat.out);
end

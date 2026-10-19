function rsa_subject_level(subj, contrastType, freqLo, freqHi, pdat, bgFitMode)
% Consolidated per-subject RSA similarity computation.
%
% Replaces Final_RSA/dcn-tutorial-master/self_rsa_vsBL.m (Ispec) and
% within_vs_betw_cat_rsa_vsBL.m (Catspec) with one parameterized
% implementation, so a new contrast type or frequency band is a config
% change, not a new script.
%
%   subj          : subject ID (numeric or string), e.g. 3109
%   contrastType  : 'Ispec'   -> within-item vs within-cat similarity
%                   'IspecPaired' -> (2026-10-19, user) Isim as 'Ispec' plus
%                                an ITEM-SPECIFIC Csim per Isim item X:
%                                presentation 1 of X vs presentation 2 of
%                                the other exemplar (1<->2) of X's category
%                                (partner any response). Both files have
%                                the same row order (one row per item, X's
%                                presentation 1 always the FIRST trial, so
%                                the time axes match) -> paired 1st-level
%                                test (depsamplesT, rsa_group_level_perAge).
%                                Own folder IspecPaired_<band>.
%                   'IspecAll'-> as 'Ispec', but WITHOUT the response
%                                filter (2026-10-18, user): ALL encoding
%                                trials of exemplars 1-2, presentations
%                                1-2, whatever the later response
%                                ('old', 'similar', 'new' or not tested);
%                                only the four experimental conditions
%                                (no target-detection or filler trials).
%                                Own output folder IspecAll_<band>, so
%                                the Ispec results are not touched.
%                   'Bsim'    -> BETWEEN-category control rows (exemplar 1 of A
%                                vs exemplar 2 of every other category, averaged), saved
%                                in the Ispec folder as <ID>_between-cat_sim_ST.mat
%                   'Catspec' -> within-cat  vs between-cat similarity
%   freqLo,freqHi : analysis sub-band in Hz, e.g. [30 100], [30 50],
%                   [60 100], [4 8]
%   pdat          : struct with fields .rsa (output, all variants); for
%                   theta/low-freq (freqHi<=20) also .tfr, .tfrBL (paths
%                   to the original pre-made files, UNCHANGED, same
%                   convention as config.pdat.* in the original Final_RSA
%                   scripts); for gamma (freqHi>20) instead .rawEnc
%                   (preprocessed Enc_allcond data) and .tfrOut (cache
%                   dir for the once-per-subject raw decomposition -- see
%                   load_bosc_correct_restrict below, 2026-09-04 redesign)
%   bgFitMode     : (optional, default 'wide') 'wide' or 'legacy' --
%                   see DESIGN DECISION below. Default is 'wide' for
%                   every variant EXCEPT Ispec broadband (30-100Hz),
%                   which submit_rsa_all_variants.sh runs as 'legacy'
%                   specifically so it reproduces the original pipeline's
%                   numbers and can serve as a validation check.
%
% DESIGN DECISION: two ways to background-correct.
%   'wide'   (new, default): BOSC is fit on the FULL taper-family block
%            -- 2-20Hz (Hanning) for bands within that range, 25-125Hz
%            (DPSS) for all gamma sub-bands -- and ONLY THEN is the
%            corrected data restricted down to [freqLo freqHi]. Needed
%            for narrow sub-bands (theta has only 3 frequency bins in the
%            low-freq TFR -- far too few to fit a stable background
%            regression on their own) and keeps every variant's
%            background fit consistent.
%   'legacy' (original self_rsa_vsBL.m behavior): BOSC is fit directly on
%            the already-restricted [freqLo freqHi] window. Fine for the
%            original single 30-100Hz broadband-gamma analysis (15 bins
%            is enough for a stable fit) but not recommended for the
%            narrower new sub-bands. Use this ONLY for the Ispec
%            broadband variant, to exactly reproduce and validate against
%            the existing published numbers.
%
% Output: saves <subj>_<nameA>_sim_ST.mat and <subj>_<nameB>_sim_ST.mat
% (same 'sim_mat' struct format as the original scripts) into
% pdat.rsa/<contrastType>_<freqLo>to<freqHi>Hz/

if nargin < 6 || isempty(bgFitMode), bgFitMode = 'wide'; end
subj = num2str(subj);

if strcmp(bgFitMode, 'legacy')
    bgFitBand = [freqLo freqHi]; % fit directly on the analysis band (original behavior)
elseif freqHi <= 20
    bgFitBand = [2 20];   % 'wide' mode, Hanning low-freq block
else
    bgFitBand = [25 125]; % 'wide' mode, DPSS high-freq block
end

% Baseline-correction METHOD (2026-09-29: back to the ORIGINAL method of
% self_rsa_vsBL.m for ALL bands): trial B's baseline spectrum is averaged
% over baseline time first, then correlated once with trial A's
% activation spectrum at each time point (tf2new = repmat(mean(tf2,2),...)).
% The "new method" (correlate with every baseline time point, then
% average the correlations; used for gamma 2026-09-XX to 09-28) is still
% in correlate_pairs_and_package, but switched off.
useNewBaselineMethod = false;

%% Load + background-correct on bgFitBand, then restrict to the analysis band
%  (identical for both modes when bgFitMode='legacy', since bgFitBand ==
%  the analysis band there -- restriction after fitting is then a no-op)

% 2026-10-19 (rerun from scratch, user): the pre-made TFR files are no
% longer on the grid -> gamma bands are decomposed from PreprocData with the
% ORIGINAL settings (rawGamma path, see load_bosc_correct_restrict).
[Data_logpow, Data_logpow_bl] = load_bosc_correct_restrict( ...
    subj, pdat, bgFitBand, [freqLo freqHi]);

%% Pairing logic + similarity computation -- the only part that differs by contrast type

switch contrastType
    case 'Ispec'
        [sim_matA, sim_matB, nameA, nameB] = pairs_within_item_vs_cat(Data_logpow, Data_logpow_bl, useNewBaselineMethod);
    case 'IspecAll'
        [sim_matA, sim_matB, nameA, nameB] = pairs_within_item_vs_cat(Data_logpow, Data_logpow_bl, useNewBaselineMethod, false);
    case 'IspecPaired'
        [sim_matA, sim_matB, nameA, nameB] = pairs_item_vs_itemcat(Data_logpow, Data_logpow_bl, useNewBaselineMethod);
    case 'Catspec'
        [sim_matA, sim_matB, nameA, nameB] = pairs_within_vs_between_cat(Data_logpow, Data_logpow_bl, useNewBaselineMethod);
    case 'Bsim'
        % between-category CONTROL pairs for the alpha x pair-type analysis
        % (2026-10-05) -- saved next to the Ispec files (same variant
        % folder), one file only
        sim_matA = pairs_between_cat_allavg(Data_logpow, Data_logpow_bl, useNewBaselineMethod);
        nameA = 'between-cat'; sim_matB = []; nameB = '';
    otherwise
        error('rsa_subject_level:contrastType', 'Unknown contrastType "%s" (expected Ispec, IspecAll, IspecPaired, Catspec or Bsim)', contrastType)
end

%% Save

bandLabel = sprintf('%dto%dHz', freqLo, freqHi);
outType = contrastType; if strcmp(contrastType, 'Bsim'), outType = 'Ispec'; end % Bsim pairs belong to the Ispec analysis; IspecAll keeps its own folder
outDir = fullfile(pdat.rsa, sprintf('%s_%s', outType, bandLabel));
if ~exist(outDir,'dir'), mkdir(outDir); end

sim_mat = sim_matA; %#ok<NASGU>
save(fullfile(outDir, [subj '_' nameA '_sim_ST.mat']), 'sim_mat', '-v7.3')
if ~isempty(sim_matB)
    sim_mat = sim_matB; %#ok<NASGU>
    save(fullfile(outDir, [subj '_' nameB '_sim_ST.mat']), 'sim_mat', '-v7.3')
end

fprintf('Saved subject %s, %s, %d-%dHz -> %s\n', subj, contrastType, freqLo, freqHi, outDir);

end


%% =========================================================
%% LOCAL FUNCTIONS
%% =========================================================

function [Data_logpow, Data_logpow_bl] = load_bosc_correct_restrict(subj, pdat, bgFitBand, analysisBand)
% 2026-09-29: ALL bands (incl. gamma) load the ORIGINAL pre-made TFR files
% (TFRdata/Time-frequencyEEGdata/<ID>_Enc_tfr_2_125Hz_ST.mat and
% TFRdat_BL/<ID>_Enc_tfrBL_2_125Hz_ST.mat) -- with bgFitMode='legacy' this
% is line-for-line self_rsa_vsBL.m (latency [0 1.7] / [-1 -0.3],
% frequency [30 100], BOSC, log10). The raw-decomposition path
% (load_bosc_correct_restrict_rawGamma, 2026-09-04 to 09-28) is no longer
% called but kept for reference.
% 2026-10-19 (rerun from scratch, user decision): the pre-made TFR files
% are no longer on the grid, so ALL gamma bands (analysis band above 20 Hz:
% Ispec, IspecAll, Bsim, Catspec) use the rawGamma path -- decomposition
% from PreprocData/<ID>_Enc_allcond.mat with the ORIGINAL settings (DPSS,
% 25:5:125 Hz, 0.4 s, +/-10 Hz = 7 tapers, 8 ms, pad 6; cached in
% pdat.tfrOut). Activation [0 1.7] / baseline [-1 -0.3], BOSC on the
% analysis band ('legacy'), log10 -- as before. Low-frequency bands (not
% used in the manuscript) would still need the pre-made files.
if analysisBand(2) > 20
    [Data_logpow, Data_logpow_bl] = load_bosc_correct_restrict_rawGamma(subj, pdat, bgFitBand, analysisBand);
else
    [Data_logpow, Data_logpow_bl] = load_bosc_correct_restrict_legacyFile(subj, pdat, bgFitBand, analysisBand);
end
end


function [Data_logpow, Data_logpow_bl] = load_bosc_correct_restrict_legacyFile(subj, pdat, bgFitBand, analysisBand)
% UNCHANGED theta/low-freq path -- load full single-trial TFR + its
% pre-stim baseline TFR from the original pre-made files, BOSC-correct
% each on bgFitBand, then restrict the corrected data to analysisBand.
% Reuses subtract_mean_noise_spectrum.m unchanged (same tool referenced
% elsewhere as mwb_subtract_mean_noise_spectrum.m).

% --- Activation (post-stimulus) data ---
disp('loading activation TFR...')
load(fullfile(pdat.tfr, [subj '_Enc_tfr_2_125Hz_ST.mat']), 'tfr')

cfg = [];
cfg.latency = [0 1.7];
cfg.frequency = bgFitBand;
tfr = ft_selectdata(cfg, tfr);

cfg = []; cfg.parameter = 'powspctrm'; cfg.operation = 'log10';
warning('off', 'stats:statrobustfit:IterationLimit')
disp('running BOSC (activation)...')
Data_c = subtract_mean_noise_spectrum(cfg, tfr);
clear tfr

ft_warning off 'FieldTrip:getdimord:warning_dimord_could_not_be_determined'
Data_logpow_wide = ft_math(cfg, Data_c);
Data_logpow_wide.trialinfo2 = Data_c.trialinfo2;
clear Data_c

cfg = []; cfg.frequency = analysisBand;
Data_logpow = ft_selectdata(cfg, Data_logpow_wide);
Data_logpow.trialinfo2 = Data_logpow_wide.trialinfo2;
clear Data_logpow_wide

% --- Baseline (pre-stimulus) data ---
disp('loading baseline TFR...')
load(fullfile(pdat.tfrBL, [subj '_Enc_tfrBL_2_125Hz_ST.mat']), 'tfr')

cfg = [];
cfg.latency = [-1 -0.3];
cfg.frequency = bgFitBand;
tfr = ft_selectdata(cfg, tfr);

cfg = []; cfg.parameter = 'powspctrm'; cfg.operation = 'log10';
warning('off', 'stats:statrobustfit:IterationLimit')
disp('running BOSC (baseline)...')
Data_c = subtract_mean_noise_spectrum(cfg, tfr);
clear tfr

Data_logpow_bl_wide = ft_math(cfg, Data_c);
Data_logpow_bl_wide.trialinfo2 = Data_c.trialinfo2;
clear Data_c

cfg = []; cfg.frequency = analysisBand;
Data_logpow_bl = ft_selectdata(cfg, Data_logpow_bl_wide);
Data_logpow_bl.trialinfo2 = Data_logpow_bl_wide.trialinfo2;
clear Data_logpow_bl_wide

end


function [Data_logpow, Data_logpow_bl] = load_bosc_correct_restrict_rawGamma(subj, pdat, bgFitBand, analysisBand)
% USED for all gamma bands since 2026-10-19 (pre-made TFR files no longer
% on the grid; see load_bosc_correct_restrict). 2026-09-29 to 10-18 gamma
% read the pre-made files instead.
%
% NEW (2026-09-04): computes the gamma single-trial TFR directly from the
% preprocessed Enc_allcond data instead of reading a pre-made file (the
% original driver script that produced *_Enc_tfr_2_125Hz_ST.mat /
% *_Enc_tfrBL_2_125Hz_ST.mat could not be located -- see session notes).
%
%   DECOMPOSITION = ORIGINAL SETTINGS (restored 2026-09-29): the high-
%   frequency half of the original driver script -- mtmconvol, DPSS,
%   foi = 25:5:125 Hz, t_ftimwin = 0.4 s, tapsmofrq = 10 Hz (+/-10 Hz,
%   7 tapers), pad = 6, toi step 8 ms on the original -2:0.008:3.5 grid.
%   (A 5 Hz / 3-taper, 30:2:100 Hz, 20 ms variant was used 2026-09-04 to
%   09-28; reverted to reproduce the original pipeline.) The original's
%   low-frequency half (Hanning, 2:2:20 Hz, 0.5 s) is NOT computed here:
%   only the gamma band is used from this path (theta/low-freq uses the
%   pre-made files, see load_bosc_correct_restrict_legacyFile).
%   toi: only the part of the original grid covering the baseline and
%   activation windows ([-1 -0.3], [0 1.7]) is computed -- mtmconvol
%   estimates each time point independently, so these values are
%   identical to the full-grid ones; it just avoids ~2x cache size.
%   With bgFitMode='legacy' and a 30-100 Hz analysis band, BOSC is fit on
%   the 15 bins 30:5:100 Hz, as in the original pipeline.
%
%   Windows [0 1.7] / [-1 -0.3] (2026-09-05, back to ORIGINAL): briefly
%   changed to [0 1.5]/[-1.5 0] (equal length) to support a shape-matched
%   cell-by-cell baseline-correction method, but reverted after the
%   shortened activation window visibly truncated the cluster (real
%   post-offset signal extends toward ~1.7s, consistent with the alpha
%   SME cluster's own extent to 1.68s -- losing it was a real cost, not
%   just a theoretical one). Baseline reverted to the original, narrower
%   [-1 -0.3] too, not just for symmetry with the activation window, but
%   because if post-stimulus representations linger past nominal offset
%   (why we need the wider activation window), the same logic suggests
%   the last ~300ms before onset may not be a clean pre-stimulus
%   reference either -- likely why the original pipeline excluded it in
%   the first place. This is fully compatible with the new baseline-
%   correction method (see correlate_pairs_and_package) since that method
%   no longer requires equal-length windows.
%
% Cached to pdat.tfrOut per subject so this expensive step only runs
% once, regardless of how many times downstream steps are rerun.
%
% ATOMIC WRITE (2026-09-05, bugfix): this cache file is shared across
% ALL frequency-band variants using bgFitMode='legacy' within 25-125Hz
% (broadband, low-gamma, high-gamma, etc.) -- if their subject-level
% array jobs are submitted around the same time, before the cache exists,
% multiple tasks for the SAME subject can race to check-then-write it
% simultaneously, corrupting the HDF5 (-v7.3) file (observed: several
% subjects failing with the same "attribute already exists" error across
% concurrently-running broadband/low-gamma/high-gamma jobs). Fixed by
% writing to a uniquely-named temp file first, then atomically renaming
% it into place -- concurrent writers can no longer interleave/corrupt
% each other; worst case is redundant computation the first time a
% subject's cache is populated, which is harmless.

% New cache name (2026-09-29): the old '_30to100Hz' caches hold the 5 Hz /
% 20 ms variant and must NOT be reused for the original settings.
cacheFile = fullfile(pdat.tfrOut, [subj '_Enc_gammaTFR_25to125Hz_orig.mat']);

if exist(cacheFile, 'file')
    disp('loading cached raw gamma TFR...')
    load(cacheFile, 'tfr')
else
    disp('loading preprocessed data (Enc_allcond)...')
    load(fullfile(pdat.rawEnc, [subj '_Enc_allcond.mat']), 'data')

    cfg = [];
    cfg.output     = 'pow';
    cfg.channel    = {'all', '-IOL','-RHEOG','-LHEOG','-A1'};
    cfg.keeptrials = 'yes';
    cfg.pad        = 6;
    cfg.method     = 'mtmconvol';
    cfg.taper      = 'dpss';                             % slepian tapers
    cfg.foi        = 25:5:125;
    cfg.t_ftimwin  = ones(size(cfg.foi)) .* 0.4;         % 400 ms window per freq
    cfg.tapsmofrq  = ones(size(cfg.foi)) .* 10;          % +/-10 Hz smoothing (7 tapers)
    toiFull        = -2:0.008:3.5;                       % original time grid
    cfg.toi        = toiFull(toiFull >= -1.01 & toiFull <= 1.71); % covers [-1 -0.3] and [0 1.7]

    disp('running TFR decomposition (original settings: DPSS, 25:5:125 Hz, 0.4 s, +/-10 Hz, 8 ms)...')
    tfr = ft_freqanalysis(cfg, data);
    tfr.trialinfo2 = data.trialinfo2;
    clear data

    if ~exist(pdat.tfrOut, 'dir'), mkdir(pdat.tfrOut); end
    tmpFile = [tempname(pdat.tfrOut) '.mat'];
    save(tmpFile, 'tfr', '-v7.3')
    movefile(tmpFile, cacheFile); % atomic rename -- avoids concurrent-writer corruption
end

% --- Activation (post-stimulus) window ---
cfg = [];
cfg.latency = [0 1.7];
cfg.frequency = bgFitBand;
tfrAct = ft_selectdata(cfg, tfr);
tfrAct.trialinfo2 = tfr.trialinfo2;

cfg = []; cfg.parameter = 'powspctrm'; cfg.operation = 'log10';
warning('off', 'stats:statrobustfit:IterationLimit')
disp('running BOSC (activation)...')
Data_c = subtract_mean_noise_spectrum(cfg, tfrAct);
clear tfrAct

ft_warning off 'FieldTrip:getdimord:warning_dimord_could_not_be_determined'
Data_logpow_wide = ft_math(cfg, Data_c);
Data_logpow_wide.trialinfo2 = Data_c.trialinfo2;
clear Data_c

cfg = []; cfg.frequency = analysisBand;
Data_logpow = ft_selectdata(cfg, Data_logpow_wide);
Data_logpow.trialinfo2 = Data_logpow_wide.trialinfo2;
clear Data_logpow_wide

% --- Baseline (pre-stimulus) window ---
cfg = [];
cfg.latency = [-1 -0.3];
cfg.frequency = bgFitBand;
tfrBLsel = ft_selectdata(cfg, tfr);
tfrBLsel.trialinfo2 = tfr.trialinfo2;
clear tfr

cfg = []; cfg.parameter = 'powspctrm'; cfg.operation = 'log10';
warning('off', 'stats:statrobustfit:IterationLimit')
disp('running BOSC (baseline)...')
Data_c = subtract_mean_noise_spectrum(cfg, tfrBLsel);
clear tfrBLsel

Data_logpow_bl_wide = ft_math(cfg, Data_c);
Data_logpow_bl_wide.trialinfo2 = Data_c.trialinfo2;
clear Data_c

cfg = []; cfg.frequency = analysisBand;
Data_logpow_bl = ft_selectdata(cfg, Data_logpow_bl_wide);
Data_logpow_bl.trialinfo2 = Data_logpow_bl_wide.trialinfo2;
clear Data_logpow_bl_wide

end


function [sim_matI, sim_matC, nameA, nameB] = pairs_within_item_vs_cat(Data_logpow, Data_logpow_bl, useNewBaselineMethod, respFilter)
% Ispec pairing logic -- ported unchanged from self_rsa_vsBL.m (audited
% and bug-fixed, see session_notes_2026-08-25_alpha_theory.txt Part 13),
% just taking already-corrected Data_logpow/Data_logpow_bl as input
% instead of loading+correcting them itself.
% respFilter (optional, default true): true = only items later called
% 'old'/'similar' (Ispec); false = no response filter (IspecAll, 2026-10-18).
if nargin < 4, respFilter = true; end

nameA = 'within-item';
nameB = 'within-cat';
pres1 = 1;
% EXEMPLARS 1-2 ONLY (final, 2026-09-30): balanced with Csim (exemplar 1
% vs 2) and with alpha_pre2post/PF (exemplar<3).
maxExemplar = 2;

%% Within-item: find stimulus pairs
pairs = [];
for stim1 = 1:size(Data_logpow.powspctrm,1)
    respOK = ~respFilter || ~isempty(find(strcmp(Data_logpow.trialinfo2{stim1,11},'old') | strcmp(Data_logpow.trialinfo2{stim1,11},'similar')));
    if ~ismember(stim1,pairs) && Data_logpow.trialinfo2{stim1,9} == pres1 && Data_logpow.trialinfo2{stim1,8} <= maxExemplar && respOK
        curr_pair = Data_logpow.trialinfo2{stim1,7};
        if strcmp(Data_logpow.trialinfo2(stim1,4),'normal') || strcmp(Data_logpow.trialinfo2(stim1,4),'high-exemplar')
            if find(strcmp(Data_logpow.trialinfo2(:,7),curr_pair) & cell2mat(Data_logpow.trialinfo2(:,9)) == 2)
                pairs(stim1,2) = find(strcmp(Data_logpow.trialinfo2(:,7),curr_pair) & cell2mat(Data_logpow.trialinfo2(:,9)) == 2);
                pairs(stim1,1) = stim1;
            end
        elseif strcmp(Data_logpow.trialinfo2(stim1,4),'high-repetition+exemplar') || strcmp(Data_logpow.trialinfo2(stim1,4),'high-repetition')
            if find(strcmp(Data_logpow.trialinfo2(:,7),curr_pair) & cell2mat(Data_logpow.trialinfo2(:,9)) == 2)
                pairs(stim1,1) = stim1;
                pairs(stim1,2) = find(strcmp(Data_logpow.trialinfo2(:,7),curr_pair) & cell2mat(Data_logpow.trialinfo2(:,9)) == 2);
            end
        end
    end
end
pairs(~all(pairs,2),:) = [];

sim_matI = correlate_pairs_and_package(Data_logpow, Data_logpow_bl, pairs, 'Computing within item similarity', useNewBaselineMethod);

%% Within-category similarity (exemplar 1 vs exemplar 2 specifically --
%  see pairs_within_cat_exemplar1v2 for the full design rationale)
sim_matC = pairs_within_cat_exemplar1v2(Data_logpow, Data_logpow_bl, useNewBaselineMethod, respFilter);

end


function [sim_matW, sim_matB, nameA, nameB] = pairs_within_vs_between_cat(Data_logpow, Data_logpow_bl, useNewBaselineMethod)
% Catspec pairing logic. Within-cat now delegates to
% pairs_within_cat_exemplar1v2 (finalized design, see that function's
% header -- no longer "ported unchanged", substantially redesigned
% 2026-08-25). Between-cat structure is still ported from
% within_vs_betw_cat_rsa_vsBL.m (bug-fixed 'omitnan' calls, see session
% notes Part 13) with one further fix: the partner's response filter
% widened from 'old'-only to 'old'-or-'similar' (see comment below).

nameA = 'within-cat';
nameB = 'between-cat';

%% Within-category similarity (exemplar 1 vs exemplar 2 -- IDENTICAL
%  computation to Ispec's cat sim, confirmed with user 2026-08-25.
%  Shared function so the two can't drift apart the way they already did
%  once (the missing 'omitnan' inconsistency found earlier).
sim_matW = pairs_within_cat_exemplar1v2(Data_logpow, Data_logpow_bl, useNewBaselineMethod);

%% Between-category similarity (exemplar 1 vs exemplar 2 of a DIFFERENT
%  category -- BOTH sides require a valid 'old'/'similar' response, no
%  relaxation here, unlike within-cat. Rationale: between-cat already
%  averages each reference item against every valid partner across all
%  other categories, so losing one partner to strict filtering barely
%  affects the average -- unlike within-cat's single, all-or-nothing
%  pair, there's no exemplar-loss cost that would justify relaxing this
%  one. Also fixed here: the original code restricted the partner to
%  'old' responses only (not 'old' or 'similar'), asymmetric with the
%  reference item's filter for no stated reason -- made symmetric.)
pairs1 = [];
pairs2 = [];
for stim1 = 1:size(Data_logpow.powspctrm,1)
    if ~ismember(stim1,pairs1) && cell2mat(Data_logpow.trialinfo2(stim1,8)) == 1 && cell2mat(Data_logpow.trialinfo2(stim1,9)) == 1 && ~isempty(find(strcmp(Data_logpow.trialinfo2{stim1,11},'old') | strcmp(Data_logpow.trialinfo2{stim1,11},'similar')))
        curr_pair = Data_logpow.trialinfo2{stim1,5};
        stim2 = find(~strcmp(Data_logpow.trialinfo2(:,5),curr_pair) & cell2mat(Data_logpow.trialinfo2(:,8)) == 2 & cell2mat(Data_logpow.trialinfo2(:,9)) == 1 & (strcmp(Data_logpow.trialinfo2(:,11),'old') | strcmp(Data_logpow.trialinfo2(:,11),'similar')));
        if ~isempty(stim2)
            for s = 1:length(stim2)
                pairs1(end+1,1) = stim1;
                pairs2(end+1,1) = stim2(s);
            end
        end
    end
end
pairsU = unique(pairs1);

corrmat = zeros(size(pairsU,1), size(Data_logpow.powspctrm,2), length(Data_logpow.time), length(Data_logpow.time));
corrmat_bl = zeros(size(corrmat));

fprintf('Computing between category similarity. Correlating pair: ');
for curr_pair = 1:size(pairsU,1)
    fprintf('%3d of %3d', curr_pair, size(pairsU,1))
    ix = find(pairs1 == pairsU(curr_pair));

    corrmat_tmp = zeros(length(ix), size(Data_logpow.powspctrm,2), length(Data_logpow.time), length(Data_logpow.time));
    corrmat_bl_tmp = zeros(size(corrmat_tmp));

    for i = 1:length(ix)
        for chan = 1:size(Data_logpow.powspctrm,2)
            tf1 = squeeze(Data_logpow.powspctrm(pairsU(curr_pair),chan,:,:));
            tf2 = squeeze(Data_logpow.powspctrm(pairs2(ix(i)),chan,:,:));
            corrmat_tmp(i,chan,:,:) = spectral_rsa([],tf1,tf2,Data_logpow.time);

            tf2 = squeeze(Data_logpow_bl.powspctrm(pairs2(ix(i)),chan,:,:));
            if useNewBaselineMethod
                r_raw = spectral_rsa([],tf1,tf2,Data_logpow_bl.time); % [nTimeAct x nTimeBL], any lengths
                corrmat_bl_tmp(i,chan,:,:) = repmat(mean(r_raw,2,'omitnan'),1,size(tf1,2));
            else
                tf2new = repmat(mean(tf2,2,'omitnan'),1,size(tf1,2));
                corrmat_bl_tmp(i,chan,:,:) = spectral_rsa([],tf1,tf2new,Data_logpow.time);
            end
        end
    end
    corrmat(curr_pair,:,:,:) = squeeze(mean(corrmat_tmp,1,'omitnan'));
    corrmat_bl(curr_pair,:,:,:) = squeeze(mean(corrmat_bl_tmp,1,'omitnan'));
    fprintf(repmat('\b',1,10));
end
fprintf(' done \n');

sim_matB = Data_logpow;
sim_matB.orig_trialinfo = Data_logpow.trialinfo;
sim_matB.freq = sim_matB.time; % actually time but called freq (FieldTrip dimord compatibility)
sim_matB = rmfield(sim_matB,'powspctrm');
sim_matB.powspctrm = corrmat;
sim_matB.powspctrmBL = corrmat_bl;
sim_matB.pairs = pairsU;
sim_matB.trialinfo2 = Data_logpow.trialinfo2(pairsU,:);
sim_matB.trialinfo = cell2mat(sim_matB.trialinfo2(:,1:3));
sim_matB.indtrials = Data_logpow.trialinfo2(sort(unique([pairs1;pairs2])),:);

sim_matB.powspctrmOrig = sim_matB.powspctrm;
sim_matB.powspctrm = sim_matB.powspctrm - sim_matB.powspctrmBL;

end


function [sim_matI, sim_matC, nameA, nameB] = pairs_item_vs_itemcat(Data_logpow, Data_logpow_bl, useNewBaselineMethod)
% IspecPaired (2026-10-19, user): one Isim and one item-specific Csim per
% item, row-aligned, for a PAIRED within-person test.
%   Isim(X) = X presentation 1 vs X presentation 2 (as Ispec: exemplars 1-2,
%             X later called 'old'/'similar', four experimental conditions)
%   Csim(X) = X presentation 1 vs presentation 2 of the OTHER exemplar of
%             X's category (exemplar 1<->2; partner any response, as the
%             at-least-one-remembered Csim rule)
% X's presentation 1 is the FIRST trial of both pairs (time1 axis), so the
% two maps of an item share their first trial and their time1 axis; the
% baseline correction uses the second trial's pre-stimulus spectrum, as
% for every pair. Items whose partner trial is missing drop out of BOTH.
nameA = 'within-item';
nameB = 'within-catItem';
ti = Data_logpow.trialinfo2;
pres = cell2mat(ti(:,9)); exem = cell2mat(ti(:,8));
resp = ti(:,11);
valid = cellfun(@(r) ischar(r) && (strcmp(r,'old') || strcmp(r,'similar')), resp);
expConds = {'normal','high-exemplar','high-repetition','high-repetition+exemplar'};
isExp = cellfun(@(c) ischar(c) && any(strcmp(c, expConds)), ti(:,4));
pairsI = zeros(0,2); pairsC = zeros(0,2);
for x1 = find(pres == 1 & exem <= 2 & valid & isExp)'
    x2 = find(strcmp(ti(:,7), ti{x1,7}) & pres == 2, 1);                       % X presentation 2
    y2 = find(strcmp(ti(:,5), ti{x1,5}) & exem == 3 - exem(x1) & pres == 2, 1); % partner presentation 2
    if isempty(x2) || isempty(y2), continue, end
    pairsI(end+1,:) = [x1 x2]; %#ok<AGROW>
    pairsC(end+1,:) = [x1 y2]; %#ok<AGROW>
end
fprintf('IspecPaired: %d items with Isim and item-specific Csim\n', size(pairsI,1));
sim_matI = correlate_pairs_and_package(Data_logpow, Data_logpow_bl, pairsI, 'Computing within item similarity (paired design)', useNewBaselineMethod);
sim_matC = correlate_pairs_and_package(Data_logpow, Data_logpow_bl, pairsC, 'Computing item-specific category similarity', useNewBaselineMethod);
end


function sim_mat = pairs_within_cat_exemplar1v2(Data_logpow, Data_logpow_bl, useNewBaselineMethod, respFilter)
% Within-category pairing. Used identically by BOTH Ispec's cat sim and
% Catspec's within-cat -- these are confirmed to be the same computation,
% hence one shared function rather than two copies that could drift.
% respFilter (optional, default true): false = no response filter
% (IspecAll, 2026-10-18); the category must then belong to one of the four
% experimental conditions (the response filter used to exclude the
% target-detection and filler trials, which were never tested).
if nargin < 4, respFilter = true; end
expConds = {'normal','high-exemplar','high-repetition','high-repetition+exemplar'};
%
% DESIGN:
%
%  1. EXEMPLAR-RESTRICTED, not "any other same-category item": Hat A must
%     be exemplar 1 (trialinfo2 col 8 == 1), Hat B must be exemplar 2
%     (col 8 == 2) SPECIFICALLY. Keeps always exactly one pair per
%     category regardless of catSize, matching item sim's structure
%     (always exactly one pair per item regardless of itemRep) -- keeps
%     the two measures on parallel footing, which matters because item
%     specificity = item sim MINUS cat sim.
%
%  UPDATE 2026-09-26: point 2a below (duos (1,2) AND (3,4), all
%  exemplars in Isim) STILL HOLDS -- constant Isim:Csim ratio (2:1 low
%  catSize, 4:2 high catSize), only independent pairs of neighbouring
%  exemplars (all presentations of exemplar k precede exemplar k+1), no
%  catSize-dependent weighting (all-pairs 1-2,1-3,...,3-4 would give 4:6
%  in high-catSize categories, shared exemplars, very unequal lags).
%  Point 2b (OR-validity) is SUPERSEDED: BOTH
%  exemplars must now have an 'old'/'similar' response. Rule used by
%  every analysis (Isim, Csim, between-cat, alpha_pre2post, PF,
%  single-trial): no 'new' responses anywhere -- forgotten items carry no
%  memory signal, and Isim pairs never contained them (both presentations
%  share the item's one test response), so OR-validity let Csim contain
%  forgotten exemplars that Isim never did. 'new' responses to old items
%  are rare (~5% in both age groups), so the trial cost is small.
%
%  2. EXEMPLAR-DUO EXTENSION + OR-VALIDITY (2026-09-XX). Two
%     combined fixes, replacing the exemplar-1-vs-2-only, single-pair-
%     per-category design used until now:
%
%     a) EXEMPLAR DUOS (1,2) AND (3,4), not just (1,2): the original
%        exemplar-1-vs-2-only restriction (itself a deliberate, pre-
%        existing fix predating this file's recent history, to stop
%        high-catSize categories silently generating more pairs than
%        low-catSize ones -- see the original design note this replaces)
%        never extended that same fix to Isim. Isim scales with exemplar
%        count (up to 4 pairs/category for high-catSize categories, 1 per
%        exemplar) while Csim was capped at exactly 1/category regardless
%        of size -- so high-catSize categories were still over-weighted
%        in Isim relative to Csim, just via a different mechanism than
%        the one already fixed. Pairing exemplar 3 against exemplar 4 as
%        a SECOND, independent duo (for categories where they exist)
%        restores a constant Isim:Csim ratio across catSize levels
%        (2 Isim : 1 Csim pair for low-catSize, 4:2 for high-catSize)
%        instead of a catSize-dependent one (2:1 vs 4:1) -- what actually
%        matters for the validity of the Isim-vs-Csim contrast is that
%        ratio being constant, not that Isim and Csim have equal absolute
%        counts. Low-catSize categories are unaffected (duo (3,4) simply
%        doesn't exist for them, so only duo (1,2) ever contributes).
%
%     b) OR-VALIDITY per duo (either member valid, not a fixed anchor):
%        a fixed-anchor rule (only exemplar 1, or only exemplar 3, gates
%        its duo) discards a duo entirely whenever the specific anchor
%        happens to be forgotten, even if its partner has a perfectly
%        valid response -- purely because of which position it occupies,
%        not for any principled reason. Critically, this doesn't reduce
%        contamination relative to OR: whenever a fixed-anchor rule DOES
%        include a duo, it already allows one confirmed-valid pattern
%        paired with one UNCHECKED partner (whatever that partner's
%        status happens to be) -- OR-validity has the exact same
%        contamination profile in every included duo, it just also
%        rescues the mirror-image case (partner valid, anchor forgotten)
%        that fixed-anchor discards for no principled reason. OR is a
%        strict superset of fixed-anchor's included duos with no new
%        failure mode, so there's no real cost to preferring it.
%        (An even stricter "both members valid" rule was also tried and
%        rejected empirically -- collapsed per-subject within-cat trial
%        counts to as low as 5 for some subjects, unusable given how thin
%        this design already is.)

pairs = [];
categories = unique(Data_logpow.trialinfo2(:,5));
exemplarDuos = [1 2]; % exemplars 1-2 only (final, 2026-09-30) -- duo (3,4) dropped, see header
for c = 1:length(categories)
    curr_cat = categories{c};
    for d = 1:size(exemplarDuos,1)
        exA = exemplarDuos(d,1); exB = exemplarDuos(d,2);
        idxA = find(strcmp(Data_logpow.trialinfo2(:,5),curr_cat) & cell2mat(Data_logpow.trialinfo2(:,8)) == exA & cell2mat(Data_logpow.trialinfo2(:,9)) == 1);
        idxB = find(strcmp(Data_logpow.trialinfo2(:,5),curr_cat) & cell2mat(Data_logpow.trialinfo2(:,8)) == exB & cell2mat(Data_logpow.trialinfo2(:,9)) == 1);
        if isempty(idxA) || isempty(idxB)
            continue % this duo doesn't exist in this category (e.g. low-catSize categories have no exemplar 3/4)
        end
        stimA = idxA(1); stimB = idxB(1);
        validA = strcmp(Data_logpow.trialinfo2{stimA,11},'old') || strcmp(Data_logpow.trialinfo2{stimA,11},'similar');
        validB = strcmp(Data_logpow.trialinfo2{stimB,11},'old') || strcmp(Data_logpow.trialinfo2{stimB,11},'similar');
        if ~respFilter % IspecAll: any response or untested, experimental conditions only
            validA = any(strcmp(Data_logpow.trialinfo2{stimA,4}, expConds));
            validB = any(strcmp(Data_logpow.trialinfo2{stimB,4}, expConds));
        end
        % 2026-10-19 (user): AT LEAST ONE exemplar remembered ('old'/'similar');
        % the partner may be 'new' or untested. Replaces the 2026-09-26
        % BOTH-rule, whose rationale assumed 'new' responses to old items
        % are rare (~5%) -- they are ~20% (YA .23, OA .20), so the BOTH-rule
        % dropped about a third of the categories. Csim is a property of
        % encoding (overlap of two exemplars' patterns), so a remembered
        % item keeps its category similarity even if the partner was later
        % forgotten. IspecAll (respFilter false): both must be experimental.
        if (respFilter && (validA || validB)) || (~respFilter && validA && validB)
            pairs(end+1,:) = [stimA, stimB]; %#ok<AGROW>
        end
    end
end

sim_mat = correlate_pairs_and_package(Data_logpow, Data_logpow_bl, pairs, 'Computing within category similarity (exemplar 1 vs 2)', useNewBaselineMethod);

end


function sim_mat = pairs_between_cat_allavg(Data_logpow, Data_logpow_bl, useNewBaselineMethod)
% BETWEEN-CATEGORY control rows (2026-10-05): per category A, exemplar 1
% (presentation 1) of A correlated with exemplar 2 (presentation 1) of
% EVERY other category, then AVERAGED over those partners -> one row per
% category, like Csim (close and distant categories both included, no
% superordinate grouping needed). Same trial types and rules as Csim
% (pairs_within_cat_exemplar1v2): only categories A with a Csim pair
% (2026-10-19: at least one of exemplar 1 and 2 'old'/'similar');
% reference = exemplar 1 of A, partners = exemplar 2 of every other
% category that has a Csim pair (any response, as in Csim). Same similarity and baseline correction as
% correlate_pairs_and_package (per partner, then averaged). No lag
% matching: temporal distance is handled by the baseline correction, as
% for Isim vs Csim. Fields as correlate_pairs_and_package; pairs(:,2) =
% NaN (many partners), partnerIdx{row} = the partner trials.
ti = Data_logpow.trialinfo2;
valid = strcmp(ti(:,11),'old') | strcmp(ti(:,11),'similar');
pres = cell2mat(ti(:,9)); exem = cell2mat(ti(:,8));
cats = unique(ti(:,5));
A1 = nan(numel(cats),1); A2 = nan(numel(cats),1);
for c = 1:numel(cats)
    i1 = find(strcmp(ti(:,5), cats{c}) & exem == 1 & pres == 1, 1);
    i2 = find(strcmp(ti(:,5), cats{c}) & exem == 2 & pres == 1, 1);
    if ~isempty(i1) && ~isempty(i2) && (valid(i1) || valid(i2)) % = Csim rule (at least one remembered, 2026-10-19)
        A1(c) = i1; A2(c) = i2;
    end
end
refC = find(~isnan(A1) & ~isnan(A2));          % = categories with a Csim pair
nCh = size(Data_logpow.powspctrm,2); nTm = length(Data_logpow.time);
corrmat = zeros(numel(refC), nCh, nTm, nTm); corrmat_bl = zeros(size(corrmat));
partnerIdx = cell(numel(refC),1);
fprintf('Computing between category similarity (all other categories, averaged). Category: ');
for r = 1:numel(refC)
    fprintf('%3d of %3d', r, numel(refC));
    c = refC(r);
    partners = A2(~isnan(A2) & (1:numel(cats))' ~= c);
    partnerIdx{r} = partners;
    acc = zeros(nCh, nTm, nTm); acc_bl = zeros(nCh, nTm, nTm);
    for k = 1:numel(partners)
        for chan = 1:nCh
            tf1 = squeeze(Data_logpow.powspctrm(A1(c),chan,:,:));
            tf2 = squeeze(Data_logpow.powspctrm(partners(k),chan,:,:));
            acc(chan,:,:) = squeeze(acc(chan,:,:)) + spectral_rsa([],tf1,tf2,Data_logpow.time);
            tf2 = squeeze(Data_logpow_bl.powspctrm(partners(k),chan,:,:));
            if useNewBaselineMethod
                r_raw = spectral_rsa([],tf1,tf2,Data_logpow_bl.time);
                acc_bl(chan,:,:) = squeeze(acc_bl(chan,:,:)) + repmat(mean(r_raw,2,'omitnan'),1,size(tf1,2));
            else
                tf2new = repmat(mean(tf2,2,'omitnan'),1,size(tf1,2));
                acc_bl(chan,:,:) = squeeze(acc_bl(chan,:,:)) + spectral_rsa([],tf1,tf2new,Data_logpow.time);
            end
        end
    end
    corrmat(r,:,:,:) = acc / numel(partners);
    corrmat_bl(r,:,:,:) = acc_bl / numel(partners);
    fprintf(repmat('\b',1,10));
end
fprintf(' done \n');

sim_mat = Data_logpow;
sim_mat.orig_trialinfo = Data_logpow.trialinfo;
sim_mat.freq = sim_mat.time; % actually time but called freq (FieldTrip dimord compatibility)
sim_mat = rmfield(sim_mat,'powspctrm');
sim_mat.powspctrmOrig = corrmat;
sim_mat.powspctrmBL   = corrmat_bl;
sim_mat.powspctrm     = corrmat - corrmat_bl;
sim_mat.pairs      = [A1(refC), nan(numel(refC),1)];
sim_mat.partnerIdx = partnerIdx;
sim_mat.nPartners  = cellfun(@numel, partnerIdx);
sim_mat.trialinfo2 = Data_logpow.trialinfo2(A1(refC),:);
sim_mat.trialinfo  = cell2mat(sim_mat.trialinfo2(:,1:3));
sim_mat.indtrials  = Data_logpow.trialinfo2(sort(unique([A1(refC); vertcat(partnerIdx{:})])),:);
fprintf('Between-category rows: %d categories, median %d partners each\n', numel(refC), median(sim_mat.nPartners));
end


function sim_mat = correlate_pairs_and_package(Data_logpow, Data_logpow_bl, pairs, progressLabel, useNewBaselineMethod)
% Shared correlate-and-package step used by both Ispec pairing types
% (within-item, within-cat) -- ported from self_rsa_vsBL.m.
%
% BASELINE-CORRECTION METHOD (2026-09-05, revised): useNewBaselineMethod
% =true (gamma) correlates activation against each RAW baseline
% timepoint individually (avoids the attenuation-due-to-unreliability
% concern of denoising the baseline before correlating -- see session
% notes Part 24 / diagnostics/test_baseline_correction_onesubj.m, "NEW
% method"), THEN averages the resulting CORRELATION VALUES across the
% baseline-time axis and broadcasts that back out to match the
% activation-vs-activation matrix's shape. This "correlate-then-average"
% ordering keeps the new method's statistical benefit without requiring
% activation/baseline windows to be the same length (an earlier version
% of this DID require equal lengths, forcing a shortened activation
% window that visibly truncated the cluster -- reverted; this version
% works with the original, unequal-length [0 1.7]/[-1 -0.3] windows).
% useNewBaselineMethod=false (theta) keeps the ORIGINAL method (average
% the raw baseline SIGNAL over time first, correlate once).

corrmat = zeros(size(pairs,1), size(Data_logpow.powspctrm,2), length(Data_logpow.time), length(Data_logpow.time));
corrmat_bl = zeros(size(corrmat));

fprintf('%s. Correlating pair: ', progressLabel);
for curr_pair = 1:size(pairs,1)
    fprintf('%3d of %3d', curr_pair, size(pairs,1))
    for chan = 1:size(Data_logpow.powspctrm,2)
        tf1 = squeeze(Data_logpow.powspctrm(pairs(curr_pair,1),chan,:,:));
        tf2 = squeeze(Data_logpow.powspctrm(pairs(curr_pair,2),chan,:,:));
        corrmat(curr_pair,chan,:,:) = spectral_rsa([],tf1,tf2,Data_logpow.time);

        tf2 = squeeze(Data_logpow_bl.powspctrm(pairs(curr_pair,2),chan,:,:));
        if useNewBaselineMethod
            r_raw = spectral_rsa([],tf1,tf2,Data_logpow_bl.time); % [nTimeAct x nTimeBL], any lengths
            corrmat_bl(curr_pair,chan,:,:) = repmat(mean(r_raw,2,'omitnan'),1,size(tf1,2));
        else
            tf2new = repmat(mean(tf2,2,'omitnan'),1,size(tf1,2));
            corrmat_bl(curr_pair,chan,:,:) = spectral_rsa([],tf1,tf2new,Data_logpow.time);
        end
    end
    fprintf(repmat('\b',1,10));
end
fprintf(' done \n');

sim_mat = Data_logpow;
sim_mat.orig_trialinfo = Data_logpow.trialinfo;
sim_mat.freq = sim_mat.time; % actually time but called freq (FieldTrip dimord compatibility)
sim_mat = rmfield(sim_mat,'powspctrm');
sim_mat.powspctrm = corrmat;
sim_mat.powspctrmBL = corrmat_bl;
sim_mat.pairs = pairs;
sim_mat.trialinfo2 = Data_logpow.trialinfo2(pairs(:,1),:);
sim_mat.itemLabel2 = Data_logpow.trialinfo2(pairs(:,2),7); % stim2's item label, per pair -- needed
% by rsa_singleTrial_extract.m's cat-sim alpha matching (averages stim1's
% and stim2's own presentation-1 power, since a within-cat pair compares
% two DIFFERENT items -- added 2026-08-25, not in the original scripts,
% which only ever kept stim1's info here).
sim_mat.trialinfo2_B = Data_logpow.trialinfo2(pairs(:,2),:); % stim2's FULL trialinfo2 row
% (2026-09-XX) -- superset of itemLabel2, added so downstream steps (the
% new category_level_paired_ttest in rsa_group_level_perAge.m) can check
% Hat B's own response for the symmetric-response-filter design, without
% needing to re-load raw trial data.
sim_mat.trialinfo = cell2mat(sim_mat.trialinfo2(:,1:3));
sim_mat.indtrials = Data_logpow.trialinfo2(sort([pairs(:,1);pairs(:,2)]),:);

sim_mat.powspctrmOrig = sim_mat.powspctrm;
sim_mat.powspctrm = sim_mat.powspctrmOrig - sim_mat.powspctrmBL;

end

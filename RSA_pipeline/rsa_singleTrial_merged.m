function rsa_singleTrial_merged(subj, freqLo, freqHi, pdat)
% Per-subject single-trial extraction for the Ispec variant: similarity
% within the Ispec cluster + log10 alpha power within the Ispec-derived
% alpha ROI, one pass per subject (grid array job).
%
% TRIAL RULE (2026-09-26, identical across ALL analyses -- Ispec,
% Catspec, alpha_pre2post, PF, single-trial): exemplars 1-2 only,
% presentations 1-2, no 'new' responses. The pairs themselves are
% defined upstream in rsa_subject_level.m (Isim: presentation 1 vs 2 of
% the same item; Csim: exemplar duos (1,2) and (3,4) of a category,
% presentation 1, BOTH exemplars 'old'/'similar'); every table here
% takes those pairs as they are.
%
% Outputs (per subject, into pdat.out/<variant>/):
%   T_item  : one row per Isim pair. Acc: old=1, similar=0. pow1st /
%             pow2nd = alpha at the item's presentation-1 / -2 trial;
%             powAVG = their mean.
%   T_cat   : one row per Csim pair. Acc = 1 only if BOTH exemplars
%             'old', else 0. pow1st / pow2nd = alpha at exemplar 1 /
%             exemplar 2 (presentation 1 each; itemLabelB = exemplar 2);
%             powAVG = their mean.
%   T_bcat  : (only if the 'Bsim' files exist) one row per category:
%             exemplar 1 of A vs exemplar 2 of every other category,
%             averaged. pow1st = alpha at exemplar 1 of A; powPartners =
%             mean alpha over the partner trials.
%   T_joint : one row per T_item item, + the Csim of its category from
%             T_cat (matched by category; NaN if the category has no Csim
%             pair). Power = pow1st (presentation 1), Power2nd, PowerAVG.
%   alphaIspec : ONE value per subject = log10 alpha power in the ROI,
%             averaged over ALL trials that enter the Ispec analysis
%             (both members of every Isim and Csim pair). Used for the
%             subject-level alpha x Ispec correlation (rsa_plots.m).
%
% Similarity mask and alpha ROI both come from the POSITIVE significant
% Ispec cluster(s) only. Alpha band = alpha_band() (8-25 Hz since 2026-10-06: the range with a
% consistent effect in the alpha_pre2post activation-vs-zero cluster).
%
%   subj          : subject ID (numeric or string)
%   freqLo,freqHi : which Ispec variant's similarity data to read (e.g.
%                   30,100 for Ispec_broadband)
%   pdat          : struct with fields .rsa, .stat, .ispecAlphaROI, .rawEnc, .out

subj = num2str(subj);
bandLabel = sprintf('%dto%dHz', freqLo, freqHi);
variantLabel = sprintf('Ispec_%s', bandLabel);
rsaDir = fullfile(pdat.rsa, variantLabel);
statDir = fullfile(pdat.stat, variantLabel);

%% Build the similarity mask -- POSITIVE clusters only, matching the original exactly

L = load(fullfile(statDir, [variantLabel '_2nd_level_cluster_stat_vsBL.mat']), 'stat2');
stat2 = L.stat2;
p_crit = stat2.cfg.alpha;

pos_signif_clust = [];
if isfield(stat2,'posclusters') && ~isempty(stat2.posclusters) % fields are absent when no cluster formed
    pos_cluster_pvals = [stat2.posclusters(:).prob];
    pos_signif_clust = find(pos_cluster_pvals < p_crit);
end
% Negative clusters deliberately NOT included -- matches
% self_sim_singlTrl_fullSample.m, where the neg-cluster loop is commented out.

% BUGFIX (2026-09-05): must iterate over the actual cluster index VALUES
% in pos_signif_clust, not over 1:count -- see rsa_group_level_fullSample.m's
% identical fix for the full explanation (silently only worked when
% significant clusters happened to be exactly {1,...,N}).
simMask = NaN(size(stat2.stat));
for sc = pos_signif_clust
    simMask(stat2.posclusterslabelmat == sc) = 1;
end

%% Extract within-item and within-cat single-trial similarity (self_sim_singlTrl_fullSample.m logic)

% ROI (2026-10-02): restrict the similarity mask to the Ispec ROI channels
% (Ispec_derive_alphaROI.m) -- the same ROI as every other analysis.
Lr = load(fullfile(pdat.ispecAlphaROI, 'Ispec_alphaROI_config.mat'), 'roiChannels');
simMask(~ismember(stat2.label, Lr.roiChannels), :, :) = NaN;

[T_item, itemTrialLabels] = extract_sim_table(fullfile(rsaDir, [subj '_within-item_sim_ST.mat']), simMask, subj, true);
[T_cat,  catTrialLabels]  = extract_sim_table(fullfile(rsaDir, [subj '_within-cat_sim_ST.mat']), simMask, subj, false);
T_cat.itemLabelB = catTrialLabels(2:2:end);   % exemplar 2's item label (trialLabels interleave A, B per kept row)

% Between-category CONTROL rows (2026-10-05; rsa_subject_level.m 'Bsim':
% exemplar 1 of category A vs exemplar 2 of EVERY other category,
% averaged -> one row per category) -- only if computed. Extracted like an
% Isim row (one reference trial: exemplar 1 of A); Acc is not used.
fBcat = fullfile(rsaDir, [subj '_between-cat_sim_ST.mat']);
haveBcat = exist(fBcat, 'file') == 2;
if haveBcat
    T_bcat = extract_sim_table(fBcat, simMask, subj, true);
end

%% Load raw EEG + wavelet TFR, restricted to the NON-CIRCULAR alpha ROI
%  (2026-09-XX): channels + time window come from the Ispec cluster
%  (Ispec_derive_alphaROI.m) -- orthogonal to the old-vs-similar response
%  distinction, unlike the retired SME-based mask (Masks_IvC.mat), which
%  was circular (SME's Item/Category conditions ARE the old/similar
%  response later predicted downstream -- see session notes). The alpha
%  BAND (alpha_band.m) comes from the separate alpha_pre2post activation-vs-zero
%  cluster. Because this ROI is a clean rectangular box (channels x band
%  x window), not an irregular per-cell significance mask like the old
%  SME one, extraction is now a plain mean over the selected data --
%  no NaN-mask multiply-and-divide-by-n_mask needed.

L = load(fullfile(pdat.ispecAlphaROI, 'Ispec_alphaROI_config.mat'), 'roiChannels');
roiChannels = L.roiChannels; % ALL channels of the Ispec cluster (Ispec_derive_alphaROI.m; core selection dropped 2026-10-02)
timeWindow = alpha_time_window(); % 0.4-1.7 s, alpha_time_window.m (2026-10-06; was the Ispec cluster's time extent)
alphaBand = alpha_band(); % defined once in alpha_band.m (8-25 Hz); wavelet foi 2:2:30 -> bins 8,10,...,24

load(fullfile(pdat.rawEnc, [subj '_Enc_allcond.mat']), 'data');
data_preprocessed = data; clear data

cfg = [];
cfg.output = 'pow';
cfg.channel = 'all';
cfg.width = 7;
cfg.keeptrials = 'yes';
cfg.method = 'wavelet';
cfg.foi = 2:2:30;
cfg.toi = -2:0.04:3.5;
cfg.pad = ceil(max(cellfun(@numel, data_preprocessed.time)) / data_preprocessed.fsample);

TFR_data_WV = ft_freqanalysis(cfg, data_preprocessed);

cfg = [];
cfg.channel = roiChannels;
cfg.latency = timeWindow;
cfg.frequency = alphaBand;
TFR_data_WV = ft_selectdata(cfg, TFR_data_WV);
TFR_data_WV.trialinfo2 = data_preprocessed.trialinfo2;

cfg = []; cfg.parameter = 'powspctrm'; cfg.operation = 'log10';
TFR_data_WV = ft_math(cfg, TFR_data_WV);
TFR_data_WV.trialinfo2 = data_preprocessed.trialinfo2;

retDat = TFR_data_WV.trialinfo2;

% ONE ROI value per trial: log10 power averaged over the ROI channels
% (all Ispec cluster channels), frequencies and time.
trlVal = mean(reshape(TFR_data_WV.powspctrm, size(TFR_data_WV.powspctrm,1), []), 2, 'omitnan');   % trials x 1

%% Add alpha power to within-item table

% Alpha of BOTH pair members and their average (2026-10-05):
%   pow1st = 1st pair member (Isim: presentation 1; Csim: exemplar 1)
%   pow2nd = 2nd pair member (Isim: presentation 2; Csim: exemplar 2)
%   powAVG = mean of the two trials (= the trials the similarity compares)
pres_ = cell2mat(retDat(:,9)); exem_ = cell2mat(retDat(:,8));
trialPow = @(lab, pr) mean(trlVal(strcmp(retDat(:,7), lab) & pres_ == pr), 'omitnan');

pow_1st = nan(height(T_item),1); pow_2nd = nan(height(T_item),1);
for t = 1:height(T_item)
    pow_1st(t) = trialPow(T_item.itemLabel{t}, 1);
    pow_2nd(t) = trialPow(T_item.itemLabel{t}, 2);
end
T_item.pow1st = pow_1st;
T_item.pow2nd = pow_2nd;
T_item.powAVG = mean([pow_1st pow_2nd], 2, 'omitnan');

%% Add alpha power to within-cat table: exemplar 1 and exemplar 2
%  (presentation 1 each) and their average
pow_1st = nan(height(T_cat),1); pow_2nd = nan(height(T_cat),1);
for t = 1:height(T_cat)
    pow_1st(t) = trialPow(T_cat.itemLabel{t}, 1);
    pow_2nd(t) = trialPow(T_cat.itemLabelB{t}, 1);
end
T_cat.pow1st = pow_1st;
T_cat.pow2nd = pow_2nd;
T_cat.powAVG = mean([pow_1st pow_2nd], 2, 'omitnan');

%% Alpha for the between-category rows: the reference trial (exemplar 1 of
%  A, presentation 1) = pow1st, and the mean over its partner trials
%  (exemplar 2, presentation 1, 'old'/'similar', all other categories)
if haveBcat
    validR = strcmp(retDat(:,11),'old') | strcmp(retDat(:,11),'similar');
    pow_1st = nan(height(T_bcat),1); pow_partners = nan(height(T_bcat),1);
    for t = 1:height(T_bcat)
        pow_1st(t) = trialPow(T_bcat.itemLabel{t}, 1);
        pow_partners(t) = mean(trlVal(exem_ == 2 & pres_ == 1 & validR & ~strcmp(retDat(:,5), T_bcat.catLabel{t})), 'omitnan');
    end
    T_bcat.pow1st = pow_1st;
    T_bcat.powPartners = pow_partners;
end

%% Per-subject alpha: mean over ALL trials entering the Ispec analysis
%  (2026-09-26) -- Isim: presentations 1 and 2 of each paired item; Csim:
%  presentation 1 of both exemplars of each pair. Union, so a trial used
%  by both measures counts once.

pres = cell2mat(retDat(:,9));
idxIspec = find((ismember(retDat(:,7), itemTrialLabels) & (pres == 1 | pres == 2)) | ...
                (ismember(retDat(:,7), catTrialLabels) & pres == 1));
alphaIspec = mean(trlVal(idxIspec), 'omitnan'); %#ok<NASGU>
nTrialsAlpha = length(idxIspec); %#ok<NASGU>

%% Build the joint table (T_joint) for the combined Power x Isim x Csim x
%  Age model (2026-10-01, simplified): one row per item = T_item (the
%  item's own Isim, Acc, pow1st) + the Csim of its category from T_cat,
%  matched by category label. Valid because with exemplars 1-2 only each
%  category has exactly ONE Csim pair per subject (checked below). Items
%  whose category has no Csim pair get Csim = NaN. (Replaces the earlier
%  per-item re-extraction from the Csim sim_mat file.)

[catU, iFirst] = unique(T_cat.catLabel);
if numel(catU) ~= height(T_cat)
    error('rsa_singleTrial_merged:dupCat', 'Subject %s: T_cat has more than one Csim pair per category', subj);
end
[hasCsim, iCat] = ismember(T_item.catLabel, catU);
Csim = nan(height(T_item),1); Csim_noBL = nan(height(T_item),1);
Csim(hasCsim)      = T_cat.Sim(iFirst(iCat(hasCsim)));
Csim_noBL(hasCsim) = T_cat.Sim_noBL(iFirst(iCat(hasCsim)));

T_joint = table(T_item.ID, T_item.catLabel, T_item.itemLabel, T_item.condition, T_item.itemRep, T_item.catSize, ...
    T_item.Sim, T_item.Sim_noBL, Csim, Csim_noBL, T_item.Acc, T_item.pow1st, T_item.pow2nd, T_item.powAVG, ...
    'VariableNames', {'ID','catLabel','itemLabel','condition','itemRep','catSize','Isim','Isim_noBL','Csim','Csim_noBL','Acc', ...
    'Power','Power2nd','PowerAVG'});   % Power = presentation 1 (as before); Power2nd = presentation 2; PowerAVG = mean

%% Save -- per-subject, so the array job can run in parallel; a follow-up
%  step concatenates all subjects into the full-sample tables (matching
%  the ORIGINAL scripts' single combined-table output).

outDir = fullfile(pdat.out, variantLabel);
if ~exist(outDir,'dir'), mkdir(outDir); end
save(fullfile(outDir, [subj '_singleTrial_within-item.mat']), 'T_item', '-v7.3');
save(fullfile(outDir, [subj '_singleTrial_within-cat.mat']), 'T_cat', '-v7.3');
save(fullfile(outDir, [subj '_singleTrial_joint.mat']), 'T_joint', '-v7.3');
if haveBcat
    save(fullfile(outDir, [subj '_singleTrial_between-cat.mat']), 'T_bcat', '-v7.3');
end
save(fullfile(outDir, [subj '_alphaIspec.mat']), 'alphaIspec', 'nTrialsAlpha');
fprintf('Saved single-trial merged tables (within-item, within-cat, joint) for subject %s -> %s\n', subj, outDir);

end


function [T, trialLabels] = extract_sim_table(loadFile, mask, subj, isItem)
% Ported from self_sim_singlTrl_fullSample.m's per-trial extraction loop
% (both within-item and within-cat blocks use identical logic apart from
% which trialinfo2 columns get pulled) -- unchanged apart from being
% split into its own function so it can be called for both conditions.
%
% ROW INCLUSION + ACC (2026-09-26, supersedes the note below where they
% differ): rows with any 'new' response are excluded (cat-sim: BOTH
% exemplars must be 'old'/'similar'), so a cat-sim pair where one
% exemplar was forgotten no longer enters as Acc=0 (or, for Hat A, gets
% silently dropped while Hat B's 'new' was kept). Acc = 1 only if both
% are 'old', else 0.
%
% ACC CODING (2026-09-XX, cat-sim only -- see session notes): the old
% coding used Hat A's (the reference trial's) response alone for BOTH
% item-sim and cat-sim, discarding Hat B's own outcome entirely for
% cat-sim -- but cat-sim's Sim is a RELATIONAL quantity (A vs B), and Hat
% A's own recognition accuracy has no strong theoretical link to how
% similar it is to a DIFFERENT item. Likely why no Power x CatSim
% interaction on memory was ever found. Switched to CONJUNCTIVE coding
% for cat-sim: Acc=1 only if BOTH Hat A and Hat B were precisely
% recognized ('old'); Acc=0 if EITHER shows a 'similar' (imprecision)
% response or Hat B was 'new' -- directly testing whether category-level
% neural overlap costs item-specific precision for the PAIR, which is
% what dedifferentiation actually predicts and connects directly to the
% existing old2sim behavioral finding. Item-sim's own Acc (a single,
% non-relational item's accuracy) is UNCHANGED.

L = load(loadFile, 'sim_mat');
sim_mat = L.sim_mat;
trlinfo = sim_mat.trialinfo2;
pow = sim_mat.powspctrm;
powNoBL = sim_mat.powspctrmOrig;
if ~isItem
    trlinfoB = sim_mat.trialinfo2_B; % Hat B's full trialinfo2 row, per pair
end

ID = []; catLabel = {}; itemLabel = {}; condition = []; itemRep = []; catSize = [];
Sim = []; Sim_noBL = []; Acc = [];
trialLabels = {}; % item labels of every trial behind the kept rows (both pair members for cat-sim)

validResp = @(r) strcmp(r,'old') || strcmp(r,'similar');
for t = 1:size(trlinfo,1)
    % Row inclusion: NO 'new' responses in the single-trial tables
    % (2026-09-26). Item-sim: the item's own response (both presentations
    % share it). Cat-sim (2026-10-19): AT LEAST ONE exemplar 'old'/'similar'
    % -- the same rule rsa_subject_level.m applies when forming the pairs,
    % so this is a safety check (no rows should be dropped here). Every
    % remembered item of the category then carries its Csim in T_joint;
    % the cat-sim row's own Acc stays conjunctive (1 only if both 'old').
    if isItem
        keepRow = validResp(trlinfo{t,11});
    else
        keepRow = validResp(trlinfo{t,11}) || validResp(trlinfoB{t,11}); % 2026-10-19: at least one exemplar remembered (= the Csim rule in rsa_subject_level.m)
    end
    if keepRow
        ID(end+1,1) = str2double(subj); %#ok<AGROW>

        tmp = squeeze(pow(t,:,:,:)) .* mask;
        Sim(end+1,1) = mean(tmp(:), 'omitnan'); %#ok<AGROW>
        tmp = squeeze(powNoBL(t,:,:,:)) .* mask;
        Sim_noBL(end+1,1) = mean(tmp(:), 'omitnan'); %#ok<AGROW>

        itemLabel{end+1,1} = trlinfo{t,7}; %#ok<AGROW> % {} not (): plain char, not a nested 1x1 cell (T_joint's unique()/strcmp need a cellstr)
        trialLabels{end+1,1} = trlinfo{t,7}; %#ok<AGROW>
        if ~isItem, trialLabels{end+1,1} = trlinfoB{t,7}; end %#ok<AGROW>
        catLabel{end+1,1} = trlinfo{t,5}; %#ok<AGROW>

        if isItem
            Acc(end+1,1) = double(strcmp(trlinfo{t,11},'old')); %#ok<AGROW> % old=1, similar=0
        else
            Acc(end+1,1) = double(strcmp(trlinfo{t,11},'old') && strcmp(trlinfoB{t,11},'old')); %#ok<AGROW> % both old=1, else (>=1 similar)=0
        end

        cond4 = trlinfo(t,4);
        if strcmp(cond4,'normal')
            condition(end+1,1) = 1; itemRep(end+1,1) = -1; catSize(end+1,1) = -1; %#ok<AGROW>
        elseif strcmp(cond4,'high-repetition')
            condition(end+1,1) = 2; itemRep(end+1,1) = 1; catSize(end+1,1) = -1; %#ok<AGROW>
        elseif strcmp(cond4,'high-exemplar')
            condition(end+1,1) = 3; itemRep(end+1,1) = -1; catSize(end+1,1) = 1; %#ok<AGROW>
        elseif strcmp(cond4,'high-repetition+exemplar')
            condition(end+1,1) = 4; itemRep(end+1,1) = 1; catSize(end+1,1) = 1; %#ok<AGROW>
        end
    end
end

if isItem
    T = table(ID, catLabel, itemLabel, condition, itemRep, catSize, Sim, Sim_noBL, Acc);
else
    T = table(ID, catLabel, itemLabel, condition, itemRep, catSize, Sim, Sim_noBL, Acc);
end

end

function rsa_group_level_perAge(group, contrastType, freqLo, freqHi, pdat, simField)
% Consolidated per-age-group step: grand average (step2) + 1st-level
% within-person paired t-test, for ONE age group.
%
% Split out from the fullSample step (rsa_group_level_fullSample.m) so the
% dependency between them is explicit in the grid submission (fullSample
% needs BOTH age groups' output from this step) rather than implicit in
% run order, which is how the original pipeline handled it (fullSample
% logic was only triggered "if group==OA", silently assuming the YA job
% had already finished and saved its file -- fragile under parallel
% submission).
%
%   group         : 'YA' or 'OA'
%   contrastType  : 'Ispec' or 'Catspec'
%   freqLo,freqHi : analysis sub-band in Hz (must match what was used for
%                   rsa_subject_level.m for these subjects)
%   pdat          : struct with fields .rsa, .ga, .stat (paths)
%   simField      : (optional, default 'powspctrm') 'powspctrm'
%                   (baseline-corrected, the primary analysis) or
%                   'powspctrmOrig' (uncorrected) -- 2026-09-XX addition
%                   so the whole group-level stream can be run on both,
%                   to see what the baseline correction actually does.
%                   The grand-average step below already saves BOTH
%                   variants regardless (sim_GA.mat / sim_GA_orig.mat) --
%                   this parameter only affects the 1st-level test, whose
%                   output filename gets an '_uncorrected' suffix when
%                   simField='powspctrmOrig'.
%
% Names A/B follow rsa_subject_level.m: Ispec -> within-item(A)/within-cat(B);
% Catspec -> within-cat(A)/between-cat(B).

if nargin < 6 || isempty(simField), simField = 'powspctrm'; end

switch contrastType
    case {'Ispec', 'IspecAll'} % IspecAll (2026-10-18): Ispec without the response filter, own folders
        nameA = 'within-item'; nameB = 'within-cat';
    case 'Catspec'
        nameA = 'within-cat';  nameB = 'between-cat';
    otherwise
        error('rsa_group_level_perAge:contrastType', 'Unknown contrastType "%s"', contrastType)
end

bandLabel = sprintf('%dto%dHz', freqLo, freqHi);
variantLabel = sprintf('%s_%s', contrastType, bandLabel);
rsaDir = fullfile(pdat.rsa, variantLabel);
gaDir  = fullfile(pdat.ga, variantLabel);
statDir = fullfile(pdat.stat, variantLabel);
if ~exist(gaDir,'dir'), mkdir(gaDir); end
if ~exist(statDir,'dir'), mkdir(statDir); end

subjects = ids_ALL(lower(group));

%% === Step 2 equivalent: grand average (BL-corrected, BL-only, and orig versions) ===
% BUGFIX (2026-08-25, carried over from Final_RSA audit): accumulator
% cells are explicitly reset before use (see step2_rsa_group_ga.m fix in
% session notes Part 13) -- no stale-carryover risk here since each name
% (A/B) is handled in its own fresh scope below rather than a shared loop.
%
% GUARD (2026-09-XX, bugfix): this section is simField-INDEPENDENT --
% grand_average_one_condition always computes and saves all three
% variants (sim_GA/_BL/_orig) regardless of which field the 1st-level
% test below will use. Only run it on the primary (simField='powspctrm')
% call. The submission script now launches the corrected and uncorrected
% per-age jobs for the same group IN PARALLEL -- running this
% unconditionally caused both jobs to write the SAME files
% simultaneously, corrupting the HDF5 (-v7.3) write (observed: "attribute
% already exists" / "object doesn't exist" errors, NOT a stale-file
% issue -- confirmed no old files were present).

if strcmp(simField, 'powspctrm')
    [gaA, gaA_BL, gaA_Orig, includedA] = grand_average_one_condition(subjects, rsaDir, nameA);
    sim = gaA; sim.subj = includedA; %#ok<NASGU>
    save(fullfile(gaDir, [group '_' nameA '_sim_GA.mat']), 'sim', '-v7.3')
    sim = gaA_BL; sim.subj = includedA; %#ok<NASGU>
    save(fullfile(gaDir, [group '_' nameA '_sim_GA_BL.mat']), 'sim', '-v7.3')
    sim = gaA_Orig; sim.subj = includedA; %#ok<NASGU>
    save(fullfile(gaDir, [group '_' nameA '_sim_GA_orig.mat']), 'sim', '-v7.3')

    [gaB, gaB_BL, gaB_Orig, includedB] = grand_average_one_condition(subjects, rsaDir, nameB);
    sim = gaB; sim.subj = includedB; %#ok<NASGU>
    save(fullfile(gaDir, [group '_' nameB '_sim_GA.mat']), 'sim', '-v7.3')
    sim = gaB_BL; sim.subj = includedB; %#ok<NASGU>
    save(fullfile(gaDir, [group '_' nameB '_sim_GA_BL.mat']), 'sim', '-v7.3')
    sim = gaB_Orig; sim.subj = includedB; %#ok<NASGU>
    save(fullfile(gaDir, [group '_' nameB '_sim_GA_orig.mat']), 'sim', '-v7.3')

    fprintf('Saved grand averages for %s, %s, %s\n', group, variantLabel, gaDir);
else
    fprintf('Skipping grand-average computation for %s, %s (simField=%s -- GA files are simField-independent and already written by the primary/corrected run)\n', group, variantLabel, simField);
end

%% === 1st-level within-person test (A vs B): indepsamplesT, all raw pairs ===
% REVERTED (2026-09-XX): back to the TRUE original method the manuscript
% is based on -- all raw item-sim pairs vs all raw cat-sim pairs,
% unpaired (indepsamplesT), faithfully restoring the original script's
% indepsamplesT branch (see independent_samples_ttest's header). Tried
% category_level_paired_ttest (still kept below, deprecated) as a more
% theoretically matched alternative; it came out empirically weaker
% (smaller, borderline cluster) for reasons diagnosed during the
% investigation (likely (a) the r=.86 correlation motivating the
% pairing's expected power gain was measured BETWEEN subjects, not
% WITHIN a subject across categories, so pairing removed little actual
% shared noise; and (b) averaging Hat A's and Hat B's own item-sim
% together can dilute the isim-vs-csim contrast itself when the two
% exemplars differ genuinely in encoding precision, not just noise). Both
% the paired and unpaired designs are legitimate ways to construct a
% per-subject summary statistic for the 2nd-level cluster test -- this is
% not a validity fix, it's reverting to the already-published, already-
% validated approach given the alternative didn't pay off empirically.

stat = independent_samples_ttest(subjects, rsaDir, nameA, nameB, simField);
statSuffix = '';
if ~strcmp(simField, 'powspctrm'), statSuffix = '_uncorrected'; end
save(fullfile(statDir, [group '_' nameA '_vs_' nameB '_1st_level_stat_vsBL' statSuffix '.mat']), 'stat', '-v7.3')
fprintf('Saved 1st-level stat (indepsamplesT, all raw pairs, simField=%s) for %s, %s -> %s\n', simField, group, variantLabel, statDir);

end


%% =========================================================
%% LOCAL FUNCTIONS
%% =========================================================

function [avg, avgBL, avgOrig, included] = grand_average_one_condition(subjects, rsaDir, condName)
% Ported from step2_rsa_group_ga.m -- averages trial-pairs within each
% subject, then stacks subjects (keepindividual) into one GA structure.

avg_cell = {}; avgBL_cell = {}; avgOrig_cell = {}; included = {};
for s = 1:length(subjects)
    curr_subj = num2str(subjects(s));
    loadFile = fullfile(rsaDir, [curr_subj '_' condName '_sim_ST.mat']);
    if ~exist(loadFile,'file')
        warning('Missing similarity file for subject %s (%s) -- skipping', curr_subj, condName);
        continue
    end
    L = load(loadFile,'sim_mat');
    sim_mat = L.sim_mat;
    sim_mat.powspctrm = single(sim_mat.powspctrm);

    simOrig = sim_mat; simOrig.powspctrm = single(simOrig.powspctrmOrig);
    simBL = sim_mat; simBL.powspctrm = single(simBL.powspctrmBL);

    ft_warning off 'FieldTrip:getdimord:warning_dimord_could_not_be_determined'
    if ~isempty(sim_mat.powspctrm)
        cfg = []; cfg.avgoverrpt = 'yes'; cfg.parameter = 'powspctrm';
        avg_cell{end+1} = ft_freqdescriptives(cfg, sim_mat); %#ok<AGROW>
        avgBL_cell{end+1} = ft_freqdescriptives(cfg, simBL); %#ok<AGROW>
        avgOrig_cell{end+1} = ft_freqdescriptives(cfg, simOrig); %#ok<AGROW>
        included{end+1} = curr_subj; %#ok<AGROW>
    end
end

cfg = []; cfg.keepindividual = 'yes';
avg = ft_freqgrandaverage(cfg, avg_cell{:});
avgBL = ft_freqgrandaverage(cfg, avgBL_cell{:});
avgOrig = ft_freqgrandaverage(cfg, avgOrig_cell{:});

end


function stat_I = first_level_paired_ttest(subjects, rsaDir, nameA, nameB) %#ok<DEFNU>
% DEPRECATED (2026-09-XX) -- NO LONGER CALLED. Kept only as an audit
% trail; do not revive without re-reading this header.
%
% Ported from step4a_sim_comparison_1st_level_vsBL.m. Investigating this
% function while designing the new category_level_paired_ttest (below)
% uncovered an important discrepancy: the ORIGINAL script
% (Final_RSA/dcn-tutorial-master/step4a_sim_comparison_1st_level_vsBL.m)
% has BOTH 'depsamplesT' (this function, pairing via "intersect on
% trialinfo2 column 1" -- exact meaning of column 1 never fully
% re-derived) and 'indepsamplesT' (all raw isim pairs vs all raw csim
% pairs, unpaired) available via a simple comment-toggle, with
% 'indepsamplesT' commented out in the copy this was ported from.
% User confirmed with certainty that the PUBLISHED manuscript numbers
% used 'indepsamplesT', not this depsamplesT/column-1-intersect version
% -- meaning that toggle must have been switched (likely during the
% user's own past exploration of the dependent-samples alternative,
% which they tried and rejected at the time specifically because it
% silently discarded Hat B's own item-sim entirely) sometime after the
% manuscript was generated, and that later state is what got carried
% into this RSA_pipeline port. Neither this function NOR the published
% indepsamplesT approach is what's used going forward -- both are
% superseded by category_level_paired_ttest, a third, deliberately-
% designed option (average Hat A's + Hat B's own item-sim per category,
% contrast against cat-sim, symmetric response filter) that fixes the
% Hat-B-exclusion problem this function had without reintroducing the
% published version's own issues (see session notes for the full
% discussion).

stats_cell = {}; dfs = [];
for s = 1:length(subjects)
    curr_subj = num2str(subjects(s));
    dat1 = load(fullfile(rsaDir, [curr_subj '_' nameA '_sim_ST.mat']));
    dat2 = load(fullfile(rsaDir, [curr_subj '_' nameB '_sim_ST.mat']));

    [~,ia,ib] = intersect(cell2mat(dat1.sim_mat.trialinfo2(:,1)), cell2mat(dat2.sim_mat.trialinfo2(:,1)));

    ft_warning off 'FieldTrip:getdimord:warning_dimord_could_not_be_determined'
    cfg = []; cfg.trials = ia;
    dat1 = ft_selectdata(cfg, dat1.sim_mat);
    cfg.trials = ib;
    dat2 = ft_selectdata(cfg, dat2.sim_mat);

    n = size(dat1.powspctrm,1);

    cfg = [];
    cfg.channel = 'all'; cfg.latency = 'all'; cfg.frequency = 'all';
    cfg.avgovertime = 'no'; cfg.avgoverfreq = 'no'; cfg.avgoverchan = 'no';
    cfg.parameter = 'powspctrm';
    cfg.method = 'analytic'; cfg.tail = 0; cfg.statistic = 'depsamplesT';
    cfg.design(1,1:n+n) = [ones(1,n) 2*ones(1,n)];
    cfg.design(2,1:n+n) = [1:n 1:n];
    cfg.ivar = 1; cfg.uvar = 2;

    stats_cell{s} = ft_freqstatistics(cfg, dat1, dat2); %#ok<AGROW>
    dfs(s) = stats_cell{s}.df; %#ok<AGROW>
end

cfg = []; cfg.keepindividual = 'yes'; cfg.parameter = 'stat';
stat_I = ft_freqgrandaverage(cfg, stats_cell{:});
stat_I.df = dfs;

end


function stat_I = category_level_paired_ttest(subjects, rsaDir, nameA, nameB, simField) %#ok<DEFNU>
% DEPRECATED (2026-09-XX) -- NO LONGER CALLED, reverted back to
% independent_samples_ttest (below), matching the manuscript's original
% method. Kept as an audit trail; do not revive without re-reading this.
%
% Tried as a more theoretically matched alternative to the published
% indepsamplesT approach (per-category pairing instead of comparing all
% raw isim pairs against all raw csim pairs, unpaired) -- it came out
% empirically weaker (smaller, borderline-significant 2nd-level cluster)
% once actually run. Likely contributors, diagnosed but not fully
% disentangled: (a) the r=.86 within-item/within-category correlation
% that motivated expecting a power gain from pairing was measured BETWEEN
% subjects (each subject's own mean), not WITHIN a subject across
% categories -- the level this design's pairing actually operates at --
% so little of that correlation's noise-cancelling benefit may have
% carried over; (b) averaging Hat A's and Hat B's own item-sim together
% can dilute the isim-vs-csim contrast itself, not just add noise, when
% the two exemplars genuinely differ in encoding precision (e.g. one
% judged 'old', the other 'similar') rather than merely being two noisy
% estimates of the same true value. Both this design and indepsamplesT
% are legitimate ways to build a per-subject summary statistic for the
% 2nd-level cluster test -- reverting was a power/consistency decision,
% not a validity fix.
%
% Original header, for reference:
% Replaces BOTH the published manuscript's method (indepsamplesT, all raw
% item-sim pairs vs all raw cat-sim pairs, unpaired) and the dormant
% depsamplesT/column-1-intersect option in the original script (never
% used for publication; see first_level_paired_ttest's header for the
% full discovery). See session notes for the complete design discussion.
%
% Design: per category, average Hat A's own item-sim (presentation
% 1-vs-2) and Hat B's own item-sim (presentation 1-vs-2) into one
% per-category item-sim value; contrast against that category's cat-sim
% value (Hat A vs Hat B, exemplar 1 vs 2). Item-sim/cat-sim pairing
% mechanics themselves are UNCHANGED from rsa_subject_level.m (still
% strictly presentation-1-vs-2 / exemplar-1-vs-2 -- deliberately NOT
% expanded to all repetitions/exemplars, given real age-differential
% repetition-suppression evidence in this dataset, ERP and alpha desync).
%
% Requires (tightening the original asymmetric within-cat rule):
%   - Hat B ALSO has a valid ('old'/'similar') response (not just Hat A)
%   - BOTH Hat A and Hat B individually have their own presentation-2
%     (so each side's own item-sim value is actually computable)
% Checked empirically (diagnostics/check_isim_csim_trial_counts.m): this
% costs ~38% of cat-sim's pairs, symmetric across age groups (no
% differential N bias). No subject exclusions applied -- full 93-subject
% sample used throughout (low-N subjects can only add noise/mask a true
% effect in this group-level permutation framework, not manufacture a
% false one).
%
%   simField : 'powspctrm' (baseline-corrected, primary) or
%              'powspctrmOrig' (uncorrected, for the corrected-vs-
%              uncorrected comparison) -- selects which field of the
%              per-subject sim_mat structs feeds this test.

stats_cell = {}; dfs = [];
for s = 1:length(subjects)
    curr_subj = num2str(subjects(s));
    L = load(fullfile(rsaDir, [curr_subj '_' nameA '_sim_ST.mat']), 'sim_mat'); datI = L.sim_mat;
    L = load(fullfile(rsaDir, [curr_subj '_' nameB '_sim_ST.mat']), 'sim_mat'); datC = L.sim_mat;

    itemLabelsI = datI.trialinfo2(:,7); % item-sim's own reference item, per row

    nCat = size(datC.trialinfo2,1);
    keepCat = false(nCat,1);
    idxItemA = nan(nCat,1);
    idxItemB = nan(nCat,1);
    for p = 1:nCat
        respB = datC.trialinfo2_B{p,11};
        if ~(strcmp(respB,'old') || strcmp(respB,'similar')), continue; end % symmetric response filter

        itemA = datC.trialinfo2{p,7};
        itemB = datC.trialinfo2_B{p,7};
        iA = find(strcmp(itemLabelsI, itemA), 1);
        iB = find(strcmp(itemLabelsI, itemB), 1);
        if isempty(iA) || isempty(iB), continue; end % both need their own presentation-2 (own item-sim value)

        keepCat(p) = true;
        idxItemA(p) = iA;
        idxItemB(p) = iB;
    end

    pIdx = find(keepCat);
    n = length(pIdx);
    if n < 2
        warning('Subject %s has fewer than 2 qualifying categories for the paired design -- skipping', curr_subj);
        continue
    end

    ft_warning off 'FieldTrip:getdimord:warning_dimord_could_not_be_determined'
    cfg = []; cfg.trials = pIdx;
    dat2 = ft_selectdata(cfg, datC); % cat-sim, restricted to qualifying categories

    avgItemSim = nan(size(dat2.(simField)));
    for k = 1:n
        avgItemSim(k,:,:,:) = (datI.(simField)(idxItemA(pIdx(k)),:,:,:) + datI.(simField)(idxItemB(pIdx(k)),:,:,:)) / 2;
    end
    dat1 = dat2; % same trial dimension/metadata as the qualifying cat-sim subset
    dat1.(simField) = avgItemSim;

    cfg = [];
    cfg.channel = 'all'; cfg.latency = 'all'; cfg.frequency = 'all';
    cfg.avgovertime = 'no'; cfg.avgoverfreq = 'no'; cfg.avgoverchan = 'no';
    cfg.parameter = simField;
    cfg.method = 'analytic'; cfg.tail = 0; cfg.statistic = 'depsamplesT';
    cfg.design(1,1:n+n) = [ones(1,n) 2*ones(1,n)];
    cfg.design(2,1:n+n) = [1:n 1:n];
    cfg.ivar = 1; cfg.uvar = 2;

    stats_cell{end+1} = ft_freqstatistics(cfg, dat1, dat2); %#ok<AGROW>
    dfs(end+1) = stats_cell{end}.df; %#ok<AGROW>
end

cfg = []; cfg.keepindividual = 'yes'; cfg.parameter = 'stat';
stat_I = ft_freqgrandaverage(cfg, stats_cell{:});
stat_I.df = dfs;

end


function stat_I = independent_samples_ttest(subjects, rsaDir, nameA, nameB, simField)
% RESTORED (2026-09-XX) -- the TRUE method the manuscript is actually
% based on: all raw item-sim pairs vs all raw cat-sim pairs, UNPAIRED
% (indepsamplesT), at each (chan,time,time) cell, per subject. Faithfully
% restores the original script's (step4a_sim_comparison_1st_level_vsBL.m)
% 'indepsamplesT' branch, which was commented out in the copy this
% pipeline was ported from in favor of the dormant depsamplesT/column-1-
% intersect option (first_level_paired_ttest, never used for publication
% -- see that function's header for the discovery). No trial-matching,
% no averaging, no symmetric response filter -- cat-sim keeps its
% original, asymmetric Hat-A-only response rule unchanged. Generic across
% contrastType (works for Ispec's within-item/within-cat or Catspec's
% within-cat/between-cat identically, matching the original script).
%
%   simField : 'powspctrm' (baseline-corrected, primary) or
%              'powspctrmOrig' (uncorrected, for the corrected-vs-
%              uncorrected comparison).

stats_cell = {}; dfs = [];
for s = 1:length(subjects)
    curr_subj = num2str(subjects(s));
    L = load(fullfile(rsaDir, [curr_subj '_' nameA '_sim_ST.mat']), 'sim_mat'); dat1 = L.sim_mat;
    L = load(fullfile(rsaDir, [curr_subj '_' nameB '_sim_ST.mat']), 'sim_mat'); dat2 = L.sim_mat;

    n1 = size(dat1.(simField),1);
    n2 = size(dat2.(simField),1);

    cfg = [];
    cfg.channel = 'all'; cfg.latency = 'all'; cfg.frequency = 'all';
    cfg.avgovertime = 'no'; cfg.avgoverfreq = 'no'; cfg.avgoverchan = 'no';
    cfg.parameter = simField;
    cfg.method = 'analytic'; cfg.tail = 0; cfg.statistic = 'indepsamplesT';
    cfg.design = zeros(1, n1+n2);
    cfg.design(1,1:n1) = 1;
    cfg.design(1,n1+1:n1+n2) = 2;
    cfg.ivar = 1;

    ft_warning off 'FieldTrip:getdimord:warning_dimord_could_not_be_determined'
    stats_cell{end+1} = ft_freqstatistics(cfg, dat1, dat2); %#ok<AGROW>
    dfs(end+1) = stats_cell{end}.df; %#ok<AGROW>
end

cfg = []; cfg.keepindividual = 'yes'; cfg.parameter = 'stat';
stat_I = ft_freqgrandaverage(cfg, stats_cell{:});
stat_I.df = dfs;

end

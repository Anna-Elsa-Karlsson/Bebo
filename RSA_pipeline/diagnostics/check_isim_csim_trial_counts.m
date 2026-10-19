function check_isim_csim_trial_counts()
% DIAGNOSTIC ONLY -- does not touch the production pipeline or its output
% files. Reports per-subject trial counts under the CURRENT (2026-09-XX,
% FINAL) cat-sim pairing rule, ported exactly from
% pairs_within_cat_exemplar1v2 in rsa_subject_level.m:
%
%   - Two exemplar duos per category: (1,2) and (3,4). Low-catSize
%     categories only ever have duo (1,2); high-catSize categories can
%     contribute both.
%   - OR-validity per duo: EITHER member of the duo having a valid
%     ('old'/'similar') response is enough to include that duo (not a
%     fixed anchor, not both required).
%
% Reports duo (1,2) and duo (3,4) counts separately, plus their sum
% (n_csim_current -- what the production pipeline actually uses), so you
% can see how much the exemplar-3/4 extension adds on top of the
% original single-duo design.
%
% n_isim_current -- item-sim pairs, unchanged rule: presentation-1 needs
% a valid 'old'/'similar' response; presentation-2 found by item label,
% no separate response check (shares presentation-1's response by
% construction, same item, one eventual test outcome). No exemplar
% restriction at all (includes exemplar 3/4 items too, unlike Csim).
%
% Pairing logic replicated EXACTLY from rsa_subject_level.m's
% pairs_within_item_vs_cat / pairs_within_cat_exemplar1v2 (trialinfo2
% column conventions: 4=condition, 5=category, 7=item label, 8=exemplar,
% 9=presentation, 11=response) -- just counting, no spectral computation,
% so no FieldTrip/TFR/BOSC dependency and no grid array job needed; this
% runs over all subjects directly in one session.
%
% Usage: matlab -r "check_isim_csim_trial_counts; exit"
% Run from anywhere -- paths are relative to this file's own location.

scriptDir = fileparts(mfilename('fullpath'));
pipelineRoot = fullfile(scriptDir, '..');
addpath(pipelineRoot); % for ids_ALL.m

rawEncDir = fullfile(pipelineRoot, '..', 'PreprocData');

ids = ids_ALL('all');
nAll = length(ids);

n_isim_current    = nan(nAll,1);
n_csim_duo12_or   = nan(nAll,1);
n_csim_duo34_or   = nan(nAll,1);
n_csim_current    = nan(nAll,1); % duo12 + duo34 -- what the production pipeline uses
AgeGroup = cell(nAll,1);

for s = 1:nAll
    subj = num2str(ids(s));
    AgeGroup{s} = ternary(ids(s) < 4000, 'YA', 'OA');

    loadFile = fullfile(rawEncDir, subj, 'eeg', [subj '_Enc_allcond.mat']);
    if ~exist(loadFile,'file')
        warning('Missing Enc_allcond for subject %s -- skipping', subj);
        continue
    end
    L = load(loadFile, 'data');
    trlinfo = L.data.trialinfo2;
    clear L

    %% Item-sim, unchanged rule (ported from pairs_within_item_vs_cat)
    pairsI = [];
    for stim1 = 1:size(trlinfo,1)
        if ~ismember(stim1,pairsI) && trlinfo{stim1,9} == 1 && ~isempty(find(strcmp(trlinfo{stim1,11},'old') | strcmp(trlinfo{stim1,11},'similar')))
            curr_pair = trlinfo{stim1,7};
            if strcmp(trlinfo(stim1,4),'normal') || strcmp(trlinfo(stim1,4),'high-exemplar')
                if find(strcmp(trlinfo(:,7),curr_pair) & cell2mat(trlinfo(:,9)) == 2)
                    pairsI(stim1,2) = find(strcmp(trlinfo(:,7),curr_pair) & cell2mat(trlinfo(:,9)) == 2); %#ok<*FNDSB>
                    pairsI(stim1,1) = stim1;
                end
            elseif strcmp(trlinfo(stim1,4),'high-repetition+exemplar') || strcmp(trlinfo(stim1,4),'high-repetition')
                if find(strcmp(trlinfo(:,7),curr_pair) & cell2mat(trlinfo(:,9)) == 2)
                    pairsI(stim1,1) = stim1;
                    pairsI(stim1,2) = find(strcmp(trlinfo(:,7),curr_pair) & cell2mat(trlinfo(:,9)) == 2);
                end
            end
        end
    end
    pairsI(~all(pairsI,2),:) = [];
    n_isim_current(s) = size(pairsI,1);

    %% Cat-sim, CURRENT rule: exemplar duos (1,2) and (3,4), OR-validity
    %  each -- ported exactly from pairs_within_cat_exemplar1v2
    categories = unique(trlinfo(:,5));
    exemplarDuos = [1 2; 3 4];
    pairsByDuo = {[], []};
    for c = 1:length(categories)
        curr_cat = categories{c};
        for d = 1:size(exemplarDuos,1)
            exA = exemplarDuos(d,1); exB = exemplarDuos(d,2);
            idxA = find(strcmp(trlinfo(:,5),curr_cat) & cell2mat(trlinfo(:,8)) == exA & cell2mat(trlinfo(:,9)) == 1);
            idxB = find(strcmp(trlinfo(:,5),curr_cat) & cell2mat(trlinfo(:,8)) == exB & cell2mat(trlinfo(:,9)) == 1);
            if isempty(idxA) || isempty(idxB)
                continue
            end
            stimA = idxA(1); stimB = idxB(1);
            validA = strcmp(trlinfo{stimA,11},'old') || strcmp(trlinfo{stimA,11},'similar');
            validB = strcmp(trlinfo{stimB,11},'old') || strcmp(trlinfo{stimB,11},'similar');
            if validA || validB
                pairsByDuo{d}(end+1,:) = [stimA, stimB]; %#ok<AGROW>
            end
        end
    end
    n_csim_duo12_or(s) = size(pairsByDuo{1},1);
    n_csim_duo34_or(s) = size(pairsByDuo{2},1);
    n_csim_current(s)  = n_csim_duo12_or(s) + n_csim_duo34_or(s);

    fprintf('Subject %s (%d/%d, %s): isim=%d, csim_duo12=%d, csim_duo34=%d, csim_total=%d\n', ...
        subj, s, nAll, AgeGroup{s}, n_isim_current(s), n_csim_duo12_or(s), n_csim_duo34_or(s), n_csim_current(s));
end

%% Summary table + save

T = table(ids(:), AgeGroup, n_isim_current, n_csim_duo12_or, n_csim_duo34_or, n_csim_current, ...
    'VariableNames', {'ID','AgeGroup','n_isim_current','n_csim_duo12_or','n_csim_duo34_or','n_csim_current'});

fprintf('\n=== SUMMARY (mean +- SD [min, max]) ===\n');
report_summary(T, 'ALL');
report_summary(T(strcmp(T.AgeGroup,'YA'),:), 'YA');
report_summary(T(strcmp(T.AgeGroup,'OA'),:), 'OA');

fprintf('\n=== SUBJECTS WITH ZERO valid within-cat pairs under the current rule ===\n');
zeroSubj = T.ID(T.n_csim_current == 0);
if isempty(zeroSubj)
    fprintf('None.\n');
else
    disp(zeroSubj')
end

outDir = fullfile(pipelineRoot, '..', 'Results', 'diagnostics');
if ~exist(outDir,'dir'), mkdir(outDir); end
writetable(T, fullfile(outDir, 'isim_csim_trial_counts.csv'));
save(fullfile(outDir, 'isim_csim_trial_counts.mat'), 'T');
fprintf('\nSaved -> %s\n', fullfile(outDir, 'isim_csim_trial_counts.csv'));

end


function report_summary(T, label)
vars = {'n_isim_current','n_csim_duo12_or','n_csim_duo34_or','n_csim_current'};
fprintf('\n-- %s (n=%d subjects) --\n', label, height(T));
for v = 1:length(vars)
    x = T.(vars{v});
    fprintf('%-20s %6.1f +- %5.1f  [%3d, %3d]\n', vars{v}, mean(x,'omitnan'), std(x,'omitnan'), min(x), max(x));
end
end


function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

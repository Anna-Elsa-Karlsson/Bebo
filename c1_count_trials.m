function c1_count_trials()
% C1_COUNT_TRIALS  Trial counts per person from the preprocessed data
% (sample check, analysis plan section 0). No EEG computation: reads only
% the trial information of PreprocData/<ID>_Enc_allcond.mat.
%
% Persons: all 96 recorded (bebo_ids('full'); incl. 3111, 3213, 4210, which were excluded
% in the original 93 and are included again since 2026-10-17; column
% 'excluded' marks them), to check the manuscript's criterion "fewer than
% 15 trials per condition of interest" (conditions = later 'old' vs later
% 'similar').
%
% Per person:
%   nTrialsData          all encoding trials left after preprocessing
%   Counts are per ITEM (presentation-1 trial), as the SME (decided
%   2026-10-17):
%   NEW rule (all exemplars, 'old'/'similar'; bebo_trials.m):
%     new_old, new_sim   items later called 'old' / 'similar'
%     new_items, new_Isim, new_Csim   items, Isim pairs, Csim pairs
%   OLD rule (exemplars 1-2, 'old'/'similar'):
%     old_old, old_sim   items later called 'old' / 'similar'
%   minNew, minOld       smaller of the two conditions (criterion: >= 15)
%
% Output: C.dir.trialLog/trial_counts_all96.csv and trial_counts_summary.txt

C = bebo_config();
if ~exist(C.dir.trialLog, 'dir'), mkdir(C.dir.trialLog); end
extra = [3111 3213 4210];
ids = bebo_ids('full');
n = numel(ids);
R = table(ids(:), repmat("", n, 1), ismember(ids(:), extra), nan(n,1), nan(n,1), nan(n,1), nan(n,1), ...
    nan(n,1), nan(n,1), nan(n,1), nan(n,1), repmat("", n, 1), ...
    'VariableNames', {'ID','group','excluded','nTrialsData','new_old','new_sim','new_items', ...
    'new_Isim','new_Csim','old_old','old_sim','note'});

for s = 1:n
    id = ids(s);
    if floor(id/1000) == 3, R.group(s) = "YA"; else, R.group(s) = "OA"; end
    f = fullfile(C.preproc, sprintf('%d_Enc_allcond.mat', id));
    if ~exist(f, 'file'), R.note(s) = "no PreprocData file"; fprintf('%d: no file\n', id); continue, end
    L = load(f, 'data');
    ti = L.data.trialinfo2; clear L
    R.nTrialsData(s) = size(ti, 1);
    [T, pI, pC] = bebo_trials(ti, C);                    % new rule
    p1 = T.pres == 1;                                     % one trial per item
    R.new_old(s) = sum(T.Acc == 1 & p1); R.new_sim(s) = sum(T.Acc == 0 & p1);
    R.new_items(s) = sum(T.pres == 1); R.new_Isim(s) = size(pI, 1); R.new_Csim(s) = size(pC, 1);
    old12 = T.exemplar <= 2;                              % old rule: exemplars 1-2
    R.old_old(s) = sum(T.Acc == 1 & old12 & p1); R.old_sim(s) = sum(T.Acc == 0 & old12 & p1);
    fprintf('%d: %d trials; items new rule old/similar %d/%d; old rule %d/%d\n', id, R.nTrialsData(s), ...
        R.new_old(s), R.new_sim(s), R.old_old(s), R.old_sim(s));
end
R.minNew = min(R.new_old, R.new_sim);
R.minOld = min(R.old_old, R.old_sim);
writetable(R, fullfile(C.dir.trialLog, 'trial_counts_all96.csv'));

fid = fopen(fullfile(C.dir.trialLog, 'trial_counts_summary.txt'), 'w');
fprintf(fid, 'Trial counts from PreprocData -- %s\nCriterion (manuscript): >= 15 trials per condition (later old / later similar); counted per item (presentation 1)\n\n', datestr(now));
for g = ["YA", "OA"]
    for ex = [false true]
        ix = R.group == g & R.excluded == ex & ~isnan(R.nTrialsData);
        if ~any(ix), continue, end
        fprintf(fid, '%s, %s (n=%d): trials in data M=%.0f [%d-%d]; min condition, new rule M=%.0f [%d-%d]; old rule M=%.0f [%d-%d]\n', ...
            g, ternary(ex, 'previously excluded', 'original 93'), sum(ix), mean(R.nTrialsData(ix)), min(R.nTrialsData(ix)), max(R.nTrialsData(ix)), ...
            mean(R.minNew(ix)), min(R.minNew(ix)), max(R.minNew(ix)), mean(R.minOld(ix)), min(R.minOld(ix)), max(R.minOld(ix)));
    end
end
fprintf(fid, '\nPersons below 15 trials in a condition (old rule or new rule):\n');
low = find(R.minOld < 15 | R.minNew < 15);
if isempty(low), fprintf(fid, '  none\n'); end
for k = low(:)'
    fprintf(fid, '  %d (%s%s): old rule %d/%d, new rule %d/%d\n', R.ID(k), R.group(k), ...
        ternary(R.excluded(k), ', previously excluded', ''), R.old_old(k), R.old_sim(k), R.new_old(k), R.new_sim(k));
end
fprintf(fid, '\nThe three previously excluded IDs:\n');
for k = find(R.excluded)'
    fprintf(fid, '  %d (%s): %s trials in data; old rule %d/%d; new rule %d/%d %s\n', R.ID(k), R.group(k), ...
        num2str(R.nTrialsData(k)), R.old_old(k), R.old_sim(k), R.new_old(k), R.new_sim(k), R.note(k));
end
fclose(fid);
fprintf('%s', fileread(fullfile(C.dir.trialLog, 'trial_counts_summary.txt')));
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

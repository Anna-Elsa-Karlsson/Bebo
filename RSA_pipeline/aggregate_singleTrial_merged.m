function aggregate_singleTrial_merged(freqLo, freqHi)
% Concatenates all subjects' (ids_ALL('all')) per-subject tables (from
% rsa_singleTrial_merged.m) into the two final combined tables, matching
% the ORIGINAL self_sim_singlTrl_fullSample.m's single-combined-table
% output structure (one row per trial per subject, stacked).
%
% Run once after the rsa_singleTrial_merged.m array job completes for all
% subjects.
%
% Usage: matlab -r "aggregate_singleTrial_merged(30,100); exit"  (or via
% the config_and_run wrapper's convention if preferred -- this one takes
% plain numeric args since it's meant to be called directly, not per-task)

pipelineRoot = fileparts(mfilename('fullpath'));
addpath(pipelineRoot); % for ids_ALL.m

bandLabel = sprintf('%dto%dHz', freqLo, freqHi);
variantLabel = sprintf('Ispec_%s', bandLabel);
inDir = fullfile(pipelineRoot, '..', 'Results', 'SingleTrialMerged', variantLabel);

subjects = ids_ALL('all');

alphaIspec = nan(length(subjects),1);   % per-subject alpha (log10) over all Ispec trials, ids_ALL order
nTrialsAlpha = nan(length(subjects),1);
T_item_all = {};
T_cat_all = {};
T_joint_all = {};
T_bcat_all = {};   % between-category control pairs (only if computed, 2026-10-05)
for s = 1:length(subjects)
    curr_subj = num2str(subjects(s));

    fItem = fullfile(inDir, [curr_subj '_singleTrial_within-item.mat']);
    if exist(fItem,'file')
        L = load(fItem, 'T_item');
        T_item_all{end+1} = L.T_item; %#ok<AGROW>
    else
        warning('Missing within-item file for subject %s -- skipping', curr_subj);
    end

    fCat = fullfile(inDir, [curr_subj '_singleTrial_within-cat.mat']);
    if exist(fCat,'file')
        L = load(fCat, 'T_cat');
        T_cat_all{end+1} = L.T_cat; %#ok<AGROW>
    else
        warning('Missing within-cat file for subject %s -- skipping', curr_subj);
    end

    fAlpha = fullfile(inDir, [curr_subj '_alphaIspec.mat']);
    if exist(fAlpha,'file')
        L = load(fAlpha, 'alphaIspec', 'nTrialsAlpha');
        alphaIspec(s) = L.alphaIspec;
        nTrialsAlpha(s) = L.nTrialsAlpha;
    else
        warning('Missing alphaIspec file for subject %s -- NaN', curr_subj);
    end

    fBcat = fullfile(inDir, [curr_subj '_singleTrial_between-cat.mat']);
    if exist(fBcat,'file')
        L = load(fBcat, 'T_bcat');
        T_bcat_all{end+1} = L.T_bcat; %#ok<AGROW>
    end

    fJoint = fullfile(inDir, [curr_subj '_singleTrial_joint.mat']);
    if exist(fJoint,'file')
        L = load(fJoint, 'T_joint');
        T_joint_all{end+1} = L.T_joint; %#ok<AGROW>
    else
        warning('Missing joint file for subject %s -- skipping', curr_subj);
    end
end

T_item = vertcat(T_item_all{:});
T_cat = vertcat(T_cat_all{:});
T_joint = vertcat(T_joint_all{:});
T_bcat = vertcat(T_bcat_all{:});
if ~isempty(T_bcat) && numel(T_bcat_all) < length(subjects)
    warning('Between-category tables found for only %d of %d subjects', numel(T_bcat_all), length(subjects));
end

% Age column, added here (not per-subject) since it's a trivial function
% of ID -- matches the ORIGINAL self_sim_singlTrl_fullSample.m convention
% exactly (curr_subj < 4000 -> YA=1, else OA=0; BEBO_behavior_glmer.R
% recodes 0->-1 itself downstream, so this must stay 1/0, not 1/-1).
% Appended as the last column rather than reordered after ID -- movevars
% is R2018a+, not available on this grid's R2017b, and column order
% doesn't matter for anything that reads by name (e.g. R).
T_item.Age = double(T_item.ID < 4000);
T_cat.Age = double(T_cat.ID < 4000);
T_joint.Age = double(T_joint.ID < 4000);
if ~isempty(T_bcat), T_bcat.Age = double(T_bcat.ID < 4000); end

save(fullfile(inDir, [variantLabel '_singleTrial_within-item_fullSample.mat']), 'T_item', '-v7.3');
writetable(T_item, fullfile(inDir, [variantLabel '_singleTrial_within-item_fullSample.txt']));
save(fullfile(inDir, [variantLabel '_singleTrial_within-cat_fullSample.mat']), 'T_cat', '-v7.3');
writetable(T_cat, fullfile(inDir, [variantLabel '_singleTrial_within-cat_fullSample.txt']));
save(fullfile(inDir, [variantLabel '_singleTrial_joint_fullSample.mat']), 'T_joint', '-v7.3');
if ~isempty(T_bcat)
    save(fullfile(inDir, [variantLabel '_singleTrial_between-cat_fullSample.mat']), 'T_bcat', '-v7.3');
end
writetable(T_joint, fullfile(inDir, [variantLabel '_singleTrial_joint_fullSample.txt']));

%% Pair-level long table (2026-10-04): one row per SIMILARITY PAIR, Isim
%  pairs (T_item) and Csim pairs (T_cat) stacked, for the model
%  Sim ~ Power * Type (does alpha desync go with MORE item-specific and
%  LESS category-like representations?). No structural NaNs: every valid
%  pair is its own row. Power = mean log10 alpha of the pair's TWO trials
%  (powAVG: Isim = presentation 1+2 of the item; Csim = exemplar 1+2 of
%  the category, presentation 1), z-scored WITHIN each person across all
%  of that person's pairs (both types on one scale). Type: Isim / Csim
%  (TypeNum: Isim = 0.5, Csim = -0.5). Lower Power = stronger desync.
%  2026-10-05: + Type 'Bsim' (between-category control rows: exemplar 1 of
%  A vs exemplar 2 of every other category, averaged; TypeNum NaN) when
%  computed. Two alpha measures:
%    Power / PowerZ       mean of the pair's two trials (Isim, Csim only;
%                         NaN for Bsim) -> main Isim vs Csim analysis
%    PowerRef / PowerRefZ 1st pair member (Isim: presentation 1; Csim and
%                         Bsim: exemplar 1 of A) -> control comparison of
%                         all three types on the same footing
%  Both z-scored within person over the rows where they are defined.
T_pair = make_pair_table(T_item, T_cat, T_bcat);
save(fullfile(inDir, [variantLabel '_singleTrial_pairs_fullSample.mat']), 'T_pair', '-v7.3');
writetable(T_pair, fullfile(inDir, [variantLabel '_singleTrial_pairs_fullSample.txt']));

ID = subjects(:); %#ok<NASGU>
save(fullfile(inDir, [variantLabel '_alphaIspec_fullSample.mat']), 'ID', 'alphaIspec', 'nTrialsAlpha');
writetable(table(subjects(:), alphaIspec, nTrialsAlpha, 'VariableNames', {'ID','alphaIspec','nTrialsAlpha'}), ...
    fullfile(inDir, [variantLabel '_alphaIspec_fullSample.txt']));

fprintf('Aggregated %d subjects (within-item: %d rows, within-cat: %d rows, joint: %d rows, pairs: %d rows) -> %s\n', ...
    length(subjects), height(T_item), height(T_cat), height(T_joint), height(T_pair), inDir);

end


function T = make_pair_table(T_item, T_cat, T_bcat)
nI = height(T_item); nC = height(T_cat);
if isempty(T_bcat), nB = 0; else, nB = height(T_bcat); end
ID       = [T_item.ID; T_cat.ID];
catLabel = [T_item.catLabel; T_cat.catLabel];
itemLabel = [T_item.itemLabel; repmat({''}, nC, 1)];   % Csim pairs: two items -> no single item label
Type     = [repmat({'Isim'}, nI, 1); repmat({'Csim'}, nC, 1)];
TypeNum  = [0.5*ones(nI,1); -0.5*ones(nC,1)];
Sim      = [T_item.Sim; T_cat.Sim];
Power    = [T_item.powAVG; T_cat.powAVG];
PowerRef = [T_item.pow1st; T_cat.pow1st];
itemRep  = [T_item.itemRep; T_cat.itemRep];
catSize  = [T_item.catSize; T_cat.catSize];
Age      = [T_item.Age; T_cat.Age];
if nB > 0   % between-category pairs (catLabel = category of exemplar A)
    ID = [ID; T_bcat.ID]; catLabel = [catLabel; T_bcat.catLabel]; itemLabel = [itemLabel; repmat({''}, nB, 1)];
    Type = [Type; repmat({'Bsim'}, nB, 1)]; TypeNum = [TypeNum; nan(nB,1)];
    Sim = [Sim; T_bcat.Sim]; Power = [Power; nan(nB,1)]; PowerRef = [PowerRef; T_bcat.pow1st];
    itemRep = [itemRep; T_bcat.itemRep]; catSize = [catSize; T_bcat.catSize]; Age = [Age; T_bcat.Age];
end
PowerZ = zwithin(Power, ID);
PowerRefZ = zwithin(PowerRef, ID);
T = table(ID, catLabel, itemLabel, Type, TypeNum, Sim, Power, PowerZ, PowerRef, PowerRefZ, itemRep, catSize, Age);
end


function z = zwithin(x, ID)
% z-score within person over the rows where x is defined
z = nan(size(x));
for id = unique(ID)'
    r = ID == id & ~isnan(x);
    z(r) = (x(r) - mean(x(r))) / std(x(r));
end
end

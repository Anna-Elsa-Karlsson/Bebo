function add_trialinfoB_one_subject(subj)
% Adds trialinfo2_B (second trial of each pair) + itemLabel2 to one
% subject's original sim_mat files and re-saves them uncompressed.
% indtrials = trialinfo2 rows of sort([pairs(:,1); pairs(:,2)]), so the
% B trial of pair p is the indtrials row with trial number pairs(p,2).

subj = num2str(subj);
rsaDir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results', 'Similarity_Data', 'Ispec_30to100Hz');

for cond = {'within-item', 'within-cat'}
    f = fullfile(rsaDir, [subj '_' cond{1} '_sim_ST.mat']);
    L = load(f, 'sim_mat'); sim_mat = L.sim_mat; clear L

    trialNr = sort([sim_mat.pairs(:,1); sim_mat.pairs(:,2)]);   % trial number of each indtrials row
    [~, kB] = ismember(sim_mat.pairs(:,2), trialNr);             % row of each B trial in indtrials

    sim_mat.trialinfo2_B = sim_mat.indtrials(kB,:);
    sim_mat.itemLabel2   = sim_mat.trialinfo2_B(:,7);

    tmp = [f(1:end-4) '_tmp.mat'];
    save(tmp, 'sim_mat', '-v7.3', '-nocompression');
    movefile(tmp, f, 'f');
    fprintf('%s %s: %d pairs, trialinfo2_B added\n', subj, cond{1}, size(sim_mat.pairs,1));
end

end

function filter_isim_exemplar12_one_subject(subj)
% Keeps only exemplar 1-2 item pairs in one subject's within-item sim_mat
% (drops exemplars 3/4, so Isim matches Csim's exemplar-1-vs-2 pairs) and
% re-saves it uncompressed. Run AFTER add_trialinfoB_one_subject.

subj = num2str(subj);
rsaDir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results', 'Similarity_Data', 'Ispec_30to100Hz');
f = fullfile(rsaDir, [subj '_within-item_sim_ST.mat']);

L = load(f, 'sim_mat'); sim_mat = L.sim_mat; clear L
keep = cell2mat(sim_mat.trialinfo2(:,8)) <= 2;   % exemplar of each pair (same item on both trials)
nBefore = numel(keep);

% fields with one row per pair
for fld = {'powspctrm','powspctrmBL','powspctrmOrig','trialinfo2','trialinfo2_B','itemLabel2','pairs','trialinfo'}
    if isfield(sim_mat, fld{1})
        x = sim_mat.(fld{1});
        sim_mat.(fld{1}) = x(keep,:,:,:);
    end
end
% indtrials: trials of all pairs -> keep exemplar 1-2 trials only
sim_mat.indtrials = sim_mat.indtrials(cell2mat(sim_mat.indtrials(:,8)) <= 2, :);

tmp = [f(1:end-4) '_tmp.mat'];
save(tmp, 'sim_mat', '-v7.3', '-nocompression');
movefile(tmp, f, 'f');
fprintf('%s within-item: kept %d of %d pairs (exemplars 1-2)\n', subj, sum(keep), nBefore);

end

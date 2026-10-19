function resave_sim_mats(toSingle)
% Loops over all subjects in ids_ALL('all') and both Ispec conditions
% (within-item, within-cat) and re-saves each subject-level sim_mat file
% without compression (see resave_one_sim_mat.m) -- faster to load in
% every later step, same content.
%
%   toSingle : (optional, default false) also convert the similarity
%              fields to single (halves file size; ~7 significant digits)
%
% Usage on the grid (from RSA_pipeline/):
%   matlab -r "resave_sim_mats; exit"          % uncompressed, keep double
%   matlab -r "resave_sim_mats(true); exit"    % uncompressed + single

if nargin < 1, toSingle = false; end

pipelineRoot = fileparts(mfilename('fullpath'));
addpath(pipelineRoot); % ids_ALL.m
rsaDir = fullfile(pipelineRoot, '..', 'Results', 'Similarity_Data', 'Ispec_30to100Hz');

subjects = ids_ALL('all');
conds = {'within-item', 'within-cat'};
tAll = tic;
for s = 1:numel(subjects)
    for c = 1:numel(conds)
        f = fullfile(rsaDir, sprintf('%d_%s_sim_ST.mat', subjects(s), conds{c}));
        fprintf('[%d/%d] ', s, numel(subjects));
        resave_one_sim_mat(f, toSingle);
    end
end
fprintf('Done: %d subjects x %d conditions in %.1f h\n', numel(subjects), numel(conds), toc(tAll)/3600);

end

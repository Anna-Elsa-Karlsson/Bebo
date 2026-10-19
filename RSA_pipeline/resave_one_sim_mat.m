function resave_one_sim_mat(f, toSingle)
% Loads one sim_mat file and saves it again WITHOUT compression (-v7.3
% -nocompression), optionally converting the three similarity fields to
% single. The new file is written under a temporary name first and only
% then moved over the original, so an interrupted job never leaves a
% half-written file behind.
%
%   f        : full path to <ID>_<condition>_sim_ST.mat
%   toSingle : (optional, default false) convert powspctrm, powspctrmBL,
%              powspctrmOrig to single

if nargin < 2, toSingle = false; end

if ~exist(f, 'file')
    warning('resave_one_sim_mat:missing', 'File not found, skipped: %s', f);
    return
end

t0 = tic;
L = load(f, 'sim_mat');
sim_mat = L.sim_mat; clear L
tLoad = toc(t0);

if toSingle
    for fld = {'powspctrm', 'powspctrmBL', 'powspctrmOrig'}
        if isfield(sim_mat, fld{1})
            sim_mat.(fld{1}) = single(sim_mat.(fld{1}));
        end
    end
end

[fdir, fname] = fileparts(f);
tmpFile = fullfile(fdir, [fname '_resave_tmp.mat']);
t1 = tic;
save(tmpFile, 'sim_mat', '-v7.3', '-nocompression');
movefile(tmpFile, f, 'f'); % replace the original only once the new file is complete
tSave = toc(t1);

d = dir(f);
note = '';
if toSingle, note = ' (single)'; end
fprintf('%s: load %.0f s, save %.0f s, new size %.2f GB%s\n', fname, tLoad, tSave, d.bytes/1e9, note);

end

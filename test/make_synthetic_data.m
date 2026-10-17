function make_synthetic_data(outBase, YA, OA, opts)
% MAKE_SYNTHETIC_DATA  Simulated BEBO encoding data for testing the pipeline.
% Writes <outBase>/PreprocData/<ID>_Enc_allcond.mat ('data', FieldTrip raw
% with trialinfo2 as in the real files), <outBase>/Performance_All.txt,
% <outBase>/Elec/elec.mat and an IDs file for bebo_ids.
%
% Design as in the experiment: 80 categories, 20 per condition (normal,
% high-exemplar, high-repetition, high-repetition+exemplar), 2/4 exemplars
% x 2/4 presentations; repeats of an item 3-10 trials apart; next exemplar
% >= 5 trials after the previous one's last presentation. In 4-exemplar
% categories 2 exemplars are tested. Responses: old/similar/new.
%
% PLANTED EFFECTS (known ground truth):
%   item code   each item has its own gamma spectral fingerprint on the 15
%               "posterior" channels, 0.2-1.2 s; amplitude opts.itemAmp(YA/OA)
%   category    exemplars of a category share part of the fingerprint
%               (opts.catShare)
%   lag drift   a slow random walk over the session in the gamma spectrum
%               of ALL trials -> similarity rises with temporal proximity,
%               regardless of item or category (opts.driftAmp)
%   line noise  50 Hz on all channels
%   alpha SME   10-20 Hz amplitude drops 0.6-1.5 s; the drop is larger for
%               trials later called 'old' (opts.smeAmp(YA/OA))

if nargin < 4, opts = struct(); end
def = struct('itemAmp', [1.0 0.6], 'catShare', 0.4, 'driftAmp', 1.0, 'smeAmp', [0.3 0.15], 'seed', 7);
f = fieldnames(def);
for k = 1:numel(f), if ~isfield(opts, f{k}), opts.(f{k}) = def.(f{k}); end, end

here = fileparts(mfilename('fullpath'));
labels = strsplit(strtrim(fileread(fullfile(here, 'channels60.txt'))), newline)';
labels = [labels; {'IOL'; 'LHEOG'; 'RHEOG'; 'A1'}];          % removed by the pipeline
nCh = numel(labels);
post = find(ismember(labels, {'P7','P5','P3','P1','Pz','P2','P4','P6','P8','PO7','PO3','POz','PO4','PO8','Oz'}));

for d = {'PreprocData', 'Elec'}, if ~exist(fullfile(outBase, d{1}), 'dir'), mkdir(fullfile(outBase, d{1})); end, end
copyfile('/Users/anna.karlsson/Documents/CL/CL/BEBO/LATEST/Elec/elec.mat', fullfile(outBase, 'Elec', 'elec.mat'));
copyfile('/Users/anna.karlsson/Documents/CL/CL/BEBO/LATEST/Elec/layout_EEG_cap_20161018.mat', fullfile(outBase, 'Elec'));
fid = fopen(fullfile(outBase, 'ids.txt'), 'w'); fprintf(fid, '%s\n%s\n', num2str(YA), num2str(OA)); fclose(fid);

fs = 250; tAx = -2:1/fs:3.5-1/fs; nS = numel(tAx);
gF = 36:4:98;                                      % gamma components
ids = [YA(:); OA(:)]; isYA = [true(numel(YA),1); false(numel(OA),1)];
perf = zeros(numel(ids), 3);

for s = 1:numel(ids)
    rng(opts.seed + ids(s), 'twister');
    g = 1 + ~isYA(s);                              % 1 = YA, 2 = OA
    % ---- design and sequence ----
    conds = repelem({'normal','high-exemplar','high-repetition','high-repetition+exemplar'}, 20);
    conds = conds(randperm(80));
    seq = {};                                       % rows: {cat, catNr, cond, exemplar, pres}
    pending = struct('cat', {}, 'list', {}, 'next', {});
    for c = 1:80
        nEx = 2 + 2*contains(conds{c}, 'exemplar'); nPr = 2 + 2*contains(conds{c}, 'repetition');
        ev = zeros(0, 2);
        for e = 1:nEx, for p = 1:nPr, ev(end+1,:) = [e p]; end, end %#ok<AGROW>
        pending(end+1) = struct('cat', c, 'list', ev, 'next', randi(650)); %#ok<AGROW>
    end
    t = 0;
    while ~isempty(pending)
        t = t + 1;
        due = find([pending.next] <= t);
        if isempty(due), seq(end+1,:) = {'filler', 0, 'filler', 0, 0}; continue, end %#ok<AGROW>
        [~, kk] = min([pending(due).next]); k = due(kk);   % most overdue first
        e = pending(k).list(1,1); p = pending(k).list(1,2); c = pending(k).cat;
        seq(end+1,:) = {sprintf('cat%02d', c), c, conds{c}, e, p}; %#ok<AGROW>
        pending(k).list(1,:) = [];
        if isempty(pending(k).list), pending(k) = [];
        elseif pending(k).list(1,1) == e, pending(k).next = t + randi([3 10]);
        else, pending(k).next = t + randi([5 12]);
        end
    end
    nTr = size(seq, 1);
    % ---- test responses (per item) ----
    resp = repmat({''}, nTr, 1); itemLab = repmat({''}, nTr, 1);
    tested = containers.Map();
    for c = 1:80
        nEx = 2 + 2*contains(conds{c}, 'exemplar');
        te = randperm(nEx, 2);
        for e = 1:nEx
            key = sprintf('cat%02d_ex%d', c, e);
            if ismember(e, te)
                u = rand; r = 'old'; if u > .45, r = 'similar'; end, if u > .9, r = 'new'; end
                tested(key) = r;
            else
                tested(key) = '';
            end
        end
    end
    for i = 1:nTr
        if seq{i,2} == 0, continue, end
        key = sprintf('%s_ex%d', seq{i,1}, seq{i,4});
        itemLab{i} = [key '.jpg']; resp{i} = tested(key);
    end
    % ---- signals ----
    fpItem = containers.Map(); fpCat = cell(80,1);
    for c = 1:80, fpCat{c} = randn(numel(post), numel(gF)); end
    drift = cumsum(randn(nTr, numel(gF)), 1); drift = drift ./ std(drift(:)) * 0.15 * opts.driftAmp;
    env = double(tAx >= 0.2 & tAx <= 1.2);
    envA = 1 - 0.5 * double(tAx >= 0.6 & tAx <= 1.5);
    trial = cell(1, nTr); ti2 = cell(nTr, 11);
    for i = 1:nTr
        X = pink_noise(nCh, nS) * 2;
        X = X + 0.5 * sin(2*pi*50*tAx + 2*pi*rand(nCh,1));                   % line noise
        if seq{i,2} > 0
            key = itemLab{i};
            if ~isKey(fpItem, key), fpItem(key) = randn(numel(post), numel(gF)); end
            W = opts.itemAmp(g) * ((1-opts.catShare) * fpItem(key) + opts.catShare * fpCat{seq{i,2}});
            W = W + drift(i, :);                                              % lag drift
            W = max(W + 1.5, 0);                                              % amplitudes
            ph = 2*pi*rand(numel(post), numel(gF));
            G = zeros(numel(post), nS);
            for q = 1:numel(gF), G = G + W(:,q) .* sin(2*pi*gF(q)*tAx + ph(:,q)); end
            X(post, :) = X(post, :) + 0.25 * G .* env;
            old = strcmp(resp{i}, 'old');
            aAmp = 1.5 * (envA - old * opts.smeAmp(g) * double(tAx >= 0.6 & tAx <= 1.5));
            X(1:60, :) = X(1:60, :) + aAmp .* sin(2*pi*(10 + 4*rand)*tAx + 2*pi*rand(60,1));
        end
        trial{i} = X;
        ti2(i,:) = {i, 0, 0, seq{i,3}, seq{i,1}, seq{i,2}, itemLab{i}, seq{i,4}, seq{i,5}, 0, resp{i}};
    end
    data = struct('label', {labels}, 'fsample', fs, 'trial', {trial}, 'time', {repmat({tAx}, 1, nTr)});
    data.trialinfo = (1:nTr)'; data.trialinfo2 = ti2;
    save(fullfile(outBase, 'PreprocData', sprintf('%d_Enc_allcond.mat', ids(s))), 'data', '-v7.3');
    perf(s,:) = [ids(s), 0.25 + 0.1*isYA(s) + 0.05*randn, 0.3 + 0.2*isYA(s) + 0.1*randn];
    fprintf('synthetic subject %d: %d trials\n', ids(s), nTr);
end
P = array2table(perf, 'VariableNames', {'ID','IM','LDI'});
writetable(P, fullfile(outBase, 'Performance_All.txt'));
end

function X = pink_noise(nCh, nS)
F = fft(randn(nCh, nS), [], 2);
f = [1, 1:floor(nS/2), floor((nS-1)/2):-1:1];
X = real(ifft(F ./ sqrt(f(1:nS)), [], 2));
X = X ./ std(X, 0, 2);
end

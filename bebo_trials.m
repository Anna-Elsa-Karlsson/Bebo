function [T, pairsI, pairsC] = bebo_trials(trialinfo2, C)
% BEBO_TRIALS  The trial rule of the analysis plan (section 0), in one place.
%
%   T      table, one row per ANALYSED trial (presentation 1 or 2 of an item
%          with an 'old'/'similar' test response), sorted by trial number:
%            row       row in the data (trialinfo2 / powspctrm)
%            trialNr   position in the encoding sequence
%            category, item, exemplar, pres, resp
%            Acc       1 = 'old', 0 = 'similar'
%            catSize   number of exemplars in the category (2 or 4)
%            itemRep   number of presentations of the item (2 or 4)
%   pairsI [nI x 2] rows of T: presentation 1 and 2 of the same item
%   pairsC [nC x 2] rows of T: presentation 1 of two exemplars of the same
%          category, ALL unique exemplar pairs (both analysed)
%   In every pair the EARLIER trial is in column 1.

ti = trialinfo2;
c  = C.ti;
n  = size(ti, 1);

resp = string(cellfun(@to_char, ti(:, c.resp), 'UniformOutput', false));
pres = cell2mat(ti(:, c.pres));
keep = ismember(resp, C.validResp) & ismember(pres, C.presentations);

cond = string(cellfun(@to_char, ti(:, c.cond), 'UniformOutput', false));
T = table();
T.row      = (1:n)';
T.trialNr  = cell2mat(ti(:, c.trialNr));
T.category = string(cellfun(@to_char, ti(:, c.category), 'UniformOutput', false));
T.item     = string(cellfun(@to_char, ti(:, c.item), 'UniformOutput', false));
T.exemplar = cell2mat(ti(:, c.exemplar));
T.pres     = pres;
T.resp     = resp;
T.Acc      = double(resp == "old");
T.catSize  = 2 + 2 * double(contains(cond, "exemplar"));
T.itemRep  = 2 + 2 * double(contains(cond, "repetition"));
T = T(keep, :);
T = sortrows(T, 'trialNr');

if numel(unique(T.trialNr)) ~= height(T)
    error('bebo_trials:trialNr', 'Duplicate trial numbers among the analysed trials');
end

% Isim: presentation 1 and 2 of the same item
items = unique(T.item);
pairsI = zeros(0, 2);
for k = 1:numel(items)
    r1 = find(T.item == items(k) & T.pres == 1);
    r2 = find(T.item == items(k) & T.pres == 2);
    if numel(r1) == 1 && numel(r2) == 1
        pairsI(end+1, :) = [r1 r2]; %#ok<AGROW>
    elseif numel(r1) > 1 || numel(r2) > 1
        error('bebo_trials:item', 'Item %s has more than one presentation-%d trial', items(k), 1 + (numel(r1) <= 1));
    end
end

% Csim: all unique pairs of presentation-1 trials of different exemplars of
% the same category
cats = unique(T.category);
pairsC = zeros(0, 2);
for k = 1:numel(cats)
    r = find(T.category == cats(k) & T.pres == 1);   % sorted by trialNr
    if numel(unique(T.exemplar(r))) ~= numel(r)
        error('bebo_trials:exemplar', 'Category %s has two presentation-1 trials of the same exemplar', cats(k));
    end
    if numel(r) >= 2
        pairsC = [pairsC; nchoosek(r, 2)]; %#ok<AGROW>
    end
end

% earlier trial first (T is sorted by trialNr, so row order = time order)
pairsI = sort(pairsI, 2);
pairsC = sort(pairsC, 2);
end

function s = to_char(x)
if ischar(x) || isstring(x), s = char(x);
elseif isnumeric(x) && isscalar(x), s = num2str(x);
else, s = '';
end
end

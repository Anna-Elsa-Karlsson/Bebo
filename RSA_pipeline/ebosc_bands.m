function B = ebosc_bands()
% Rhythm bands for eBOSC (2026-10-05): what counts as an ALPHA vs a BETA
% RHYTHM (episode mean frequency), and the ranges left out of the 1/f
% background fit (around the empirical peak in each range). Separate from
% alpha_band() (the 8-25 Hz POWER/desync band of the main analyses, taken
% from the desync cluster): rhythms are reported separately so that faster
% (beta) episodes are not counted as alpha. Change here only; changing it
% requires rerunning the eBOSC subject step.
B.names  = {'alpha', 'beta'};
B.ranges = [8 14; 14 25];          % Hz; an episode belongs to a band if lo <= freq < hi
B.exclBG = [3 7; 8 14; 14 25];     % theta, alpha, beta peaks removed before the background fit
end

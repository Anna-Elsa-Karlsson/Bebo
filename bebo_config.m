function C = bebo_config()
% BEBO_CONFIG  All paths and fixed analysis settings in ONE place.
% Every step of the pipeline reads its settings from here, so the scripts
% cannot drift apart from each other or from the analysis plan
% (BEBO/analysis_plan_2026-10-17_FINAL.txt -- section numbers below refer
% to that file).
%
% Folder layout (grid):  ~/BEBO2/bebo-pipeline/   <- this repository
%                        ~/BEBO2/PreprocData/      <- input, <ID>_Enc_allcond.mat
%                        ~/BEBO2/Elec/elec.mat     <- electrode positions
%                        ~/BEBO2/Performance_All.txt
%                        ~/BEBO2/Results_v2/       <- ALL output (new, empty)

C.root     = fileparts(mfilename('fullpath'));    % this repository
C.base     = fileparts(C.root);                   % parent folder (~/BEBO2)
if ~isempty(getenv('BEBO_BASE')), C.base = getenv('BEBO_BASE'); end  % test runs only
C.preproc  = fullfile(C.base, 'PreprocData');
C.elecFile = fullfile(C.base, 'Elec', 'elec.mat');
C.layoutFile = fullfile(C.base, 'Elec', 'layout_EEG_cap_20161018.mat');   % topoplots only
C.perfFile = fullfile(C.base, 'Performance_All.txt');
C.out      = fullfile(C.base, 'Results_v2');

% Output subfolders
C.dir.gammaSubj = fullfile(C.out, 'gamma_subject');   % per-person RSA (s1)
C.dir.alphaSubj = fullfile(C.out, 'alpha_subject');   % per-person power (s2)
C.dir.trialLog  = fullfile(C.out, 'trial_log');       % per-person trial/pair counts + lags
C.dir.group     = fullfile(C.out, 'group');           % cluster tests (g1, g2)
C.dir.tables    = fullfile(C.out, 'tables');          % per-person + single-trial tables (g3)
C.dir.stats     = fullfile(C.out, 'stats');           % statistics text files (g4, R)
C.dir.figs      = fullfile(C.out, 'figures');         % plot scripts only

% FieldTrip (grid version first; local copy for testing)
C.fieldtrip = {'/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/', ...
               '/Users/anna.karlsson/Documents/MATLAB/fieldtrip/'};

% Channels: 60 EEG channels (EOG and mastoid removed)
C.channel = {'all', '-IOL', '-RHEOG', '-LHEOG', '-A1'};

%% 0. Trial rule
C.validResp = {'old', 'similar'};   % 'new' (and untested) excluded
C.presentations = [1 2];            % analysed trials = presentations 1 and 2
% trialinfo2 columns used (confirmed on the original data):
C.ti.trialNr  = 1;    % position in the encoding sequence (lag = difference)
C.ti.cond     = 4;    % 'normal' | 'high-exemplar' | 'high-repetition' | 'high-repetition+exemplar'
C.ti.category = 5;    % category label
C.ti.item     = 7;    % item (image) label
C.ti.exemplar = 8;    % exemplar number within the category
C.ti.pres     = 9;    % presentation number of the item
C.ti.resp     = 11;   % test response: 'old' | 'similar' | 'new' | (untested)

%% 1a. Gamma RSA
C.gamma.dftfreq   = [50 100];           % line noise removed before the TFR
C.gamma.foi       = 35:2:100;           % Hz (35 Hz bin spans 30-40 Hz)
C.gamma.tapsmofrq = 5;                  % +/- Hz -> 3 DPSS tapers with 0.4 s
C.gamma.t_ftimwin = 0.4;                % s
C.gamma.toi       = 0:0.02:1.7;         % s, 20 ms steps
C.gamma.pad       = 6;                  % s
C.gamma.boscBand  = [35 100];           % aperiodic fit on the analysis band
% Reference matched on lag and presentation: ALL unrelated pairs with the
% same lag and presentation structure
C.lag.minRef   = 20;     % if fewer than 20 exist at the exact lag ...
C.lag.maxWiden = 2;      % ... add lags +/-1, then +/-2

%% 1b. Alpha SME + single-trial power
C.alpha.foi    = 2:2:30;                % Hz
C.alpha.width  = 7;                     % Morlet cycles
C.alpha.toi    = -0.28:0.04:1.68;       % s; on the sample grid, covers -0.3 to 1.7
C.alpha.pad    = 6;                     % s (> 5.5 s epochs; 1/6 Hz grid -> exact 2:2:30 Hz)
C.alpha.glmBand = [8 25];               % Hz, a priori band for the glmer (2)

%% Cluster permutation tests (1a, 1b)
C.clus.clusteralpha     = 0.05;
C.clus.minnbchan        = 1;
C.clus.numrandomization = 1000;
C.clus.alpha            = 0.025;        % two-tailed
C.clus.neighbourdist    = 55;           % mm
C.clus.smeFreq          = [2 30];
C.clus.smeLatency       = [-0.3 1.7];

%% Random seeds (everything random is reproducible)
C.seed.clusterI  = 101;    % Ispec permutation test
C.seed.clusterS  = 102;    % SME permutation test
C.seed.splitHalf = 103;    % split-half category split

end

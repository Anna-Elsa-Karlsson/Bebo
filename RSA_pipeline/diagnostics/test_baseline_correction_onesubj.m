function test_baseline_correction_onesubj(subj, freqLo, freqHi, bgFitMode)
% DIAGNOSTIC ONLY -- does not touch the production pipeline or its output
% files. Compares, for ONE subject's within-item pairs, two things at
% once (2026-08-27, see session notes):
%
%   (1) METHOD: how the "baseline" reference term is computed --
%       OLD (current production code): average the baseline vector across
%         time FIRST, then correlate that single denoised vector against
%         activation. Concern: this reduces noise on only one side of the
%         correlation (attenuation-due-to-unreliability), which could
%         inflate the baseline term relative to the real activation-
%         activation correlation (always single-timepoint-vs-single-
%         timepoint, equally noisy on both sides).
%       NEW (proposed fix): correlate activation against EVERY individual
%         raw baseline timepoint separately (matched noise level), THEN
%         average the resulting correlation VALUES. Same "flat across
%         retrieval-time" output shape as OLD, just averaged in the other
%         order.
%
%   (2) WINDOW: which pre-stimulus window counts as "baseline" --
%       NARROW: [-1 -0.3] s, the current production window.
%       WIDE:   [-1.5 0] s, added after realizing the actual trial timing
%         (1500ms stimulus, 1500-2000ms jittered ITI) leaves much more
%         usable clean pre-stimulus time than the narrow window uses.
%
% Reports all 2x2 combinations (narrow/wide x old/new) plus the real
% activation-activation correlation and each combination's resulting
% corrected similarity (real - baseline), across all within-item pairs x
% channels x time-time cells for this one subject. Saves the flattened
% comparison arrays to
% Results/BaselineDiagnostic/<subj>_baseline_comparison.mat.
%
% Usage: matlab -r "test_baseline_correction_onesubj('3109','30','100','legacy'); exit"
% Defaults to Ispec broadband / legacy bgFitMode (the variant meant to
% match the already-published, known-significant result) if freq/mode
% args are omitted.
%
% NOTE: assumes the precomputed '_Enc_tfrBL_2_125Hz_ST.mat' file actually
% covers at least [-1.5 0]s pre-stimulus -- if it doesn't extend that far,
% the wide-window ft_selectdata call below will error, telling us the
% file's actual available range.

if nargin < 1 || isempty(subj), error('subj is required, e.g. ''3109'''); end
if nargin < 2 || isempty(freqLo), freqLo = 30; end
if nargin < 3 || isempty(freqHi), freqHi = 100; end
if nargin < 4 || isempty(bgFitMode), bgFitMode = 'legacy'; end
if ischar(freqLo), freqLo = str2double(freqLo); end
if ischar(freqHi), freqHi = str2double(freqHi); end
subj = num2str(subj);

pipelineRoot = fileparts(mfilename('fullpath'));
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/');
ft_defaults

pdat.tfr   = fullfile(pipelineRoot, '..', 'TFRdata', 'Time-frequencyEEGdata');
pdat.tfrBL = fullfile(pipelineRoot, '..', 'TFRdat_BL');

if strcmp(bgFitMode, 'legacy')
    bgFitBand = [freqLo freqHi];
elseif freqHi <= 20
    bgFitBand = [2 20];
else
    bgFitBand = [25 125];
end
analysisBand = [freqLo freqHi];

fprintf('Loading + BOSC-correcting subject %s (%d-%dHz, bgFitMode=%s)...\n', subj, freqLo, freqHi, bgFitMode);
Data_logpow = load_activation(subj, pdat, bgFitBand, analysisBand);

% Trimmed activation, [0 1.5]s -- matches stimulus-on-screen duration AND
% (deliberately) the wide baseline window's duration ([-1.5 0]s = 1.5s
% too), so a literal cell-by-cell subtraction between the raw activ-vs-
% activ matrix and the raw activ-vs-baseline matrix becomes shape-valid,
% not just a pooled distributional comparison (2026-08-27, per user
% request). Post-hoc latency restriction on the already-BOSC-corrected
% Data_logpow -- no re-fitting, just cutting off the last 0.2s.
cfg = []; cfg.latency = [0 1.5];
Data_logpow_trim = ft_selectdata(cfg, Data_logpow);
Data_logpow_trim.trialinfo2 = Data_logpow.trialinfo2;

tfr_bl_raw = load_baseline_raw(subj, pdat);
fprintf('Correcting baseline, narrow window [-1 -0.3]s...\n');
Data_logpow_bl_narrow = correct_baseline_window(tfr_bl_raw, bgFitBand, analysisBand, [-1 -0.3]);
fprintf('Correcting baseline, wide window [-1.5 0]s...\n');
Data_logpow_bl_wide = correct_baseline_window(tfr_bl_raw, bgFitBand, analysisBand, [-1.5 0]);
clear tfr_bl_raw

nTimeTrim = length(Data_logpow_trim.time);
nTimeBL_wide_check = length(Data_logpow_bl_wide.time);
if nTimeTrim ~= nTimeBL_wide_check
    error(['Trimmed activation [0 1.5]s has %d timepoints but wide baseline ' ...
        '[-1.5 0]s has %d -- the two TFR files use different sampling grids, ' ...
        'so a literal cell-by-cell subtraction is not valid as-is.'], nTimeTrim, nTimeBL_wide_check);
end

%% Within-item pairs -- identical logic to rsa_subject_level.m's pairs_within_item_vs_cat

pres1 = 1;
pairs = [];
for stim1 = 1:size(Data_logpow.powspctrm,1)
    if ~ismember(stim1,pairs) && Data_logpow.trialinfo2{stim1,9} == pres1 && ~isempty(find(strcmp(Data_logpow.trialinfo2{stim1,11},'old') | strcmp(Data_logpow.trialinfo2{stim1,11},'similar')))
        curr_pair = Data_logpow.trialinfo2{stim1,7};
        if strcmp(Data_logpow.trialinfo2(stim1,4),'normal') || strcmp(Data_logpow.trialinfo2(stim1,4),'high-exemplar')
            if find(strcmp(Data_logpow.trialinfo2(:,7),curr_pair) & cell2mat(Data_logpow.trialinfo2(:,9)) == 2)
                pairs(stim1,2) = find(strcmp(Data_logpow.trialinfo2(:,7),curr_pair) & cell2mat(Data_logpow.trialinfo2(:,9)) == 2);
                pairs(stim1,1) = stim1;
            end
        elseif strcmp(Data_logpow.trialinfo2(stim1,4),'high-repetition+exemplar') || strcmp(Data_logpow.trialinfo2(stim1,4),'high-repetition')
            if find(strcmp(Data_logpow.trialinfo2(:,7),curr_pair) & cell2mat(Data_logpow.trialinfo2(:,9)) == 2)
                pairs(stim1,1) = stim1;
                pairs(stim1,2) = find(strcmp(Data_logpow.trialinfo2(:,7),curr_pair) & cell2mat(Data_logpow.trialinfo2(:,9)) == 2);
            end
        end
    end
end
pairs(~all(pairs,2),:) = [];

nPairs = size(pairs,1);
nChan = size(Data_logpow.powspctrm,2);
nTimeEnc = length(Data_logpow.time);
nTimeBL_narrow = length(Data_logpow_bl_narrow.time);
nTimeBL_wide = length(Data_logpow_bl_wide.time);
fprintf('Subject %s: %d within-item pairs, %d channels, %d encoding timepoints (baseline: %d narrow / %d wide timepoints)\n', ...
    subj, nPairs, nChan, nTimeEnc, nTimeBL_narrow, nTimeBL_wide);

%% Compute all 2x2 (window x method) baseline variants for every pair x channel,
%  PLUS the fully-raw (no averaging on either side, at any point)
%  activation-vs-baseline distribution, to directly check whether raw
%  single-timepoint-vs-single-timepoint activ-bl correlations already
%  look comparable in scale to raw single-timepoint-vs-single-timepoint
%  activ-activ correlations (the "real" row) -- i.e. is there already an
%  activ-vs-bl difference even before any averaging/denoising is applied
%  to either side?

nCellsPerCorr = nTimeEnc * nTimeEnc;
totalCells = nPairs * nChan * nCellsPerCorr;
real_all       = nan(totalCells,1,'single');
bl_narrowOld   = nan(totalCells,1,'single');
bl_narrowNew   = nan(totalCells,1,'single');
bl_wideOld     = nan(totalCells,1,'single');
bl_wideNew     = nan(totalCells,1,'single');
idx = 0;

totalCells_narrowRaw = nPairs * nChan * nTimeEnc * nTimeBL_narrow;
totalCells_wideRaw   = nPairs * nChan * nTimeEnc * nTimeBL_wide;
bl_narrowRaw = nan(totalCells_narrowRaw,1,'single');
bl_wideRaw   = nan(totalCells_wideRaw,1,'single');
idxNR = 0;
idxWR = 0;

% Literal, shape-matched subtraction: trimmed [0 1.5]s activation vs
% activation (real_trim), and trimmed activation vs raw wide baseline
% (bl_wideRaw_trim) -- both [nTimeTrim x nTimeTrim], so literalCorrected
% = real_trim - bl_wideRaw_trim is a genuine cell-by-cell subtraction,
% not a pooled-distribution comparison like the RAW rows above.
totalCells_trim = nPairs * nChan * nTimeTrim * nTimeTrim;
real_trim_all      = nan(totalCells_trim,1,'single');
bl_wideRaw_trim    = nan(totalCells_trim,1,'single');
idxTrim = 0;

for curr_pair = 1:nPairs
    for chan = 1:nChan
        tf1 = squeeze(Data_logpow.powspctrm(pairs(curr_pair,1),chan,:,:));
        tf2 = squeeze(Data_logpow.powspctrm(pairs(curr_pair,2),chan,:,:));
        r_real = spectral_rsa([],tf1,tf2,Data_logpow.time);

        tf2_narrow = squeeze(Data_logpow_bl_narrow.powspctrm(pairs(curr_pair,2),chan,:,:));
        tf2_wide   = squeeze(Data_logpow_bl_wide.powspctrm(pairs(curr_pair,2),chan,:,:));

        % narrow window, OLD method (avg vector, correlate once)
        tf2new = repmat(mean(tf2_narrow,2,'omitnan'),1,size(tf1,2));
        r_narrowOld = spectral_rsa([],tf1,tf2new,Data_logpow.time);

        % narrow window: fully raw activ-vs-baseline correlation matrix,
        % NO averaging on either side -- [nTimeEnc x nTimeBL_narrow].
        % NEW method reuses this by averaging across the baseline-time
        % dimension (dim 2) afterward.
        r_raw_narrow = spectral_rsa([],tf1,tf2_narrow,Data_logpow.time);
        r_narrowNew = repmat(mean(r_raw_narrow,2,'omitnan'),1,size(tf1,2));

        % wide window, OLD method
        tf2new = repmat(mean(tf2_wide,2,'omitnan'),1,size(tf1,2));
        r_wideOld = spectral_rsa([],tf1,tf2new,Data_logpow.time);

        % wide window: fully raw activ-vs-baseline correlation matrix
        r_raw_wide = spectral_rsa([],tf1,tf2_wide,Data_logpow.time);
        r_wideNew = repmat(mean(r_raw_wide,2,'omitnan'),1,size(tf1,2));

        n = numel(r_real);
        real_all(idx+1:idx+n)     = r_real(:);
        bl_narrowOld(idx+1:idx+n) = r_narrowOld(:);
        bl_narrowNew(idx+1:idx+n) = r_narrowNew(:);
        bl_wideOld(idx+1:idx+n)   = r_wideOld(:);
        bl_wideNew(idx+1:idx+n)   = r_wideNew(:);
        idx = idx + n;

        nNR = numel(r_raw_narrow);
        bl_narrowRaw(idxNR+1:idxNR+nNR) = r_raw_narrow(:);
        idxNR = idxNR + nNR;

        nWR = numel(r_raw_wide);
        bl_wideRaw(idxWR+1:idxWR+nWR) = r_raw_wide(:);
        idxWR = idxWR + nWR;

        % Literal shape-matched subtraction: trimmed activation vs
        % trimmed activation, and trimmed activation vs raw wide baseline
        tf1_trim = squeeze(Data_logpow_trim.powspctrm(pairs(curr_pair,1),chan,:,:));
        tf2_trim = squeeze(Data_logpow_trim.powspctrm(pairs(curr_pair,2),chan,:,:));
        r_real_trim = spectral_rsa([],tf1_trim,tf2_trim,Data_logpow_trim.time);
        r_bl_wide_trim = spectral_rsa([],tf1_trim,tf2_wide,Data_logpow_trim.time);

        nT = numel(r_real_trim);
        real_trim_all(idxTrim+1:idxTrim+nT)   = r_real_trim(:);
        bl_wideRaw_trim(idxTrim+1:idxTrim+nT) = r_bl_wide_trim(:);
        idxTrim = idxTrim + nT;
    end
    if mod(curr_pair,10)==0, fprintf('  processed %d/%d pairs\n', curr_pair, nPairs); end
end

%% Summary

corr_narrowOld = real_all - bl_narrowOld;
corr_narrowNew = real_all - bl_narrowNew;
corr_wideOld   = real_all - bl_wideOld;
corr_wideNew   = real_all - bl_wideNew;
literalCorrected_wide = real_trim_all - bl_wideRaw_trim;

fprintf('\n=== BASELINE CORRECTION COMPARISON (subject %s, %d-%dHz, %s) ===\n', subj, freqLo, freqHi, bgFitMode);
fprintf('%-30s %10s %10s %10s\n', '', 'mean', 'SD', 'median');
report_row('real (activ-activ)', real_all);
report_row('baseline narrow-RAW (no avg either side)', bl_narrowRaw);
report_row('baseline wide-RAW (no avg either side)',   bl_wideRaw);
report_row('baseline narrow-OLD', bl_narrowOld);
report_row('baseline narrow-NEW', bl_narrowNew);
report_row('baseline wide-OLD',   bl_wideOld);
report_row('baseline wide-NEW',   bl_wideNew);
report_row('corrected narrow-OLD', corr_narrowOld);
report_row('corrected narrow-NEW', corr_narrowNew);
report_row('corrected wide-OLD',   corr_wideOld);
report_row('corrected wide-NEW',   corr_wideNew);

fprintf('\n--- Noise-comparability check (no averaging on either side, at any point) ---\n');
fprintf('Mean difference: real - baseline-narrow-RAW = %.4f\n', mean(real_all,'omitnan') - mean(bl_narrowRaw,'omitnan'));
fprintf('Mean difference: real - baseline-wide-RAW   = %.4f\n', mean(real_all,'omitnan') - mean(bl_wideRaw,'omitnan'));
fprintf('(SD comparison -- real SD=%.4f vs narrow-RAW SD=%.4f vs wide-RAW SD=%.4f; similar SDs would confirm comparable noise levels)\n', ...
    std(real_all,'omitnan'), std(bl_narrowRaw,'omitnan'), std(bl_wideRaw,'omitnan'));

fprintf('\n--- Method effect (averaging order) vs window effect ---\n');
fprintf('Mean baseline differences (method, within narrow window): OLD-NEW = %.4f\n', mean(bl_narrowOld - bl_narrowNew, 'omitnan'));
fprintf('Mean baseline differences (method, within wide window):   OLD-NEW = %.4f\n', mean(bl_wideOld - bl_wideNew, 'omitnan'));
fprintf('Mean baseline differences (window, within OLD method):    narrow-wide = %.4f\n', mean(bl_narrowOld - bl_wideOld, 'omitnan'));
fprintf('Mean baseline differences (window, within NEW method):    narrow-wide = %.4f\n', mean(bl_narrowNew - bl_wideNew, 'omitnan'));

fprintf('\n--- LITERAL cell-by-cell subtraction, shape-matched ([0 1.5]s activation vs [-1.5 0]s baseline, %d timepoints each) ---\n', nTimeTrim);
report_row('real (trimmed activ-activ)', real_trim_all);
report_row('baseline (trimmed activ vs raw wide bl)', bl_wideRaw_trim);
report_row('literal corrected (real - bl, per cell)', literalCorrected_wide);

outDir = fullfile(pipelineRoot, 'Results', 'BaselineDiagnostic');
if ~exist(outDir,'dir'), mkdir(outDir); end
save(fullfile(outDir, [subj '_baseline_comparison.mat']), ...
    'real_all','bl_narrowRaw','bl_wideRaw','bl_narrowOld','bl_narrowNew','bl_wideOld','bl_wideNew', ...
    'corr_narrowOld','corr_narrowNew','corr_wideOld','corr_wideNew', ...
    'real_trim_all','bl_wideRaw_trim','literalCorrected_wide','pairs','-v7.3');
fprintf('\nSaved to %s\n', fullfile(outDir, [subj '_baseline_comparison.mat']));

end


function report_row(label, x)
fprintf('%-30s %10.4f %10.4f %10.4f\n', label, mean(x,'omitnan'), std(x,'omitnan'), median(x,'omitnan'));
end


function Data_logpow = load_activation(subj, pdat, bgFitBand, analysisBand)
% IDENTICAL to the activation-loading half of rsa_subject_level.m's
% local load_bosc_correct_restrict -- copied verbatim so this diagnostic
% uses exactly the same activation data prep as production.

disp('loading activation TFR...')
load(fullfile(pdat.tfr, [subj '_Enc_tfr_2_125Hz_ST.mat']), 'tfr')

cfg = [];
cfg.latency = [0 1.7];
cfg.frequency = bgFitBand;
tfr = ft_selectdata(cfg, tfr);

cfg = []; cfg.parameter = 'powspctrm'; cfg.operation = 'log10';
warning('off', 'stats:statrobustfit:IterationLimit')
disp('running BOSC (activation)...')
Data_c = subtract_mean_noise_spectrum(cfg, tfr);
clear tfr

ft_warning off 'FieldTrip:getdimord:warning_dimord_could_not_be_determined'
Data_logpow_wide = ft_math(cfg, Data_c);
Data_logpow_wide.trialinfo2 = Data_c.trialinfo2;
clear Data_c

cfg = []; cfg.frequency = analysisBand;
Data_logpow = ft_selectdata(cfg, Data_logpow_wide);
Data_logpow.trialinfo2 = Data_logpow_wide.trialinfo2;

end


function tfr_bl_raw = load_baseline_raw(subj, pdat)
% Loads the precomputed baseline TFR file ONCE, unclipped -- both the
% narrow and wide window corrections below select their own latency
% range from this same loaded structure, so the file is only read from
% disk one time.

disp('loading baseline TFR (raw, unclipped)...')
load(fullfile(pdat.tfrBL, [subj '_Enc_tfrBL_2_125Hz_ST.mat']), 'tfr')
tfr_bl_raw = tfr;

end


function Data_logpow_bl = correct_baseline_window(tfr_bl_raw, bgFitBand, analysisBand, baseLatency)
% Same BOSC-correction steps as rsa_subject_level.m's baseline half, just
% parameterized on baseLatency so it can be called once for the narrow
% [-1 -0.3]s window and once for the wide [-1.5 0]s window.

cfg = [];
cfg.latency = baseLatency;
cfg.frequency = bgFitBand;
tfr = ft_selectdata(cfg, tfr_bl_raw);

cfg = []; cfg.parameter = 'powspctrm'; cfg.operation = 'log10';
warning('off', 'stats:statrobustfit:IterationLimit')
disp(['running BOSC (baseline, latency ' mat2str(baseLatency) ')...'])
Data_c = subtract_mean_noise_spectrum(cfg, tfr);
clear tfr

ft_warning off 'FieldTrip:getdimord:warning_dimord_could_not_be_determined'
Data_logpow_bl_wide = ft_math(cfg, Data_c);
Data_logpow_bl_wide.trialinfo2 = Data_c.trialinfo2;
clear Data_c

cfg = []; cfg.frequency = analysisBand;
Data_logpow_bl = ft_selectdata(cfg, Data_logpow_bl_wide);
Data_logpow_bl.trialinfo2 = Data_logpow_bl_wide.trialinfo2;

end

function PF_group_level(pdat)
% Group-level PF step: aggregates per-subject raw spectra
% (PF_subject_level.m), background-corrects, detects peaks (Corcoran et
% al. 2018 method), runs group stats, and makes plots. Ported unchanged
% in method from BEBO_PF_analysis.m Sections 2-6 -- only the search band
% source has changed (loaded from PF_roi_config.mat, derived from the
% alpha_pre2post cluster, instead of the old hardcoded SME-derived 6-16Hz).
%
%   pdat : struct with fields .out (per-subject *_PFraw.mat files +
%          this step's output), .roiDir (contains PF_roi_config.mat),
%          .plotDir

ids = ids_ALL('all');
nAll = length(ids);
ID_ya = ids_ALL('ya'); nYA = length(ID_ya);
ID_oa = ids_ALL('oa'); nOA = length(ID_oa);
idx_YA = 1:nYA;
idx_OA = nYA+1:nAll;
age_vec = [ones(nYA,1); 2*ones(nOA,1)];

col_YA = [30 144 255]/255;
col_OA = [80  80  80]/255;

L = load(fullfile(pdat.roiDir, 'PF_roi_config.mat'), 'freqBand');
searchBand = L.freqBand; % = alpha_band() (8-16 Hz), set in PF_derive_roi.m -> flanks [2 7] and [17 20] Hz

foiGrid = 2:0.25:20;
nFreq = length(foiGrid);
% Flank bands (background fit) exclude the search band with a 1Hz buffer
% on each side, clipped to foiGrid's own extent -- same buffer logic as
% the original design (search=[6,16], flank=[2,5] U [17,20]).
flankBands = [2, max(2, searchBand(1)-1); min(20, searchBand(2)+1), 20];

sgFrameLen = 11;
sgOrder = 5;

%% Load + stack per-subject raw spectra
spec_raw = nan(nAll, nFreq);
nTrials = nan(nAll,1);
for j = 1:nAll
    subj = num2str(ids(j));
    f = fullfile(pdat.out, [subj '_PFraw.mat']);
    if ~exist(f,'file')
        warning('Missing PF raw spectrum for subject %s -- skipping', subj);
        continue
    end
    L = load(f, 'spec_raw', 'nTrials');
    spec_raw(j,:) = L.spec_raw;
    nTrials(j) = L.nTrials;
end

%% Background (aperiodic) correction -- fit only on flanking bins
spec_corrected = nan(size(spec_raw));
bgSlope = nan(nAll,1);
bgIntercept = nan(nAll,1);
flankMask = (foiGrid >= flankBands(1,1) & foiGrid <= flankBands(1,2)) | ...
            (foiGrid >= flankBands(2,1) & foiGrid <= flankBands(2,2));

for j = 1:nAll
    if all(isnan(spec_raw(j,:))), continue; end
    logPow = log(spec_raw(j,:));
    logFreq = log(foiGrid);
    b = robustfit(logFreq(flankMask), logPow(flankMask));
    bgIntercept(j) = b(1);
    bgSlope(j) = b(2);
    bgFit = b(1) + b(2) * logFreq;
    spec_corrected(j,:) = exp(logPow - bgFit);
end

%% Peak / CoG detection (Corcoran et al., 2018)
searchIdx = foiGrid >= searchBand(1) & foiGrid <= searchBand(2);
freqSearch = foiGrid(searchIdx);

PF_peak = nan(nAll,1);
PF_CoG = nan(nAll,1);
PF_final = nan(nAll,1);
peakDetected = false(nAll,1);
peakProminence = nan(nAll,1);

for j = 1:nAll
    y = spec_corrected(j,:);
    if all(isnan(y)), continue; end

    ySmooth = sgolayfilt(y, sgOrder, sgFrameLen);
    ySearch = ySmooth(searchIdx);

    noiseFloor = std(y(flankMask) - 1);
    minProm = 2 * noiseFloor;

    % local_findpeaks instead of findpeaks (2026-09-26): findpeaks loads
    % libmwdastudio.so, which fails on grid nodes without libXcomposite
    % ("Can't load ... libmwdastudio.so"). Same peak + prominence
    % definition as findpeaks.
    [pks, locs, proms] = local_findpeaks(ySearch, freqSearch, minProm); %#ok<ASGLU>

    if ~isempty(pks)
        [~, iBest] = max(proms);
        peakDetected(j) = true;
        peakProminence(j) = proms(iBest);

        [~, iFull] = min(abs(foiGrid - locs(iBest)));
        if iFull > 1 && iFull < length(foiGrid)
            y_1 = ySmooth(iFull-1); y0 = ySmooth(iFull); y_p1 = ySmooth(iFull+1);
            df = foiGrid(2) - foiGrid(1);
            offset = 0.5 * (y_1 - y_p1) / (y_1 - 2*y0 + y_p1);
            PF_peak(j) = foiGrid(iFull) + offset * df;
        else
            PF_peak(j) = locs(iBest);
        end
    end

    PF_CoG(j) = sum(freqSearch .* ySearch) / sum(ySearch);

    if peakDetected(j)
        PF_final(j) = PF_peak(j);
    else
        PF_final(j) = PF_CoG(j);
    end
end

%% Save per-subject table
method = repmat({'CoG'}, nAll, 1);
method(peakDetected) = {'peak'};
AgeGroup = repmat({'YA'}, nAll, 1);
AgeGroup(idx_OA) = {'OA'};

T = table(ids', AgeGroup, nTrials, PF_final, method, PF_peak, PF_CoG, ...
    peakProminence, bgSlope, bgIntercept, ...
    'VariableNames', {'ID','AgeGroup','nTrials','PF_final','method','PF_peak','PF_CoG', ...
                       'peakProminence','bgSlope','bgIntercept'});
if ~exist(pdat.out,'dir'), mkdir(pdat.out); end
writetable(T, fullfile(pdat.out,'BEBO_PF_estimates.csv'));
save(fullfile(pdat.out,'BEBO_PF_estimates.mat'), 'T', 'spec_raw', 'spec_corrected', 'foiGrid', 'searchBand');
fprintf('Per-subject PF table saved to %s\n', fullfile(pdat.out,'BEBO_PF_estimates.csv'));

%% Diagnostics + group stats
fid_pf = fopen(fullfile(pdat.out,'BEBO_PF_stats.txt'),'w');
fprintf(fid_pf,'BEBO Peak Alpha Frequency -- Statistical Results\nGenerated: %s\nSearch band (fixed alpha band, alpha_band.m): [%.1f %.1f]Hz\n\n', datestr(now), searchBand(1), searchBand(2));

[~, p_trl, ~, stats_trl] = ttest2(nTrials(idx_YA), nTrials(idx_OA));
fprintf(fid_pf,'=== DIAGNOSTIC: TRIAL COUNTS ===\n\n');
fprintf(fid_pf,'Younger: M=%.1f SD=%.1f (min=%d)  Older: M=%.1f SD=%.1f (min=%d)\n', ...
    mean(nTrials(idx_YA)), std(nTrials(idx_YA)), min(nTrials(idx_YA)), ...
    mean(nTrials(idx_OA)), std(nTrials(idx_OA)), min(nTrials(idx_OA)));
fprintf(fid_pf,'t(%d)=%.3f, p=%.4f\n\n', stats_trl.df, stats_trl.tstat, p_trl);

detTbl = [sum(peakDetected(idx_YA)), sum(~peakDetected(idx_YA)); ...
          sum(peakDetected(idx_OA)), sum(~peakDetected(idx_OA))];
[~, p_det] = fishertest(detTbl);
fprintf(fid_pf,'=== DIAGNOSTIC: PEAK-DETECTION RATE ===\n\n');
fprintf(fid_pf,'Younger: %d/%d peak-detected (%.1f%%)  Older: %d/%d peak-detected (%.1f%%)\n', ...
    sum(peakDetected(idx_YA)), nYA, 100*mean(peakDetected(idx_YA)), ...
    sum(peakDetected(idx_OA)), nOA, 100*mean(peakDetected(idx_OA)));
fprintf(fid_pf,'Fisher exact test: p=%.4f\n\n', p_det);

[~, p_pf, ~, stats_pf] = ttest2(PF_final(idx_YA), PF_final(idx_OA));
d_pf = (mean(PF_final(idx_YA),'omitnan') - mean(PF_final(idx_OA),'omitnan')) / ...
       sqrt(((nYA-1)*var(PF_final(idx_YA),'omitnan') + (nOA-1)*var(PF_final(idx_OA),'omitnan')) / (nAll-2));
[p_lev, ~] = vartestn([PF_final(idx_YA); PF_final(idx_OA)], age_vec, 'TestType','LeveneAbsolute','Display','off');
[~, p_pf_w, ~, stats_pf_w] = ttest2(PF_final(idx_YA), PF_final(idx_OA), 'Vartype','unequal');

fprintf(fid_pf,'=== PEAK ALPHA FREQUENCY: YOUNGER VS OLDER (all subjects, peak+CoG mixed) ===\n\n');
fprintf(fid_pf,'Younger: M=%.3f SD=%.3f  Older: M=%.3f SD=%.3f\n', ...
    mean(PF_final(idx_YA),'omitnan'), std(PF_final(idx_YA),'omitnan'), ...
    mean(PF_final(idx_OA),'omitnan'), std(PF_final(idx_OA),'omitnan'));
fprintf(fid_pf,'Standard t-test: t(%d)=%.3f, p=%.4f, Cohen d=%.3f\n', stats_pf.df, stats_pf.tstat, p_pf, d_pf);
fprintf(fid_pf,'Levene''s test (variance homogeneity): p=%.4f\n', p_lev);
fprintf(fid_pf,'Welch''s t-test (unequal variances): t(%.1f)=%.3f, p=%.4f\n\n', stats_pf_w.df, stats_pf_w.tstat, p_pf_w);

idx_YA_peak = idx_YA(peakDetected(idx_YA));
idx_OA_peak = idx_OA(peakDetected(idx_OA));
nYA_p = length(idx_YA_peak); nOA_p = length(idx_OA_peak);
age_vec_peak = [ones(nYA_p,1); 2*ones(nOA_p,1)]; %#ok<NASGU>

[~, p_pf_c, ~, stats_pf_c] = ttest2(PF_peak(idx_YA_peak), PF_peak(idx_OA_peak));
d_pf_c = (mean(PF_peak(idx_YA_peak)) - mean(PF_peak(idx_OA_peak))) / ...
    sqrt(((nYA_p-1)*var(PF_peak(idx_YA_peak)) + (nOA_p-1)*var(PF_peak(idx_OA_peak))) / (nYA_p+nOA_p-2));
[p_lev_c, ~] = vartestn([PF_peak(idx_YA_peak); PF_peak(idx_OA_peak)], [ones(nYA_p,1);2*ones(nOA_p,1)], 'TestType','LeveneAbsolute','Display','off');
[~, p_pf_w_c, ~, stats_pf_w_c] = ttest2(PF_peak(idx_YA_peak), PF_peak(idx_OA_peak), 'Vartype','unequal');

fprintf(fid_pf,'=== PRIMARY ANALYSIS: PEAK-DETECTED SUBJECTS ONLY (n=%d YA, n=%d OA) ===\n\n', nYA_p, nOA_p);
fprintf(fid_pf,'Younger: M=%.3f SD=%.3f  Older: M=%.3f SD=%.3f\n', ...
    mean(PF_peak(idx_YA_peak)), std(PF_peak(idx_YA_peak)), mean(PF_peak(idx_OA_peak)), std(PF_peak(idx_OA_peak)));
fprintf(fid_pf,'Standard t-test: t(%d)=%.3f, p=%.4f, Cohen d=%.3f\n', stats_pf_c.df, stats_pf_c.tstat, p_pf_c, d_pf_c);
fprintf(fid_pf,'Levene''s test (variance homogeneity): p=%.4f\n', p_lev_c);
fprintf(fid_pf,'Welch''s t-test (unequal variances): t(%.1f)=%.3f, p=%.4f\n\n', stats_pf_w_c.df, stats_pf_w_c.tstat, p_pf_w_c);

[~, p_cog, ~, stats_cog] = ttest2(PF_CoG(idx_YA), PF_CoG(idx_OA));
fprintf(fid_pf,'=== ROBUSTNESS CHECK: CoG-ONLY (bypassing peak/CoG fallback) ===\n\n');
fprintf(fid_pf,'Younger: M=%.3f SD=%.3f  Older: M=%.3f SD=%.3f\n', ...
    mean(PF_CoG(idx_YA)), std(PF_CoG(idx_YA)), mean(PF_CoG(idx_OA)), std(PF_CoG(idx_OA)));
fprintf(fid_pf,'t(%d)=%.3f, p=%.4f\n\n', stats_cog.df, stats_cog.tstat, p_cog);

fclose(fid_pf);
fprintf('Stats saved to: %s\n', fullfile(pdat.out,'BEBO_PF_stats.txt'));

%% Plots
if ~exist(pdat.plotDir,'dir'), mkdir(pdat.plotDir); end

fig1 = figure('Units','centimeters','Position',[2 2 12 14],'Color','w','Visible','off');
ax = axes(fig1);
raincloud_compare_PF(PF_peak(idx_YA_peak), PF_peak(idx_OA_peak), col_YA, col_OA, ax);
title(ax,'Peak alpha frequency (peak-detected subjects only)'); ylabel(ax,'Frequency (Hz)');
set(ax,'XTick',[1 2],'XTickLabel',{'Younger','Older'});
savefig_all_PF(fig1, fullfile(pdat.plotDir,'FigPF_Raincloud'));

fig2 = figure('Units','centimeters','Position',[2 2 16 12],'Color','w','Visible','off');
ax = axes(fig2); hold(ax,'on');
m_YA = mean(spec_corrected(idx_YA,:), 1, 'omitnan');
se_YA = std(spec_corrected(idx_YA,:), 0, 1, 'omitnan') / sqrt(nYA);
m_OA = mean(spec_corrected(idx_OA,:), 1, 'omitnan');
se_OA = std(spec_corrected(idx_OA,:), 0, 1, 'omitnan') / sqrt(nOA);
yl = [0, max([m_YA+se_YA, m_OA+se_OA])*1.15];
patch(ax, [searchBand(1) searchBand(2) searchBand(2) searchBand(1)], ...
    [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9], 'EdgeColor','none','FaceAlpha',0.6, 'HandleVisibility','off');
fill_between_PF(ax, foiGrid, m_YA-se_YA, m_YA+se_YA, col_YA, 0.2);
fill_between_PF(ax, foiGrid, m_OA-se_OA, m_OA+se_OA, col_OA, 0.2);
h_YA = plot(ax, foiGrid, m_YA, '-', 'Color', col_YA, 'LineWidth', 1.5);
h_OA = plot(ax, foiGrid, m_OA, '-', 'Color', col_OA, 'LineWidth', 1.5);
% R2017b: xline (R2018b+) not available -- plain vertical lines instead
plot(ax, mean(PF_peak(idx_YA_peak))*[1 1], yl, '--', 'Color', col_YA, 'LineWidth', 1, 'HandleVisibility','off');
plot(ax, mean(PF_peak(idx_OA_peak))*[1 1], yl, '--', 'Color', col_OA, 'LineWidth', 1, 'HandleVisibility','off');
xlim(ax, [foiGrid(1) foiGrid(end)]); ylim(ax, yl); box(ax,'off');
xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Background-corrected power (a.u.)');
title(ax, sprintf('Group-average spectra (%.1f-%.1f Hz search band shaded)', searchBand(1), searchBand(2)));
legend([h_YA h_OA], {'Younger','Older'}, 'Location','best','Box','off','FontSize',8);
savefig_all_PF(fig2, fullfile(pdat.plotDir,'FigPF_GroupSpectra'));

fprintf('All PF figures saved to: %s\n', pdat.plotDir);

end


function [pks, locs, proms] = local_findpeaks(y, x, minProm)
% Replacement for findpeaks(y, x, 'MinPeakProminence', minProm) with
% findpeaks' definitions:
%  - peak: a sample larger than its left neighbour and not smaller than
%    its right one; a flat top counts once (its first sample) if the
%    signal drops after it. Endpoints are never peaks.
%  - prominence: height above the higher of the two lowest points
%    reached on either side before hitting a higher sample (or the
%    signal end).
y = y(:)'; x = x(:)';
n = length(y);
pks = []; locs = []; proms = [];
i = 2;
while i < n
    if y(i) > y(i-1)
        j = i;
        while j < n && y(j+1) == y(i), j = j + 1; end % flat top
        if j < n && y(j+1) < y(i)
            k = i - 1;
            while k >= 1 && y(k) <= y(i), k = k - 1; end
            leftMin = min(y(k+1:i));
            k = j + 1;
            while k <= n && y(k) <= y(i), k = k + 1; end
            rightMin = min(y(j:k-1));
            p = y(i) - max(leftMin, rightMin);
            if p >= minProm
                pks(end+1) = y(i); locs(end+1) = x(i); proms(end+1) = p; %#ok<AGROW>
            end
        end
        i = j + 1;
    else
        i = i + 1;
    end
end
end

function raincloud_compare_PF(data_YA, data_OA, col_YA, col_OA, ax)
hold(ax,'on');
raincloud_single_PF(data_YA(:), 1, col_YA, ax);
raincloud_single_PF(data_OA(:), 2, col_OA, ax);
xlim(ax,[0.3 2.7]); box(ax,'off');
end

function raincloud_single_PF(data, xpos, col, ax)
data = data(~isnan(data));
[f, xi] = ksdensity(data, 'NumPoints', 200);
f = f / max(f) * 0.35;
patch(ax, [xpos + f, xpos*ones(1,length(f))], [xi, fliplr(xi)], col, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
q1 = prctile(data,25); q3 = prctile(data,75); med = median(data);
wlo = max(min(data), q1-1.5*(q3-q1)); whi = min(max(data), q3+1.5*(q3-q1));
bw = 0.07;
patch(ax, xpos+[-bw bw bw -bw], [q1 q1 q3 q3], col, 'FaceAlpha', 0.85, 'EdgeColor', col*0.6, 'LineWidth', 0.8);
line(ax, xpos+[-bw bw], [med med], 'Color','w', 'LineWidth', 2);
line(ax, [xpos xpos], [wlo q1], 'Color', col*0.6, 'LineWidth', 1);
line(ax, [xpos xpos], [q3 whi], 'Color', col*0.6, 'LineWidth', 1);
jitter = (rand(length(data),1)-0.5)*0.15;
scatter(ax, xpos-0.22+jitter, data, 18, col, 'filled', 'MarkerFaceAlpha', 0.45, 'MarkerEdgeColor','none');
end

function fill_between_PF(ax, x, lo, hi, col, alpha)
patch(ax, [x, fliplr(x)], [lo, fliplr(hi)], col, 'FaceAlpha', alpha, 'EdgeColor','none', 'HandleVisibility','off');
end

function savefig_all_PF(fig, filepath)
% exportgraphics (R2020a+) is not available on the grid's MATLAB R2017b --
% print() instead, same fix as rsa_plots.m's savefig_all.
print(fig, [filepath '.pdf'], '-dpdf');
print(fig, [filepath '.png'], '-dpng', '-r300');
print(fig, [filepath '.tif'], '-dtiff', '-r300');
fprintf('Saved: %s (.pdf/.png/.tif)\n', filepath);
end

function ebosc_plots(pdat)
% eBOSC figure -- plotting only: loads Results_eBOSC/eBOSC_results.mat
% (computed by ebosc_results.m) and draws it. No computation (2026-10-05).
%   Row 1: spectrum above background (log10(power/background), all time
%          points, mean over the Ispec cluster channels): individual
%          subjects (thin) + group mean, YA | OA; bands shaded
%   Row 2: alpha rhythm: frequency | abundance | amplitude, YA vs OA
%   Row 3: beta rhythm:  frequency | abundance | amplitude, YA vs OA
%   Row 4: alpha frequency x Ispec | beta frequency x Ispec (regression
%          line per age group)
%   Row 5: group-mean peak frequency topographies: alpha YA | alpha OA |
%          beta YA | beta OA (* = Ispec cluster channels)

R = load(fullfile(pdat.res, 'eBOSC_results.mat'));
T = R.T; sp = R.spec; g = R.grp; B = R.bands;
iY = R.idx_YA; iO = R.idx_OA;
col_YA = [30 144 255]/255; col_OA = [80 80 80]/255;
bandCol = [0.92 0.92 0.92; 0.97 0.93 0.86];

fig = figure('Units','centimeters','Position',[2 2 28 44],'Color','w','Visible','off');
nR = 5;

dat = {sp.above_YA, sp.above_OA}; cols = {col_YA, col_OA}; ttl = {'Younger adults','Older adults'};
allv = [sp.above_YA(:); sp.above_OA(:)]; allv = allv(~isnan(allv));
yl = [min(allv) max(allv)] + [-0.05 0.05] * (max(allv) - min(allv));   % full range: no individual peak cut off
for gI = 1:2
    ax = subplot(nR,2,gI); hold(ax,'on');
    for b = 1:numel(B.names)
        patch(ax, [B.ranges(b,:) fliplr(B.ranges(b,:))], [yl(1) yl(1) yl(2) yl(2)], bandCol(b,:), 'EdgeColor','none','HandleVisibility','off');
    end
    plot(ax, sp.F, dat{gI}', '-', 'Color', [cols{gI} 0.25], 'LineWidth', 0.5);
    plot(ax, sp.F, mean(dat{gI},1,'omitnan'), '-', 'Color', cols{gI}*0.7, 'LineWidth', 2.5);
    plot(ax, [sp.F(1) sp.F(end)], [0 0], ':', 'Color', [0.4 0.4 0.4]);
    set(ax, 'XScale','log', 'XLim', [2 40], 'XTick', [2 4 8 14 25 40]); ylim(ax, yl); box(ax,'off');
    xlabel(ax,'Frequency (Hz)'); ylabel(ax,'log_{10}(power / background)');
    title(ax, {ttl{gI}, sprintf('individual spectra (n=%d); shaded: %s', size(dat{gI},1), strjoin(B.names, ', '))}, 'FontSize', 9);
end

rc = {'freq','abundance','amp'};
yl_ = {'Frequency (Hz)','Abundance','Amplitude above BG (log_{10})'};
for b = 1:numel(B.names)
    nm = B.names{b};
    for m = 1:3
        ax = subplot(nR,3,3*b + m);
        v = T.([nm '_' rc{m}]);
        raincloud_compare(v(iY), v(iO), col_YA, col_OA, ax);
        set(ax,'XTick',[1 2],'XTickLabel',{'Younger','Older'}); ylabel(ax, yl_{m});
        gg = g.(nm).(rc{m});
        title(ax, {sprintf('%s (%g-%g Hz)', nm, B.ranges(b,1), B.ranges(b,2)), ...
            sprintf('t(%.1f)=%.2f, p_{Holm}=%.3f', gg.df, gg.t, gg.p_Holm)}, 'FontSize', 9);
    end
end

for b = 1:numel(B.names)
    nm = B.names{b}; c = R.corrF.(nm);
    ax = subplot(nR,2,6+b); hold(ax,'on');
    x = T.([nm '_freq']);
    scatter(ax, x(iO), T.Ispec(iO), 30, col_OA, 'x', 'LineWidth', 1.2);
    scatter(ax, x(iY), T.Ispec(iY), 30, col_YA, 'o', 'filled', 'MarkerFaceAlpha', 0.6);
    regline(ax, x(iY), T.Ispec(iY), col_YA); regline(ax, x(iO), T.Ispec(iO), col_OA);
    xlabel(ax, sprintf('%s frequency (Hz)', nm)); ylabel(ax,'Item specificity'); box(ax,'off');
    title(ax, sprintf('%s frequency x Ispec: partial r=%.2f, p=%.3f', nm, c.r_part, c.p_part), 'FontSize', 9);
    legend(ax, {'Older','Younger'}, 'Location','best','Box','off','FontSize',8);
end

% Row 5: peak frequency topographies per band and group (shared colour
% scale within a band)
if exist('ft_topoplotER','file') == 2 && isfield(R, 'topo')
    pan = 0;
    for b = 1:numel(B.names)
        nm = B.names{b};
        v = [R.topo.(nm).PF_YA; R.topo.(nm).PF_OA]; v = v(~isnan(v));
        zl = [min(v) max(v)]; if diff(zl) == 0, zl = zl + [-0.5 0.5]; end
        for gI = 1:2
            pan = pan + 1;
            ax = subplot(nR,4,16+pan);
            if gI == 1, val = R.topo.(nm).PF_YA; gl = 'YA'; else, val = R.topo.(nm).PF_OA; gl = 'OA'; end
            dummy = struct('label', {R.topo.label}, 'avg', val, 'time', 0, 'dimord', 'chan_time');
            cfg = []; cfg.layout = pdat.layout; cfg.marker = 'off'; cfg.style = 'fill'; cfg.comment = 'no';
            cfg.zlim = zl; cfg.axes = ax; cfg.colormap = parula(256);
            cfg.highlight = 'on'; cfg.highlightchannel = R.clusterChans; cfg.highlightsymbol = '*';
            cfg.highlightsize = 5; cfg.highlightcolor = [0 0 0];
            ft_topoplotER(cfg, dummy);
            caxis(ax, zl); colormap(ax, parula(256)); cb = colorbar(ax); ylabel(cb, 'Hz');
            title(ax, sprintf('%s PF, %s', nm, gl), 'FontSize', 9);
        end
    end
end

if ~exist(pdat.fig,'dir'), mkdir(pdat.fig); end
fp = fullfile(pdat.fig, 'Fig_eBOSC');
set(fig, 'Units','centimeters'); pos = get(fig,'Position');
set(fig, 'PaperUnits','centimeters', 'PaperSize', pos(3:4), 'PaperPosition', [0 0 pos(3:4)]);
print(fig, [fp '.pdf'], '-dpdf');
print(fig, [fp '.png'], '-dpng', '-r300');
fprintf('Saved: %s (.pdf/.png)\n', fp);
end


function regline(ax, x, y, col)
ok = ~isnan(x) & ~isnan(y);
if sum(ok) < 3, return; end
p = polyfit(x(ok), y(ok), 1); xl = [min(x(ok)) max(x(ok))];
plot(ax, xl, polyval(p, xl), '-', 'Color', col, 'LineWidth', 1.2, 'HandleVisibility', 'off');
end

function raincloud_compare(data_YA, data_OA, col_YA, col_OA, ax)
hold(ax,'on');
raincloud_single(data_YA(:), 1, col_YA, ax);
raincloud_single(data_OA(:), 2, col_OA, ax);
xlim(ax,[0.3 2.7]); box(ax,'off');
end

function raincloud_single(data, xpos, col, ax)
data = data(~isnan(data));
if numel(data) < 3, return; end
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

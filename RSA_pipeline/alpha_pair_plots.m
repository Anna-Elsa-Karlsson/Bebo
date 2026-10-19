function alpha_pair_plots(pdat)
% PLOT step for alpha_pair_results.m -- loads alpha_pair_results.mat and
% draws, no computation (2026-10-04; control figure 2026-10-05).
%   Fig_AlphaPair          MAIN block (Isim vs Csim, alpha = mean of the
%                          pair's two trials)
%   Fig_AlphaPair_control  CONTROL block (Isim, Csim, Bsim; alpha = 1st
%                          pair member), only if it was computed
% Per figure: row 1 similarity per alpha quintile (within-person z) per
% type, YA | OA (mean +/- SE); row 2 per-person slopes, YA vs OA; row 3
% (if any) slope differences. Titles: one-sample t vs 0, all participants.
%   pdat : .out (where the results are), .fig

R = load(fullfile(pdat.out, 'alpha_pair_results.mat'));
if ~exist(pdat.fig, 'dir'), mkdir(pdat.fig); end
draw_block(R.blocks.main, R, 'Alpha power (mean of the pair''s two trials, z within person)', ...
    fullfile(pdat.fig, 'Fig_AlphaPair'));
if isfield(R.blocks, 'control')
    draw_block(R.blocks.control, R, 'Alpha power (1st pair member, z within person)', ...
        fullfile(pdat.fig, 'Fig_AlphaPair_control'));
end
end


function draw_block(B, R, xlab, filepath)
col_YA = [30 144 255]/255; col_OA = [80 80 80]/255;
colT = struct('Isim', [0 0 0], 'Csim', [120 70 215]/255, 'Bsim', [0.6 0.6 0.6]);
labT = struct('Isim', 'Item similarity', 'Csim', 'Category similarity', 'Bsim', 'Between-category (control)');
nT = numel(B.types);
isDiff = strncmp({B.meas.name}, 'diff_', 5);
rows = {find(~isDiff), find(isDiff)};
nRows = 1 + sum(~cellfun(@isempty, rows));

fig = figure('Units','centimeters','Position',[2 2 28 9.5*nRows],'Color','w','Visible','off');
grp = {'YA','OA'}; ttl = {'Younger adults','Older adults'};
for g = 1:2
    ax = subplot(nRows,2,g); hold(ax,'on');
    b = B.grpBin.(grp{g});
    for ty = 1:nT
        c = colT.(B.types{ty});
        errorbar(ax, b.pow_M(:,ty), b.sim_M(:,ty), b.sim_SE(:,ty), '-o', 'Color', c, 'MarkerFaceColor', c, 'LineWidth', 1.3);
    end
    xlabel(ax, {xlab, 'lower = stronger desync'}); ylabel(ax, 'Similarity (r)');
    title(ax, ttl{g}); box(ax,'off');
    legend(ax, cellfun(@(t) labT.(t), B.types, 'uni', 0), 'Location','best','Box','off','FontSize',8);
end
rowI = 1;
for rI = 1:2
    if isempty(rows{rI}), continue; end
    rowI = rowI + 1;
    for j = 1:numel(rows{rI})
        m = rows{rI}(j);
        ax = subplot(nRows,3,3*(rowI-1) + j); hold(ax,'on');
        v = B.meas(m).v;
        raincloud_single(v(R.idx_YA), 1, col_YA, ax);
        raincloud_single(v(R.idx_OA), 2, col_OA, ax);
        plot(ax, [0.3 2.7], [0 0], ':', 'Color', [0.5 0.5 0.5]);
        xlim(ax,[0.3 2.7]); set(ax,'XTick',[1 2],'XTickLabel',{'Younger','Older'}); box(ax,'off');
        g = B.G.(B.meas(m).name).all;
        title(ax, {strrep(B.meas(m).label, '->', '\rightarrow'), sprintf('all: t(%d)=%.2f, p=%.3f', g.df, g.t, g.p)}, 'FontSize', 9);
        if j == 1
            if rI == 1, ylabel(ax, 'Slope (similarity per SD alpha)'); else, ylabel(ax, 'Slope difference'); end
        end
    end
end
savefig_all(fig, filepath);
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

function savefig_all(fig, filepath)
set(fig, 'Units','centimeters'); pos = get(fig,'Position');
set(fig, 'PaperUnits','centimeters', 'PaperSize', pos(3:4), 'PaperPosition', [0 0 pos(3:4)]);
print(fig, [filepath '.pdf'], '-dpdf');
print(fig, [filepath '.png'], '-dpng', '-r300');
fprintf('Saved: %s (.pdf/.png)\n', filepath);
end

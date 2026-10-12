% Supplementary, robustness-oriented companion to plotCoupledMeshOrder.m.
%
% Motivation: the four coupled-mesh levels are four INDEPENDENT optimisation
% runs. At a given outer-iteration index they sit in different states, and
% because T_j = O(h_j) they also accumulate different pseudo-times. The
% per-iteration order estimator used by the main figure is therefore
% noise-dominated for the many iterations whose drift is near the floor.
%
% This script keeps the main figure untouched and adds:
%   (1) order from the ratio of per-level medians,
%   (2) order restricted to iterations matched at a COMMON accumulated
%       pseudo-time (exact integer matching when one exists),
%   (3) a diagnostic panel of median drift versus accumulated pseudo-time.
%
% Outputs:
%   Output/CoupledMeshRobustSummary.csv
%   figure/CoupledMeshOrderRobust.pdf / .png

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
% The study is selected through environment variables so that one script serves any
% coupled-mesh data set (e.g. the fixed-T study as well as the original one).
inputName = getenv('COUPLED_CSV');
if isempty(inputName)
    inputName = 'CoupledMeshStudy.csv';
end
tag = getenv('COUPLED_TAG');
if isempty(tag)
    tag = 'CoupledMesh';
end
dataFile = fullfile(projectRoot, 'Output', inputName);
summaryFile = fullfile(projectRoot, 'Output', [tag 'RobustSummary.csv']);
figureDir = fullfile(projectRoot, 'figure');

data = readtable(dataFile, 'VariableNamingRule', 'preserve');
if ~any(strcmp(data.Properties.VariableNames, 'transport_relative_drift'))
    data.transport_relative_drift = data.transport_abs_drift ./ data.volume_before;
end
if ~any(strcmp(data.Properties.VariableNames, 'step_relative_drift'))
    data.step_relative_drift = data.step_abs_drift ./ data.volume_before;
end

scales = unique(data.mesh_scale);
nLevels = numel(scales);
if nLevels < 2
    error('At least two mesh levels are required.');
end

levelH = zeros(nLevels, 1);
levelT = zeros(nLevels, 1);
for k = 1:nLevels
    rows = data(data.mesh_scale == scales(k), :);
    levelH(k) = rows.hmax(1);
    levelT(k) = rows.T(1);
end
[levelH, order] = sort(levelH, 'descend');
scales = scales(order);
levelT = levelT(order);

iterationCount = sum(data.mesh_scale == scales(1));
drift = nan(iterationCount, nLevels);
pseudoTime = nan(iterationCount, nLevels);
for k = 1:nLevels
    rows = sortrows(data(data.mesh_scale == scales(k), :), 'iteration');
    if height(rows) ~= iterationCount
        error('Coupled mesh levels have inconsistent row counts.');
    end
    drift(:, k) = rows.transport_relative_drift;
    pseudoTime(:, k) = rows.pseudo_time_after;
end

levelMedian = median(drift, 1, 'omitnan')';
maxPseudoTime = max(pseudoTime, [], 1)';

% ---- estimator A: median of per-iteration orders (what the main figure uses)
perIterOrder = nan(iterationCount, nLevels - 1);
for k = 1:nLevels - 1
    perIterOrder(:, k) = log(drift(:, k) ./ drift(:, k + 1)) / log(levelH(k) / levelH(k + 1));
end

% ---- estimator B: order from the ratio of per-level medians
ratioMedianOrder = zeros(nLevels - 1, 1);
for k = 1:nLevels - 1
    ratioMedianOrder(k) = log(levelMedian(k) / levelMedian(k + 1)) / log(levelH(k) / levelH(k + 1));
end

% ---- estimator C: matched at a common accumulated pseudo-time
% Find the largest g such that every T_j / g is an integer, then the smallest
% common pseudo-time step. With T = [.075 .06 .05 .0375] this yields
% g = .0025, multipliers [30 24 20 15], LCM 120, step .3, strides [4 5 6 8].
g = NaN;
for n = 1:100000
    candidate = min(levelT) / n;
    q = levelT / candidate;
    if all(abs(q - round(q)) < 1e-9)
        g = candidate;
        break;
    end
end
if isnan(g)
    g = min(levelT) / 100000;
end
multipliers = round(levelT / g);
L = multipliers(1);
for k = 2:numel(multipliers)
    L = lcm(L, multipliers(k));
end
pseudoTimeStep = L * g;
strides = L ./ multipliers;              % iterations per level for one common step
maxMatch = floor((iterationCount - 1) / max(strides));
matchOrder = nan(maxMatch, nLevels - 1);
matchedDrift = nan(maxMatch, nLevels);
if maxMatch >= 2
    for k = 1:nLevels
        idx = (1:maxMatch)' * strides(k);        % pseudo_time_before values
        matchedDrift(:, k) = drift(idx + 1, k);  % iteration index is zero-based
    end
    for k = 1:nLevels - 1
        matchOrder(:, k) = log(matchedDrift(:, k) ./ matchedDrift(:, k + 1)) / ...
            log(levelH(k) / levelH(k + 1));
    end
end

pair = (1:nLevels - 1)';
hCoarse = levelH(1:end - 1);
hFine = levelH(2:end);
orderRatioMedian = ratioMedianOrder;
orderPerIterMedian = zeros(nLevels - 1, 1);
orderPerIterQ25 = zeros(nLevels - 1, 1);
orderPerIterQ75 = zeros(nLevels - 1, 1);
orderMatchedMedian = nan(nLevels - 1, 1);
orderMatchedQ25 = nan(nLevels - 1, 1);
orderMatchedQ75 = nan(nLevels - 1, 1);
for k = 1:nLevels - 1
    values = perIterOrder(:, k);
    values = values(isfinite(values));
    orderPerIterMedian(k) = median(values);
    orderPerIterQ25(k) = prctile(values, 25);
    orderPerIterQ75(k) = prctile(values, 75);
    if maxMatch >= 2
        mv = matchOrder(:, k);
        mv = mv(isfinite(mv));
        orderMatchedMedian(k) = median(mv);
        orderMatchedQ25(k) = prctile(mv, 25);
        orderMatchedQ75(k) = prctile(mv, 75);
    end
end

summary = table(pair, hCoarse, hFine, orderRatioMedian, ...
    orderPerIterMedian, orderPerIterQ25, orderPerIterQ75, ...
    orderMatchedMedian, orderMatchedQ25, orderMatchedQ75, ...
    'VariableNames', {'pair', 'h_coarse', 'h_fine', ...
    'order_from_median_ratio', 'order_per_iteration_median', ...
    'order_per_iteration_q25', 'order_per_iteration_q75', ...
    'order_pseudotime_matched_median', 'order_pseudotime_matched_q25', ...
    'order_pseudotime_matched_q75'});
writetable(summary, summaryFile);

blue = [0.000, 0.447, 0.741];
red = [0.635, 0.078, 0.184];
green = [0.196, 0.553, 0.204];
gray = [0.25, 0.25, 0.25];

fig = figure('Units', 'centimeters', 'Position', [2, 2, 18, 8.8], 'Color', 'w');
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(layout);
errorbar(ax1, pair - 0.12, orderPerIterMedian, ...
    orderPerIterMedian - orderPerIterQ25, orderPerIterQ75 - orderPerIterMedian, ...
    'o', 'Color', blue, 'MarkerFaceColor', blue, 'LineWidth', 1.1, ...
    'MarkerSize', 5, 'CapSize', 5, 'DisplayName', 'Per-iteration median (IQR)');
hold(ax1, 'on');
plot(ax1, pair, orderRatioMedian, 's-', 'Color', red, 'MarkerFaceColor', red, ...
    'LineWidth', 1.2, 'MarkerSize', 5, 'DisplayName', 'Ratio of level medians');
if maxMatch >= 2
    errorbar(ax1, pair + 0.12, orderMatchedMedian, ...
        orderMatchedMedian - orderMatchedQ25, orderMatchedQ75 - orderMatchedMedian, ...
        'd', 'Color', green, 'MarkerFaceColor', green, 'LineWidth', 1.1, ...
        'MarkerSize', 5, 'CapSize', 5, 'DisplayName', 'Pseudo-time matched (IQR)');
end
yline(ax1, 2/3, '--', 'Color', gray, 'LineWidth', 1.15, 'DisplayName', 'Reference slope $2/3$');
yline(ax1, 1, ':', 'Color', [0.55, 0.55, 0.55], 'LineWidth', 1.0, 'DisplayName', 'First-order guide');
set(ax1, 'XTick', pair, 'TickLabelInterpreter', 'latex', ...
    'XTickLabel', arrayfun(@(k) sprintf('$h_%d\\to h_%d$', k, k + 1), pair, 'UniformOutput', false));
xlim(ax1, [0.55, nLevels - 0.45]);
% Reserve vertical headroom so the five-entry legend cannot cover the error bars,
% and keep the panel tag out of the legend's corner.
highs = [orderPerIterQ75; orderMatchedQ75];
lows = [orderPerIterQ25; orderMatchedQ25];
hi = max(highs(isfinite(highs)));
lo = min(lows(isfinite(lows)));
if ~isfinite(hi), hi = max(orderPerIterMedian); end
if ~isfinite(lo), lo = min(orderPerIterMedian); end
span = max(hi - lo, 1);
ylim(ax1, [lo - 0.15 * span, hi + 0.50 * span]);
grid(ax1, 'on');
xlabel(ax1, 'Adjacent mesh pair', 'Interpreter', 'latex');
ylabel(ax1, 'Coupled refinement order', 'Interpreter', 'latex');
legend(ax1, 'Location', 'northwest', 'Interpreter', 'latex');
text(ax1, 0.96, 0.95, '(a)', 'Units', 'normalized', 'FontWeight', 'bold', ...
    'FontSize', 10, 'HorizontalAlignment', 'right');

ax2 = nexttile(layout);
colors = lines(nLevels);
for k = 1:nLevels
    plot(ax2, pseudoTime(:, k), drift(:, k), '-', 'Color', colors(k, :), ...
        'LineWidth', 1.0, 'DisplayName', sprintf('$n=%d$  ($T=%.4g$)', scales(k), levelT(k)));
    hold(ax2, 'on');
end
set(ax2, 'YScale', 'log');
grid(ax2, 'on');
xlabel(ax2, 'Accumulated pseudo-time', 'Interpreter', 'latex');
ylabel(ax2, 'Relative transport drift', 'Interpreter', 'latex');
% The curves decay away from the origin, so the lower-left corner stays clear.
legend(ax2, 'Location', 'southwest', 'Interpreter', 'latex');
text(ax2, 0.035, 0.94, '(b)', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 10);

set([ax1, ax2], 'FontName', 'Times New Roman', 'FontSize', 8.5, ...
    'LineWidth', 0.8, 'Box', 'on');
% Suppress interactive toolbars, otherwise exportgraphics embeds an axes toolbar
% overlay in the vector output.
set(fig, 'ToolBar', 'none');
set([ax1, ax2], 'Toolbar', []);
drawnow;
exportgraphics(fig, fullfile(figureDir, [tag 'OrderRobust.pdf']), 'ContentType', 'vector');
exportgraphics(fig, fullfile(figureDir, [tag 'OrderRobust.png']), 'Resolution', 300);

fprintf('levels=%d  common pseudo-time step=%g  strides=[%s]  matched points=%d\n', ...
    nLevels, pseudoTimeStep, num2str(strides), maxMatch);
fprintf('max accumulated pseudo-time per level: [%s]\n', num2str(maxPseudoTime', '%.4g '));
fprintf('Saved robust summary: %s\n', summaryFile);
fprintf('Saved robust figure in: %s\n', figureDir);

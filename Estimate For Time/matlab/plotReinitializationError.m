scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
baseFile = fullfile(projectRoot, 'Output', 'StepDiagnostics.csv');
levelSummaryFile = fullfile(projectRoot, 'Output', ...
    'ReinitializationErrorSummary.csv');
orderSummaryFile = fullfile(projectRoot, 'Output', ...
    'ReinitializationOrderSummary.csv');
figureDir = fullfile(projectRoot, 'figure');

base = readtable(baseFile, 'VariableNamingRule', 'preserve');
baseReinit = base(logical(base.did_reinitialize), :);
baseRelativeJump = baseReinit.reinit_abs_jump ./ baseReinit.volume_before;

scales = [2; 3; 5; 7];
nLevels = numel(scales);
levelH = zeros(nLevels, 1);
levelMedian = zeros(nLevels, 1);
levelQ25 = zeros(nLevels, 1);
levelQ75 = zeros(nLevels, 1);
reinitMatrix = [];

for k = 1:nLevels
    fileName = fullfile(projectRoot, 'Output', ...
        sprintf('CoupledMeshFixedT_n%d.csv', scales(k)));
    rows = readtable(fileName, 'VariableNamingRule', 'preserve');
    rows = sortrows(rows, 'iteration');
    reinitRows = rows(logical(rows.did_reinitialize), :);
    values = reinitRows.reinit_abs_jump ./ reinitRows.volume_before;
    if k == 1
        reinitMatrix = zeros(numel(values), nLevels);
        reinitIterations = reinitRows.iteration;
    elseif numel(values) ~= size(reinitMatrix, 1) || ...
            any(reinitRows.iteration ~= reinitIterations)
        error('Reinitialization schedules differ between mesh levels.');
    end
    reinitMatrix(:, k) = values;
    levelH(k) = rows.hmax(1);
    levelMedian(k) = median(values);
    levelQ25(k) = prctile(values, 25);
    levelQ75(k) = prctile(values, 75);
end

[levelH, order] = sort(levelH, 'descend');
scales = scales(order);
levelMedian = levelMedian(order);
levelQ25 = levelQ25(order);
levelQ75 = levelQ75(order);
reinitMatrix = reinitMatrix(:, order);

nPairs = nLevels - 1;
pairMedian = zeros(nPairs, 1);
pairQ25 = zeros(nPairs, 1);
pairQ75 = zeros(nPairs, 1);
orderFromMedianRatio = zeros(nPairs, 1);
for k = 1:nPairs
    denominator = log(levelH(k) / levelH(k + 1));
    values = log(reinitMatrix(:, k) ./ reinitMatrix(:, k + 1)) / denominator;
    values = values(isfinite(values));
    pairMedian(k) = median(values);
    pairQ25(k) = prctile(values, 25);
    pairQ75(k) = prctile(values, 75);
    orderFromMedianRatio(k) = log(levelMedian(k) / levelMedian(k + 1)) / denominator;
end

levelSummary = table((1:nLevels)', scales, levelH, levelMedian, levelQ25, ...
    levelQ75, 'VariableNames', {'level', 'mesh_scale', 'hmax', ...
    'median_relative_reinit_jump', 'q25_relative_reinit_jump', ...
    'q75_relative_reinit_jump'});
pairSummary = table((1:nPairs)', levelH(1:end-1), levelH(2:end), ...
    orderFromMedianRatio, pairMedian, pairQ25, pairQ75, ...
    'VariableNames', {'pair', 'h_coarse', 'h_fine', ...
    'order_from_median_ratio', 'per_event_order_median', ...
    'per_event_order_q25', 'per_event_order_q75'});
writetable(levelSummary, levelSummaryFile);
writetable(pairSummary, orderSummaryFile);

blue = [0.000, 0.447, 0.741];
red = [0.635, 0.078, 0.184];
green = [0.196, 0.553, 0.204];
gray = [0.25, 0.25, 0.25];

fig = figure('Units', 'centimeters', 'Position', [2, 2, 18, 14], ...
    'Color', 'w');
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax1 = nexttile(layout, [1, 2]);
semilogy(ax1, baseReinit.iteration, baseRelativeJump, 'o-', ...
    'Color', red, 'MarkerFaceColor', red, 'LineWidth', 1.0, ...
    'MarkerSize', 4, 'DisplayName', 'Reinitialization jump');
hold(ax1, 'on');
yline(ax1, median(baseRelativeJump), '--', 'Color', gray, ...
    'LineWidth', 1.1, 'DisplayName', 'Median');
grid(ax1, 'on');
xlabel(ax1, 'Outer iteration', 'Interpreter', 'latex');
ylabel(ax1, 'Relative reinitialization jump', 'Interpreter', 'latex');
legend(ax1, 'Location', 'best', 'Interpreter', 'latex');
text(ax1, 0.018, 0.91, '(a)', 'Units', 'normalized', ...
    'FontWeight', 'bold', 'FontSize', 10);

ax2 = nexttile(layout);
errorbar(ax2, levelH, levelMedian, levelMedian-levelQ25, ...
    levelQ75-levelMedian, 'o-', 'Color', blue, ...
    'MarkerFaceColor', blue, 'LineWidth', 1.2, 'MarkerSize', 5, ...
    'CapSize', 6, 'DisplayName', 'Median (IQR)');
set(ax2, 'XScale', 'log', 'YScale', 'log');
grid(ax2, 'on');
xlabel(ax2, 'Measured mesh size $h_{\max}$', 'Interpreter', 'latex');
ylabel(ax2, 'Relative reinitialization jump', 'Interpreter', 'latex');
legend(ax2, 'Location', 'best', 'Interpreter', 'latex');
text(ax2, 0.035, 0.91, '(b)', 'Units', 'normalized', ...
    'FontWeight', 'bold', 'FontSize', 10);

ax3 = nexttile(layout);
pairs = (1:nPairs)';
errorbar(ax3, pairs-0.06, pairMedian, pairMedian-pairQ25, ...
    pairQ75-pairMedian, 'd', 'Color', green, ...
    'MarkerFaceColor', green, 'LineWidth', 1.1, 'MarkerSize', 5, ...
    'CapSize', 5, 'DisplayName', 'Per-event median (IQR)');
hold(ax3, 'on');
plot(ax3, pairs+0.06, orderFromMedianRatio, 's-', 'Color', red, ...
    'MarkerFaceColor', red, 'LineWidth', 1.1, 'MarkerSize', 5, ...
    'DisplayName', 'Ratio of level medians');
yline(ax3, 0, '-', 'Color', [0.65, 0.65, 0.65], ...
    'HandleVisibility', 'off');
yline(ax3, 2/3, '--', 'Color', gray, 'LineWidth', 1.0, ...
    'DisplayName', 'Reference slope $2/3$');
set(ax3, 'XTick', pairs, ...
    'XTickLabel', {'$h_1\to h_2$', '$h_2\to h_3$', '$h_3\to h_4$'}, ...
    'TickLabelInterpreter', 'latex');
xlim(ax3, [0.6, 3.4]);
grid(ax3, 'on');
xlabel(ax3, 'Adjacent mesh pair', 'Interpreter', 'latex');
ylabel(ax3, 'Reinitialization-jump order', 'Interpreter', 'latex');
legend(ax3, 'Location', 'best', 'Interpreter', 'latex');
text(ax3, 0.035, 0.91, '(c)', 'Units', 'normalized', ...
    'FontWeight', 'bold', 'FontSize', 10);

set([ax1, ax2, ax3], 'FontName', 'Times New Roman', 'FontSize', 8.5, ...
    'LineWidth', 0.8, 'Box', 'on');
set(fig, 'ToolBar', 'none');
set([ax1, ax2, ax3], 'Toolbar', []);
drawnow;
exportgraphics(fig, fullfile(figureDir, 'ReinitializationError.pdf'), ...
    'ContentType', 'vector');
exportgraphics(fig, fullfile(figureDir, 'ReinitializationError.png'), ...
    'Resolution', 300);

fprintf('Base-grid reinitialization events: %d\n', height(baseReinit));
fprintf('Base-grid median relative jump: %.8e\n', median(baseRelativeJump));
fprintf('Saved level summary: %s\n', levelSummaryFile);
fprintf('Saved order summary: %s\n', orderSummaryFile);
fprintf('Saved reinitialization figure in: %s\n', figureDir);

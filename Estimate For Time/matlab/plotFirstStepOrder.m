% First-step (iteration 0) volume-drift analysis for the coupled-mesh study.
%
% Why only the first step: the mesh levels are independent optimisation runs. From
% iteration 1 onward they sit in different states, so iteration-matched comparisons
% mix discretisation error with trajectory divergence. The very first step is the only
% one taken from a common state -- the same analytic initial shape sampled on each
% mesh -- so it is the defensible basis for an order estimate.
%
% Input : Output/CoupledMeshFirstStep.csv  (one row per mesh level)
% Output: Output/FirstStepOrderSummary.csv
%         figure/FirstStepOrder.pdf / .png

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
% The study is selected through environment variables so that one script serves both
% the T = O(h) study and the fixed-T study.
inputName = getenv('FIRSTSTEP_CSV');
if isempty(inputName)
    inputName = 'CoupledMeshFirstStep.csv';
end
tag = getenv('FIRSTSTEP_TAG');
if isempty(tag)
    tag = 'FirstStep';
end
dataFile = fullfile(projectRoot, 'Output', inputName);
summaryFile = fullfile(projectRoot, 'Output', [tag 'OrderSummary.csv']);
figureDir = fullfile(projectRoot, 'figure');

data = readtable(dataFile, 'VariableNamingRule', 'preserve');
data = sortrows(data, 'hmax', 'descend');

h = data.hmax;
T = data.T;
E = data.transport_relative_drift;
nLevels = height(data);
if nLevels < 2
    error('At least two mesh levels are required.');
end

pair = (1:nLevels - 1)';
hCoarse = h(1:end - 1);
hFine = h(2:end);
tCoarse = T(1:end - 1);
tFine = T(2:end);
eCoarse = E(1:end - 1);
eFine = E(2:end);
orderMeasured = log(eCoarse ./ eFine) ./ log(hCoarse ./ hFine);
orderIfProportionalToT = log(tCoarse ./ tFine) ./ log(hCoarse ./ hFine);
tIsConstant = all(abs(T / T(1) - 1) < 1e-12);

summary = table(pair, hCoarse, hFine, tCoarse, tFine, eCoarse, eFine, ...
    orderMeasured, orderIfProportionalToT, ...
    'VariableNames', {'pair', 'h_coarse', 'h_fine', 'T_coarse', 'T_fine', ...
    'E_coarse', 'E_fine', 'order_measured', 'order_if_E_proportional_to_T'});
writetable(summary, summaryFile);

blue = [0.000, 0.447, 0.741];
red = [0.635, 0.078, 0.184];
gray = [0.25, 0.25, 0.25];

fig = figure('Units', 'centimeters', 'Position', [2, 2, 18, 8.8], 'Color', 'w');
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% ---- (a) first-step drift versus mesh size, with reference behaviours
ax1 = nexttile(layout);
hRef = logspace(log10(min(h)) - 0.03, log10(max(h)) + 0.03, 50);
plot(ax1, hRef, E(1) * (hRef / h(1)).^(2/3), '--', 'Color', gray, ...
    'LineWidth', 1.15, 'DisplayName', 'Slope $2/3$ reference');
hold(ax1, 'on');
if ~tIsConstant
    plot(ax1, h, E(1) * (T / T(1)), ':s', 'Color', red, 'LineWidth', 1.15, ...
        'MarkerSize', 4.5, 'MarkerFaceColor', red, ...
        'DisplayName', 'Expected if $E \propto T$');
end
loglog(ax1, h, E, 'o-', 'Color', blue, 'MarkerFaceColor', blue, ...
    'LineWidth', 1.3, 'MarkerSize', 5.5, 'DisplayName', 'Measured first-step drift');
set(ax1, 'XScale', 'log', 'YScale', 'log');
% The curves sit in the upper half, so the lower-right corner stays clear of them.
ylim(ax1, [min(E) * 0.45, max(E) * 3]);
grid(ax1, 'on');
xlabel(ax1, 'Measured mesh size $h_{\max}$', 'Interpreter', 'latex');
ylabel(ax1, 'First-step relative volume drift', 'Interpreter', 'latex');
legend(ax1, 'Location', 'southeast', 'Interpreter', 'latex');
text(ax1, 0.035, 0.94, '(a)', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 10);

% ---- (b) adjacent-pair order actually observed
ax2 = nexttile(layout);
b1 = bar(ax2, pair - 0.17, orderMeasured, 0.30, 'FaceColor', blue, 'EdgeColor', 'none');
hold(ax2, 'on');
barHandles = b1;
barLabels = {'Measured order'};
if ~tIsConstant
    b2 = bar(ax2, pair + 0.17, orderIfProportionalToT, 0.30, ...
        'FaceColor', red, 'EdgeColor', 'none');
    barHandles = [b1, b2];
    barLabels = {'Measured order', 'Expected if $E \propto T$'};
end
refLine = yline(ax2, 2/3, '--', 'Color', gray, 'LineWidth', 1.15);
orderLine = yline(ax2, 1, ':', 'Color', [0.55, 0.55, 0.55], 'LineWidth', 1.0);
set(ax2, 'XTick', pair, 'TickLabelInterpreter', 'latex', ...
    'XTickLabel', arrayfun(@(k) sprintf('$h_%d\\to h_%d$', k, k + 1), pair, 'UniformOutput', false));
xlim(ax2, [0.4, nLevels - 0.6]);
% Headroom plus a top-right legend so it cannot crowd the tall middle bars.
ylim(ax2, [0, 2.9]);
grid(ax2, 'on');
xlabel(ax2, 'Adjacent mesh pair', 'Interpreter', 'latex');
ylabel(ax2, 'First-step scaling order', 'Interpreter', 'latex');
legend(ax2, [barHandles, refLine, orderLine], ...
    [barLabels, {'Reference slope $2/3$', 'First-order guide'}], ...
    'Location', 'northeast', 'Interpreter', 'latex');
text(ax2, 0.035, 0.94, '(b)', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 10);

set([ax1, ax2], 'FontName', 'Times New Roman', 'FontSize', 8.5, ...
    'LineWidth', 0.8, 'Box', 'on');
% Suppress interactive toolbars, otherwise exportgraphics embeds an axes toolbar
% overlay in the vector output.
set(fig, 'ToolBar', 'none');
set([ax1, ax2], 'Toolbar', []);
drawnow;
exportgraphics(fig, fullfile(figureDir, [tag 'Order.pdf']), 'ContentType', 'vector');
exportgraphics(fig, fullfile(figureDir, [tag 'Order.png']), 'Resolution', 300);

fprintf('T constant across levels: %d\n', tIsConstant);
fprintf('Saved first-step order summary: %s\n', summaryFile);
fprintf('Saved first-step order figure in: %s\n', figureDir);
disp(summary);

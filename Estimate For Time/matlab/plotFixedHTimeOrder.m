scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
stepFile = fullfile(projectRoot, 'Output', 'StepDiagnostics.csv');
orderFile = fullfile(projectRoot, 'Output', 'StepOrder.csv');
summaryFile = fullfile(projectRoot, 'Output', 'FixedHTimeStudy.csv');
figureDir = fullfile(projectRoot, 'figure');

step = readtable(stepFile, 'VariableNamingRule', 'preserve');
order = readtable(orderFile, 'VariableNamingRule', 'preserve');
if max(order.level) < 3
    error('Four levels T, T/2, T/4, T/8 are required.');
end

V0 = step.volume_before(1);
levels = (0:3)';
duration = zeros(4, 1);
medianTransportError = zeros(4, 1);
q25TransportError = zeros(4, 1);
q75TransportError = zeros(4, 1);
medianTransportOrder = nan(4, 1);
q25TransportOrder = nan(4, 1);
q75TransportOrder = nan(4, 1);
medianReinitOrder = nan(4, 1);
q25ReinitOrder = nan(4, 1);
q75ReinitOrder = nan(4, 1);

for k = 1:4
    rows = order(order.level == levels(k), :);
    duration(k) = rows.T(1);
    errors = rows.transport_abs_drift / V0;
    medianTransportError(k) = median(errors);
    q25TransportError(k) = prctile(errors, 25);
    q75TransportError(k) = prctile(errors, 75);
    if levels(k) > 0
        rates = rows.transport_observed_rate(isfinite(rows.transport_observed_rate));
        reinitRates = rows.step_observed_rate(logical(rows.did_reinitialize));
        reinitRates = reinitRates(isfinite(reinitRates));
        medianTransportOrder(k) = median(rates);
        q25TransportOrder(k) = prctile(rates, 25);
        q75TransportOrder(k) = prctile(rates, 75);
        medianReinitOrder(k) = median(reinitRates);
        q25ReinitOrder(k) = prctile(reinitRates, 25);
        q75ReinitOrder(k) = prctile(reinitRates, 75);
    end
end

summary = table(levels, duration, medianTransportError, ...
    q25TransportError, q75TransportError, medianTransportOrder, ...
    q25TransportOrder, q75TransportOrder, medianReinitOrder, ...
    q25ReinitOrder, q75ReinitOrder);
writetable(summary, summaryFile);

blue = [0.000, 0.447, 0.741];
red = [0.635, 0.078, 0.184];
gray = [0.45, 0.45, 0.45];
fig = figure('Units', 'centimeters', 'Position', [2, 2, 18, 8.8], ...
    'Color', 'w');
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(layout);
x = (1:4)';
errorbar(ax1, x, medianTransportError, ...
    medianTransportError-q25TransportError, ...
    q75TransportError-medianTransportError, 'o-', ...
    'Color', blue, 'MarkerFaceColor', blue, 'LineWidth', 1.25, ...
    'MarkerSize', 5, 'CapSize', 6);
set(ax1, 'YScale', 'log', 'XTick', x, ...
    'XTickLabel', {'$T$', '$T/2$', '$T/4$', '$T/8$'}, ...
    'TickLabelInterpreter', 'latex');
xlim(ax1, [0.65, 4.35]);
grid(ax1, 'on');
xlabel(ax1, 'Single-step transport interval', 'Interpreter', 'latex');
ylabel(ax1, 'Median relative volume drift (IQR)', 'Interpreter', 'latex');
text(ax1, 0.035, 0.94, '(a)', 'Units', 'normalized', ...
    'FontWeight', 'bold', 'FontSize', 10);

ax2 = nexttile(layout);
pairs = (1:3)';
errorbar(ax2, pairs-0.07, medianTransportOrder(2:4), ...
    medianTransportOrder(2:4)-q25TransportOrder(2:4), ...
    q75TransportOrder(2:4)-medianTransportOrder(2:4), 'o', ...
    'Color', blue, 'MarkerFaceColor', blue, 'LineWidth', 1.2, ...
    'MarkerSize', 5, 'CapSize', 6, ...
    'DisplayName', 'Transport drift order');
hold(ax2, 'on');
errorbar(ax2, pairs+0.07, medianReinitOrder(2:4), ...
    medianReinitOrder(2:4)-q25ReinitOrder(2:4), ...
    q75ReinitOrder(2:4)-medianReinitOrder(2:4), 's', ...
    'Color', red, 'MarkerFaceColor', red, 'LineWidth', 1.1, ...
    'MarkerSize', 5, 'CapSize', 6, ...
    'DisplayName', 'Reinitialization-step order');
yline(ax2, 1/3, '--', 'Color', [0.1, 0.1, 0.1], 'LineWidth', 1.1, ...
    'DisplayName', 'Time-dominated bound rate $1/3$');
yline(ax2, 1, ':', 'Color', gray, 'LineWidth', 1.1, ...
    'DisplayName', 'First-order guide');
set(ax2, 'XTick', pairs, ...
    'XTickLabel', {'$T\to T/2$', '$T/2\to T/4$', '$T/4\to T/8$'}, ...
    'TickLabelInterpreter', 'latex');
xlim(ax2, [0.6, 3.4]);
grid(ax2, 'on');
xlabel(ax2, 'Adjacent interval pair', 'Interpreter', 'latex');
ylabel(ax2, 'Single-step scaling order', 'Interpreter', 'latex');
legend(ax2, 'Location', 'best', 'Interpreter', 'latex');
text(ax2, 0.035, 0.94, '(b)', 'Units', 'normalized', ...
    'FontWeight', 'bold', 'FontSize', 10);

set([ax1, ax2], 'FontName', 'Times New Roman', 'FontSize', 8.5, ...
    'LineWidth', 0.8, 'Box', 'on');
exportgraphics(fig, fullfile(figureDir, 'FixedHTimeOrder.pdf'), ...
    'ContentType', 'vector');
exportgraphics(fig, fullfile(figureDir, 'FixedHTimeOrder.png'), ...
    'Resolution', 300);

fprintf('Saved fixed-h summary: %s\n', summaryFile);
fprintf('Saved fixed-h order figure in: %s\n', figureDir);

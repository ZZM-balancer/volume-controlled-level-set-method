scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
dataFile = fullfile(projectRoot, 'Output', 'CoupledMeshStudy.csv');
levelSummaryFile = fullfile(projectRoot, 'Output', 'CoupledMeshLevelSummary.csv');
orderSummaryFile = fullfile(projectRoot, 'Output', 'CoupledMeshOrderSummary.csv');
figureDir = fullfile(projectRoot, 'figure');

data = readtable(dataFile, 'VariableNamingRule', 'preserve');
scales = unique(data.mesh_scale);
if numel(scales) ~= 4
    error('Exactly four mesh levels are required.');
end

levelH = zeros(4, 1);
for k = 1:4
    rows = data(data.mesh_scale == scales(k), :);
    levelH(k) = rows.hmax(1);
end
[levelH, sortIndex] = sort(levelH, 'descend');
scales = scales(sortIndex);

iterationCount = sum(data.mesh_scale == scales(1));
if iterationCount ~= 300
    error('Each coupled mesh level must contain 300 outer iterations.');
end

transportError = zeros(iterationCount, 4);
stepError = zeros(iterationCount, 4);
isReinit = false(iterationCount, 1);
levelT = zeros(4, 1);
levelGamma = zeros(4, 1);
for k = 1:4
    rows = sortrows(data(data.mesh_scale == scales(k), :), 'iteration');
    if height(rows) ~= iterationCount
        error('Coupled mesh levels have inconsistent row counts.');
    end
    transportError(:, k) = rows.transport_relative_drift;
    stepError(:, k) = rows.step_relative_drift;
    levelT(k) = rows.T(1);
    levelGamma(k) = rows.gamma(1);
    if k == 1
        isReinit = logical(rows.did_reinitialize);
    elseif any(logical(rows.did_reinitialize) ~= isReinit)
        error('Reinitialization schedules differ between mesh levels.');
    end
end

medianTransportError = median(transportError, 1)';
q25TransportError = prctile(transportError, 25, 1)';
q75TransportError = prctile(transportError, 75, 1)';
medianReinitStepError = median(stepError(isReinit, :), 1)';
q25ReinitStepError = prctile(stepError(isReinit, :), 25, 1)';
q75ReinitStepError = prctile(stepError(isReinit, :), 75, 1)';

transportOrder = zeros(iterationCount, 3);
reinitOrder = zeros(sum(isReinit), 3);
for k = 1:3
    denominator = log(levelH(k) / levelH(k+1));
    transportOrder(:, k) = log(transportError(:, k) ./ ...
        transportError(:, k+1)) / denominator;
    reinitOrder(:, k) = log(stepError(isReinit, k) ./ ...
        stepError(isReinit, k+1)) / denominator;
end

medianTransportOrder = zeros(3, 1);
q25TransportOrder = zeros(3, 1);
q75TransportOrder = zeros(3, 1);
medianReinitOrder = zeros(3, 1);
q25ReinitOrder = zeros(3, 1);
q75ReinitOrder = zeros(3, 1);
for k = 1:3
    values = transportOrder(:, k);
    values = values(isfinite(values));
    medianTransportOrder(k) = median(values);
    q25TransportOrder(k) = prctile(values, 25);
    q75TransportOrder(k) = prctile(values, 75);
    values = reinitOrder(:, k);
    values = values(isfinite(values));
    medianReinitOrder(k) = median(values);
    q25ReinitOrder(k) = prctile(values, 25);
    q75ReinitOrder(k) = prctile(values, 75);
end

levelSummary = table((1:4)', scales, levelH, levelT, levelGamma, ...
    levelT./levelH, medianTransportError, q25TransportError, ...
    q75TransportError, medianReinitStepError, q25ReinitStepError, ...
    q75ReinitStepError, 'VariableNames', {'level', 'mesh_scale', 'hmax', ...
    'T', 'gamma', 'T_over_hmax', 'median_transport_error', ...
    'q25_transport_error', 'q75_transport_error', ...
    'median_reinit_step_error', 'q25_reinit_step_error', ...
    'q75_reinit_step_error'});
pairSummary = table((1:3)', levelH(1:3), levelH(2:4), ...
    medianTransportOrder, q25TransportOrder, q75TransportOrder, ...
    medianReinitOrder, q25ReinitOrder, q75ReinitOrder, ...
    'VariableNames', {'pair', 'h_coarse', 'h_fine', ...
    'median_transport_order', 'q25_transport_order', ...
    'q75_transport_order', 'median_reinit_order', ...
    'q25_reinit_order', 'q75_reinit_order'});
writetable(levelSummary, levelSummaryFile);
writetable(pairSummary, orderSummaryFile);

blue = [0.000, 0.447, 0.741];
red = [0.635, 0.078, 0.184];
gray = [0.25, 0.25, 0.25];
reference = medianTransportError(1) * (levelH / levelH(1)).^(2/3);

fig = figure('Units', 'centimeters', 'Position', [2, 2, 18, 8.8], ...
    'Color', 'w');
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(layout);
errorbar(ax1, levelH, medianTransportError, ...
    medianTransportError-q25TransportError, ...
    q75TransportError-medianTransportError, 'o-', ...
    'Color', blue, 'MarkerFaceColor', blue, 'LineWidth', 1.25, ...
    'MarkerSize', 5, 'CapSize', 6, 'DisplayName', 'Transport drift');
hold(ax1, 'on');
errorbar(ax1, levelH, medianReinitStepError, ...
    medianReinitStepError-q25ReinitStepError, ...
    q75ReinitStepError-medianReinitStepError, 's--', ...
    'Color', red, 'MarkerFaceColor', red, 'LineWidth', 1.1, ...
    'MarkerSize', 5, 'CapSize', 6, ...
    'DisplayName', 'Reinitialization-step drift');
plot(ax1, levelH, reference, '--', 'Color', gray, 'LineWidth', 1.15, ...
    'DisplayName', 'Reference slope $2/3$');
set(ax1, 'XScale', 'log', 'YScale', 'log');
grid(ax1, 'on');
xlabel(ax1, 'Measured mesh size $h_{\max}$', 'Interpreter', 'latex');
ylabel(ax1, 'Median relative one-step drift (IQR)', 'Interpreter', 'latex');
legend(ax1, 'Location', 'best', 'Interpreter', 'latex');
text(ax1, 0.035, 0.94, '(a)', 'Units', 'normalized', ...
    'FontWeight', 'bold', 'FontSize', 10);

ax2 = nexttile(layout);
pairs = (1:3)';
errorbar(ax2, pairs-0.07, medianTransportOrder, ...
    medianTransportOrder-q25TransportOrder, ...
    q75TransportOrder-medianTransportOrder, 'o', ...
    'Color', blue, 'MarkerFaceColor', blue, 'LineWidth', 1.2, ...
    'MarkerSize', 5, 'CapSize', 6, 'DisplayName', 'Transport drift order');
hold(ax2, 'on');
errorbar(ax2, pairs+0.07, medianReinitOrder, ...
    medianReinitOrder-q25ReinitOrder, ...
    q75ReinitOrder-medianReinitOrder, 's', ...
    'Color', red, 'MarkerFaceColor', red, 'LineWidth', 1.1, ...
    'MarkerSize', 5, 'CapSize', 6, ...
    'DisplayName', 'Reinitialization-step order');
yline(ax2, 2/3, '--', 'Color', gray, 'LineWidth', 1.15, ...
    'DisplayName', 'Reference slope $2/3$');
yline(ax2, 1, ':', 'Color', [0.55, 0.55, 0.55], 'LineWidth', 1.0, ...
    'DisplayName', 'First-order guide');
set(ax2, 'XTick', pairs, ...
    'XTickLabel', {'$h_1\to h_2$', '$h_2\to h_3$', '$h_3\to h_4$'}, ...
    'TickLabelInterpreter', 'latex');
xlim(ax2, [0.6, 3.4]);
grid(ax2, 'on');
xlabel(ax2, 'Adjacent mesh pair', 'Interpreter', 'latex');
ylabel(ax2, 'Coupled refinement order', 'Interpreter', 'latex');
legend(ax2, 'Location', 'best', 'Interpreter', 'latex');
text(ax2, 0.035, 0.94, '(b)', 'Units', 'normalized', ...
    'FontWeight', 'bold', 'FontSize', 10);

set([ax1, ax2], 'FontName', 'Times New Roman', 'FontSize', 8.5, ...
    'LineWidth', 0.8, 'Box', 'on');
exportgraphics(fig, fullfile(figureDir, 'CoupledMeshOrder.pdf'), ...
    'ContentType', 'vector');
exportgraphics(fig, fullfile(figureDir, 'CoupledMeshOrder.png'), ...
    'Resolution', 300);

fprintf('Saved coupled mesh level summary: %s\n', levelSummaryFile);
fprintf('Saved coupled mesh order summary: %s\n', orderSummaryFile);
fprintf('Saved coupled mesh order figure in: %s\n', figureDir);

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
stepFile = fullfile(projectRoot, 'Output', 'StepDiagnostics.csv');
orderFile = fullfile(projectRoot, 'Output', 'StepOrder.csv');
summaryFile = fullfile(projectRoot, 'Output', 'StepSummary.csv');
contributionFile = fullfile(projectRoot, 'Output', 'ContributionSummary.csv');
temporalOrderSummaryFile = fullfile(projectRoot, 'Output', ...
    'TemporalOrderSummary.csv');
timeLevelSummaryFile = fullfile(projectRoot, 'Output', ...
    'TimeLevelErrorSummary.csv');
figureDir = fullfile(projectRoot, 'figure');

step = readtable(stepFile, 'VariableNamingRule', 'preserve');
order = readtable(orderFile, 'VariableNamingRule', 'preserve');

n = step.iteration;
V0 = step.volume_before(1);
stepRelative = step.step_abs_drift / V0;
transportRelative = step.transport_abs_drift / V0;
reinitRelative = step.reinit_abs_jump / V0;
isReinit = logical(step.did_reinitialize);

ordinary = ~isReinit;
groups = {true(height(step), 1), ordinary, isReinit};
groupNames = {'all_steps'; 'ordinary_steps'; 'reinitialization_steps'};

meanStep = zeros(3, 1);
medianStep = zeros(3, 1);
p95Step = zeros(3, 1);
maxStep = zeros(3, 1);
for k = 1:3
    values = stepRelative(groups{k});
    meanStep(k) = mean(values);
    medianStep(k) = median(values);
    p95Step(k) = prctile(values, 95);
    maxStep(k) = max(values);
end

summary = table(groupNames, meanStep, medianStep, p95Step, maxStep, ...
    'VariableNames', {'group', 'mean_step_relative', 'median_step_relative', ...
    'p95_step_relative', 'max_step_relative'});
writetable(summary, summaryFile);

signedTransport = sum(step.transport_signed_drift) / V0;
signedReinit = sum(step.reinit_signed_jump) / V0;
signedTotal = sum(step.step_signed_drift) / V0;
absoluteTransport = sum(step.transport_abs_drift) / V0;
absoluteReinit = sum(step.reinit_abs_jump) / V0;
absoluteSourceTotal = absoluteTransport + absoluteReinit;
contributionMetric = {
    'signed_transport_over_V0';
    'signed_reinitialization_over_V0';
    'signed_total_over_V0';
    'absolute_transport_over_V0';
    'absolute_reinitialization_over_V0';
    'transport_share_of_absolute_sources';
    'reinitialization_share_of_absolute_sources';
    'number_of_reinitialization_steps'
};
contributionValue = [
    signedTransport;
    signedReinit;
    signedTotal;
    absoluteTransport;
    absoluteReinit;
    absoluteTransport / absoluteSourceTotal;
    absoluteReinit / absoluteSourceTotal;
    sum(isReinit)
];
contribution = table(contributionMetric, contributionValue, ...
    'VariableNames', {'metric', 'value'});
writetable(contribution, contributionFile);

fprintf('\nSingle-step relative-volume statistics\n');
disp(summary);
fprintf('Signed transport contribution / V0: %.8e\n', ...
    signedTransport);
fprintf('Signed reinitialization contribution / V0: %.8e\n', ...
    signedReinit);
fprintf('Final relative volume change: %.8e\n', ...
    signedTotal);
fprintf('Absolute transport accumulation / V0: %.8e\n', ...
    absoluteTransport);
fprintf('Absolute reinitialization accumulation / V0: %.8e\n', ...
    absoluteReinit);

blue = [0.000, 0.447, 0.741];
orange = [0.850, 0.325, 0.098];
red = [0.635, 0.078, 0.184];
green = [0.466, 0.674, 0.188];
gray = [0.25, 0.25, 0.25];

fig1 = figure('Units', 'centimeters', 'Position', [2, 2, 18, 14], ...
    'Color', 'w');
layout = tiledlayout(fig1, 2, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax1 = nexttile(layout);
semilogy(ax1, n, stepRelative, '-', 'Color', blue, 'LineWidth', 1.25, ...
    'DisplayName', '$E_n^{\mathrm{step}}/V_\gamma^0$');
hold(ax1, 'on');
semilogy(ax1, n, transportRelative, '--', 'Color', orange, ...
    'LineWidth', 1.15, 'DisplayName', '$E_n^{\mathrm{tr}}/V_\gamma^0$');
scatter(ax1, n(isReinit), reinitRelative(isReinit), 24, red, 'filled', ...
    'MarkerEdgeColor', 'w', 'LineWidth', 0.4, ...
    'DisplayName', '$E_n^{\mathrm{re}}/V_\gamma^0$');
grid(ax1, 'on');
ylabel(ax1, 'Relative one-step change', 'Interpreter', 'latex');
legend(ax1, 'Location', 'best', 'Interpreter', 'latex', 'NumColumns', 3);
text(ax1, 0.012, 0.90, '(a)', 'Units', 'normalized', ...
    'FontWeight', 'bold', 'FontSize', 10);

ax2 = nexttile(layout);
plot(ax2, n, cumsum(step.step_signed_drift) / V0, '-', ...
    'Color', blue, 'LineWidth', 1.35, ...
    'DisplayName', '$\sum_{k\leq n}\Delta V_k^{\mathrm{step}}/V_\gamma^0$');
hold(ax2, 'on');
plot(ax2, n, cumsum(step.transport_signed_drift) / V0, '--', ...
    'Color', orange, 'LineWidth', 1.2, ...
    'DisplayName', '$\sum_{k\leq n}\Delta V_k^{\mathrm{tr}}/V_\gamma^0$');
plot(ax2, n, cumsum(step.reinit_signed_jump) / V0, ':', ...
    'Color', red, 'LineWidth', 1.4, ...
    'DisplayName', '$\sum_{k\leq n}\Delta V_k^{\mathrm{re}}/V_\gamma^0$');
yline(ax2, 0, '-', 'Color', gray, 'HandleVisibility', 'off');
grid(ax2, 'on');
xlabel(ax2, 'Outer iteration $n$', 'Interpreter', 'latex');
ylabel(ax2, 'Cumulative signed change', 'Interpreter', 'latex');
legend(ax2, 'Location', 'best', 'Interpreter', 'latex', 'NumColumns', 1);
text(ax2, 0.012, 0.90, '(b)', 'Units', 'normalized', ...
    'FontWeight', 'bold', 'FontSize', 10);

set([ax1, ax2], 'FontName', 'Times New Roman', 'FontSize', 9, ...
    'LineWidth', 0.8, 'Box', 'on');
exportgraphics(fig1, fullfile(figureDir, 'StepVolumeDecomposition.pdf'), ...
    'ContentType', 'vector');
exportgraphics(fig1, fullfile(figureDir, 'StepVolumeDecomposition.png'), ...
    'Resolution', 300);

maxLevel = max(order.level);
if maxLevel < 3
    error('Four time levels T, T/2, T/4, T/8 are required.');
end
levels = (0:3)';
levelLabels = {'$T$', '$T/2$', '$T/4$', '$T/8$'};
pairLabels = {'$T\rightarrow T/2$', '$T/2\rightarrow T/4$', ...
    '$T/4\rightarrow T/8$'};

medianTransportError = zeros(4, 1);
q25TransportError = zeros(4, 1);
q75TransportError = zeros(4, 1);
medianReinitStepError = zeros(4, 1);
q25ReinitStepError = zeros(4, 1);
q75ReinitStepError = zeros(4, 1);
for k = 1:4
    levelRows = order(order.level == levels(k), :);
    transportValues = levelRows.transport_abs_drift / V0;
    reinitValues = levelRows.step_abs_drift(...
        logical(levelRows.did_reinitialize)) / V0;
    medianTransportError(k) = median(transportValues);
    q25TransportError(k) = prctile(transportValues, 25);
    q75TransportError(k) = prctile(transportValues, 75);
    medianReinitStepError(k) = median(reinitValues);
    q25ReinitStepError(k) = prctile(reinitValues, 25);
    q75ReinitStepError(k) = prctile(reinitValues, 75);
end

timeLevelDuration = zeros(4, 1);
for k = 1:4
    levelRows = order(order.level == levels(k), :);
    timeLevelDuration(k) = levelRows.T(1);
end
timeLevelErrorSummary = table(levels, timeLevelDuration, ...
    medianTransportError, q25TransportError, q75TransportError, ...
    medianReinitStepError, q25ReinitStepError, q75ReinitStepError);
writetable(timeLevelErrorSummary, timeLevelSummaryFile);

if ~ismember('velocity_component_linf', step.Properties.VariableNames)
    error(['StepDiagnostics.csv lacks velocity_component_linf. ', ...
        'Rerun RunFullStepStudy.ps1 before plotting.']);
end

hMax = 0.0250722;
coarseT = zeros(3, 1);
fineT = zeros(3, 1);
medianTransportOrder = zeros(3, 1);
q25TransportOrder = zeros(3, 1);
q75TransportOrder = zeros(3, 1);
medianReinitStepOrder = zeros(3, 1);
q25ReinitStepOrder = zeros(3, 1);
q75ReinitStepOrder = zeros(3, 1);
medianPrimaryTheoryOrder = zeros(3, 1);
medianSecondaryTheoryOrder = zeros(3, 1);

for k = 1:3
    rateRows = order(order.level == k, :);
    transportRates = rateRows.transport_observed_rate;
    transportRates = transportRates(isfinite(transportRates));
    reinitRates = rateRows.step_observed_rate(...
        logical(rateRows.did_reinitialize));
    reinitRates = reinitRates(isfinite(reinitRates));

    medianTransportOrder(k) = median(transportRates);
    q25TransportOrder(k) = prctile(transportRates, 25);
    q75TransportOrder(k) = prctile(transportRates, 75);
    medianReinitStepOrder(k) = median(reinitRates);
    q25ReinitStepOrder(k) = prctile(reinitRates, 25);
    q75ReinitStepOrder(k) = prctile(reinitRates, 75);

    coarseT(k) = 2 * rateRows.T(1);
    fineT(k) = rateRows.T(1);
    [matched, locations] = ismember(rateRows.iteration, step.iteration);
    if ~all(matched)
        error('Could not match shadow rows to StepDiagnostics iterations.');
    end
    theta = step.velocity_component_linf(locations);
    primaryCoarse = (hMax^2 + hMax .* theta .* coarseT(k)).^(1/3);
    primaryFine = (hMax^2 + hMax .* theta .* fineT(k)).^(1/3);
    secondaryCoarse = coarseT(k) .* ...
        (hMax^2 + theta.^2 .* coarseT(k)^2).^(1/3);
    secondaryFine = fineT(k) .* ...
        (hMax^2 + theta.^2 .* fineT(k)^2).^(1/3);
    medianPrimaryTheoryOrder(k) = median(log2(primaryCoarse ./ primaryFine));
    medianSecondaryTheoryOrder(k) = median(...
        log2(secondaryCoarse ./ secondaryFine));
end

temporalOrderSummary = table(coarseT, fineT, medianTransportOrder, ...
    q25TransportOrder, q75TransportOrder, medianReinitStepOrder, ...
    q25ReinitStepOrder, q75ReinitStepOrder, ...
    medianPrimaryTheoryOrder, medianSecondaryTheoryOrder);
writetable(temporalOrderSummary, temporalOrderSummaryFile);

fig2 = figure('Units', 'centimeters', 'Position', [2, 2, 18, 8.8], ...
    'Color', 'w');
orderLayout = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax3 = nexttile(orderLayout);
xLevels = (1:4)';
errorbar(ax3, xLevels-0.07, medianTransportError, ...
    medianTransportError-q25TransportError, ...
    q75TransportError-medianTransportError, 'o-', ...
    'Color', blue, 'MarkerFaceColor', blue, 'LineWidth', 1.2, ...
    'MarkerSize', 5, 'CapSize', 6, ...
    'DisplayName', 'Transport error');
hold(ax3, 'on');
errorbar(ax3, xLevels+0.07, medianReinitStepError, ...
    medianReinitStepError-q25ReinitStepError, ...
    q75ReinitStepError-medianReinitStepError, 's--', ...
    'Color', red, 'MarkerFaceColor', red, 'LineWidth', 1.1, ...
    'MarkerSize', 5, 'CapSize', 6, ...
    'DisplayName', 'Full-step error at reinitialization');
set(ax3, 'YScale', 'log', 'XTick', xLevels, 'XTickLabel', levelLabels, ...
    'TickLabelInterpreter', 'latex');
xlim(ax3, [0.6, 4.4]);
grid(ax3, 'on');
xlabel(ax3, 'Shadow transport duration', 'Interpreter', 'latex');
ylabel(ax3, 'Median relative error (IQR)', 'Interpreter', 'latex');
legend(ax3, 'Location', 'southwest', 'Interpreter', 'latex');
text(ax3, 0.035, 0.94, '(a)', 'Units', 'normalized', ...
    'FontWeight', 'bold', 'FontSize', 10);

ax4 = nexttile(orderLayout);
xPairs = (1:3)';
errorbar(ax4, xPairs-0.10, medianTransportOrder, ...
    medianTransportOrder-q25TransportOrder, ...
    q75TransportOrder-medianTransportOrder, 'o', ...
    'Color', blue, 'MarkerFaceColor', blue, 'LineWidth', 1.2, ...
    'MarkerSize', 5, 'CapSize', 6, ...
    'DisplayName', 'Transport numerical order');
hold(ax4, 'on');
errorbar(ax4, xPairs+0.10, medianReinitStepOrder, ...
    medianReinitStepOrder-q25ReinitStepOrder, ...
    q75ReinitStepOrder-medianReinitStepOrder, 's', ...
    'Color', red, 'MarkerFaceColor', red, 'LineWidth', 1.1, ...
    'MarkerSize', 5, 'CapSize', 6, ...
    'DisplayName', 'Reinitialization-step numerical order');
plot(ax4, xPairs, medianPrimaryTheoryOrder, 'd--', ...
    'Color', [0.10, 0.10, 0.10], 'MarkerFaceColor', 'w', ...
    'LineWidth', 1.15, 'MarkerSize', 5, ...
    'DisplayName', 'Fixed-$h$ primary bound term');
plot(ax4, xPairs, medianSecondaryTheoryOrder, '^:', ...
    'Color', [0.45, 0.45, 0.45], 'MarkerFaceColor', 'w', ...
    'LineWidth', 1.15, 'MarkerSize', 5, ...
    'DisplayName', 'Fixed-$h$ secondary bound term');
yline(ax4, 0, '-', 'Color', [0.7, 0.7, 0.7], ...
    'HandleVisibility', 'off');
set(ax4, 'XTick', xPairs, 'XTickLabel', pairLabels, ...
    'TickLabelInterpreter', 'latex');
xlim(ax4, [0.6, 3.4]);
grid(ax4, 'on');
xlabel(ax4, 'Adjacent time-level pair', 'Interpreter', 'latex');
ylabel(ax4, 'Time-interval scaling order', 'Interpreter', 'latex');
legend(ax4, 'Location', 'best', 'Interpreter', 'latex');
text(ax4, 0.035, 0.94, '(b)', 'Units', 'normalized', ...
    'FontWeight', 'bold', 'FontSize', 10);

set([ax3, ax4], 'FontName', 'Times New Roman', 'FontSize', 8.5, ...
    'LineWidth', 0.8, 'Box', 'on');
exportgraphics(fig2, fullfile(figureDir, 'StepObservedOrder.pdf'), ...
    'ContentType', 'vector');
exportgraphics(fig2, fullfile(figureDir, 'StepObservedOrder.png'), ...
    'Resolution', 300);

fprintf('Saved summary: %s\n', summaryFile);
fprintf('Saved contribution summary: %s\n', contributionFile);
fprintf('Saved temporal-order summary: %s\n', temporalOrderSummaryFile);
fprintf('Saved time-level error summary: %s\n', timeLevelSummaryFile);
fprintf('Saved publication figures in: %s\n', figureDir);

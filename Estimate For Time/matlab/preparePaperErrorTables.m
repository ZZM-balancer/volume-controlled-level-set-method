function preparePaperErrorTables()
%PREPAREPAPERERRORTABLES Build publication-ready transport/reinitialization tables.
% Transport errors contain only the volume change before reinitialization.
% Reinitialization errors contain only the jump introduced by reinitialization.

rootDir = fileparts(fileparts(mfilename('fullpath')));
outputDir = fullfile(rootDir, 'Output');

%% Fixed-h temporal refinement: T, T/2, T/4, T/8.
stepOrder = readtable(fullfile(outputDir, 'StepOrder.csv'));
stepDiag = readtable(fullfile(outputDir, 'StepDiagnostics.csv'));
levels = unique(stepOrder.level, 'sorted');
nLevels = numel(levels);

timeLevel = zeros(nLevels, 1);
Tvalue = zeros(nLevels, 1);
timeCount = zeros(nLevels, 1);
timeMax = zeros(nLevels, 1);
timeRms = zeros(nLevels, 1);
timeEnd = zeros(nLevels, 1);
timeMedian = zeros(nLevels, 1);
timeMaximumIteration = zeros(nLevels, 1);

for k = 1:nLevels
    rows = stepOrder.level == levels(k);
    sub = stepOrder(rows, :);
    [isFound, idx] = ismember(sub.iteration, stepDiag.iteration);
    assert(all(isFound), 'A StepOrder iteration is absent from StepDiagnostics.');
    before = stepDiag.volume_before(idx);
    err = sub.transport_abs_drift ./ before;
    [timeMax(k), timeRms(k), timeEnd(k), timeMedian(k), timeMaximumIteration(k)] = ...
        metrics(err, sub.iteration);
    timeLevel(k) = levels(k);
    Tvalue(k) = sub.T(1);
    timeCount(k) = numel(err);
end

timeMetrics = table(timeLevel, Tvalue, timeCount, timeMaximumIteration, timeMax, timeRms, ...
    timeEnd, timeMedian, 'VariableNames', {'level','T','sample_count', ...
    'maximum_iteration','maximum','rms','endpoint','median'});
writetable(timeMetrics, fullfile(outputDir, 'PaperTimeTransportMetrics.csv'));

[timeOrder, timePairs] = adjacentOrders(Tvalue, timeMax, timeRms, timeEnd, timeMedian);
timeOrders = table(timePairs(:,1), timePairs(:,2), Tvalue(1:end-1), Tvalue(2:end), ...
    timeOrder(:,1), timeOrder(:,2), timeOrder(:,3), timeOrder(:,4), ...
    'VariableNames', {'coarse_level','fine_level','T_coarse','T_fine', ...
    'maximum_order','rms_order','endpoint_order','median_order'});
writetable(timeOrders, fullfile(outputDir, 'PaperTimeTransportOrders.csv'));

%% Four-mesh, fixed-T, 300-step transport and reinitialization statistics.
meshScales = [2; 3; 5; 7];
nMeshes = numel(meshScales);
hmax = zeros(nMeshes,1);
meshCount = zeros(nMeshes,1);
meshMax = zeros(nMeshes,1);
meshRms = zeros(nMeshes,1);
meshEnd = zeros(nMeshes,1);
meshMedian = zeros(nMeshes,1);
reinitCount = zeros(nMeshes,1);
reinitMax = zeros(nMeshes,1);
reinitRms = zeros(nMeshes,1);
reinitEnd = zeros(nMeshes,1);
reinitMedian = zeros(nMeshes,1);
reinitShare = zeros(nMeshes,1);
reinitSignedNet = zeros(nMeshes,1);
meshMaximumIteration = zeros(nMeshes,1);
initialVelocityLinf = zeros(nMeshes,1);
medianVelocityLinf = zeros(nMeshes,1);

for k = 1:nMeshes
    fileName = sprintf('CoupledMeshFixedT_n%d.csv', meshScales(k));
    data = readtable(fullfile(outputDir, fileName));
    hmax(k) = data.hmax(1);

    transportErr = data.transport_abs_drift ./ data.volume_before;
    [meshMax(k), meshRms(k), meshEnd(k), meshMedian(k), meshMaximumIteration(k)] = ...
        metrics(transportErr, data.iteration);
    meshCount(k) = numel(transportErr);
    initialVelocityLinf(k) = data.velocity_component_linf(1);
    medianVelocityLinf(k) = median(data.velocity_component_linf);

    reRows = data.did_reinitialize == 1;
    reErr = data.reinit_abs_jump(reRows) ./ data.volume_before(reRows);
    [reinitMax(k), reinitRms(k), reinitEnd(k), reinitMedian(k)] = ...
        metrics(reErr, data.iteration(reRows));
    reinitCount(k) = numel(reErr);

    absTransport = sum(data.transport_abs_drift);
    absReinit = sum(data.reinit_abs_jump);
    reinitShare(k) = absReinit / (absTransport + absReinit);
    reinitSignedNet(k) = sum(data.reinit_signed_jump) / data.volume_before(1);
end

meshMetrics = table(meshScales, hmax, meshCount, meshMaximumIteration, meshMax, meshRms, ...
    meshEnd, meshMedian, 'VariableNames', {'boundary_scale','hmax','step_count', ...
    'maximum_iteration','maximum','rms','endpoint','median'});
writetable(meshMetrics, fullfile(outputDir, 'PaperMeshTransportMetrics.csv'));

initialToMedianVelocityRatio = initialVelocityLinf ./ medianVelocityLinf;
initialStepDiagnostic = table(meshScales, hmax, meshMaximumIteration, ...
    initialVelocityLinf, medianVelocityLinf, initialToMedianVelocityRatio, ...
    'VariableNames', {'boundary_scale','hmax','maximum_error_iteration', ...
    'initial_velocity_linf','median_velocity_linf','initial_to_median_velocity_ratio'});
writetable(initialStepDiagnostic, fullfile(outputDir, 'PaperInitialStepDiagnostic.csv'));

[meshOrder, meshPairs] = adjacentOrders(hmax, meshMax, meshRms, meshEnd, meshMedian);
meshOrders = table(meshPairs(:,1), meshPairs(:,2), hmax(1:end-1), hmax(2:end), ...
    meshOrder(:,1), meshOrder(:,2), meshOrder(:,3), meshOrder(:,4), ...
    'VariableNames', {'coarse_scale','fine_scale','h_coarse','h_fine', ...
    'maximum_order','rms_order','endpoint_order','median_order'});
writetable(meshOrders, fullfile(outputDir, 'PaperMeshTransportOrders.csv'));

reinitMetrics = table(meshScales, hmax, reinitCount, reinitMax, reinitRms, ...
    reinitEnd, reinitMedian, reinitShare, reinitSignedNet, ...
    'VariableNames', {'boundary_scale','hmax','event_count','maximum','rms', ...
    'endpoint','median','absolute_contribution_share','signed_net_over_initial_volume'});
writetable(reinitMetrics, fullfile(outputDir, 'PaperReinitializationMetrics.csv'));

[reinitOrder, reinitPairs] = adjacentOrders(hmax, reinitMax, reinitRms, ...
    reinitEnd, reinitMedian);
reinitOrders = table(reinitPairs(:,1), reinitPairs(:,2), hmax(1:end-1), hmax(2:end), ...
    reinitOrder(:,1), reinitOrder(:,2), reinitOrder(:,3), reinitOrder(:,4), ...
    'VariableNames', {'coarse_scale','fine_scale','h_coarse','h_fine', ...
    'maximum_order','rms_order','endpoint_order','median_order'});
writetable(reinitOrders, fullfile(outputDir, 'PaperReinitializationOrders.csv'));

fprintf('Wrote seven paper tables to %s\n', outputDir);
end

function [maximumValue, rmsValue, endpointValue, medianValue, maximumIteration] = metrics(err, iteration)
assert(~isempty(err), 'Metric input must not be empty.');
[maximumValue, maximumIndex] = max(err);
maximumIteration = iteration(maximumIndex);
rmsValue = sqrt(mean(err.^2));
[~, endpointIndex] = max(iteration);
endpointValue = err(endpointIndex);
medianValue = median(err);
end

function [orders, pairs] = adjacentOrders(scale, maximumValue, rmsValue, endpointValue, medianValue)
values = [maximumValue, rmsValue, endpointValue, medianValue];
nPairs = numel(scale) - 1;
orders = zeros(nPairs, size(values,2));
pairs = [(1:nPairs)', (2:nPairs+1)'];
for k = 1:nPairs
    orders(k,:) = log(values(k,:) ./ values(k+1,:)) ./ log(scale(k) / scale(k+1));
end
end

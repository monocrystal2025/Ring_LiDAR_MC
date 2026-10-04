function [summary, comparisons] = multiuav_summarize(results, cfg)
%MULTIUAV_SUMMARIZE Convert scene-level results into paper-ready metrics.

arguments
    results (1,:) cell
    cfg (1,1) struct
end

rowCount = numel(results) * numel(cfg.scene.timelyThresholds_s);
Region = strings(rowCount, 1);
Trajectory = strings(rowCount, 1);
Scenario = strings(rowCount, 1);
Beam = strings(rowCount, 1);
BeamLabel = strings(rowCount, 1);
TargetCount = zeros(rowCount, 1);
wD_mrad = zeros(rowCount, 1);
wd_mrad = zeros(rowCount, 1);
Threshold_s = zeros(rowCount, 1);
PAll = zeros(rowCount, 1);
MeanRecall = zeros(rowCount, 1);
PRecall80 = zeros(rowCount, 1);
MeanCompleteTime_s = nan(rowCount, 1);
MedianCompleteTime_s = nan(rowCount, 1);
CompletionConditionalN = zeros(rowCount, 1);
MeanRepeatedEncounters = zeros(rowCount, 1);
MeanResolvedHits = zeros(rowCount, 1);
MeanMaxTargetsPerPulse = zeros(rowCount, 1);
MeanMultiTargetPulses = zeros(rowCount, 1);
RangeAmbiguityRate = zeros(rowCount, 1);
MeanExitWithoutDetection = zeros(rowCount, 1);
CycleTime_s = zeros(rowCount, 1);
SceneN = zeros(rowCount, 1);

row = 0;
for resultIndex = 1:numel(results)
    result = results{resultIndex};
    M = result.targetCount;
    for threshold = cfg.scene.timelyThresholds_s
        if threshold > cfg.system.maxObservation_s + eps
            continue;
        end
        row = row + 1;
        timely = result.detect & result.firstTime_s <= threshold;
        recall = sum(timely, 2) ./ M;
        allCaptured = all(timely, 2);
        completeTime = max(result.firstTime_s, [], 2);
        completeTime(~allCaptured) = NaN;
        validComplete = completeTime(isfinite(completeTime));
        resolvedTotal = sum(result.resolvedHitCount, 'all');
        ambiguousTotal = sum(result.ambiguousHitCount, 'all');

        Region(row) = result.region;
        Trajectory(row) = result.trajectory;
        Scenario(row) = result.scenario;
        Beam(row) = result.beam;
        BeamLabel(row) = result.beamLabel;
        TargetCount(row) = M;
        wD_mrad(row) = result.wD_mrad;
        wd_mrad(row) = result.wd_mrad;
        Threshold_s(row) = threshold;
        PAll(row) = mean(allCaptured);
        MeanRecall(row) = mean(recall);
        PRecall80(row) = mean(recall >= 0.8);
        if ~isempty(validComplete)
            MeanCompleteTime_s(row) = mean(validComplete);
            MedianCompleteTime_s(row) = median(validComplete);
        end
        CompletionConditionalN(row) = numel(validComplete);
        MeanRepeatedEncounters(row) = mean(max(result.encounterCount - 1, 0), 'all');
        MeanResolvedHits(row) = mean(result.resolvedHitCount, 'all');
        MeanMaxTargetsPerPulse(row) = mean(result.maxResolvedPerPulse);
        MeanMultiTargetPulses(row) = mean(result.multiTargetPulseCount);
        RangeAmbiguityRate(row) = ambiguousTotal / ...
            max(resolvedTotal + ambiguousTotal, 1);
        MeanExitWithoutDetection(row) = mean(result.exitWithoutDetection, 'all');
        CycleTime_s(row) = result.cycleTime_s;
        SceneN(row) = size(result.detect, 1);
    end
end

summary = table(Region(1:row), Trajectory(1:row), Scenario(1:row), ...
    Beam(1:row), BeamLabel(1:row), TargetCount(1:row), ...
    wD_mrad(1:row), wd_mrad(1:row), Threshold_s(1:row), ...
    PAll(1:row), MeanRecall(1:row), PRecall80(1:row), ...
    MeanCompleteTime_s(1:row), MedianCompleteTime_s(1:row), ...
    CompletionConditionalN(1:row), MeanRepeatedEncounters(1:row), ...
    MeanResolvedHits(1:row), MeanMaxTargetsPerPulse(1:row), ...
    MeanMultiTargetPulses(1:row), RangeAmbiguityRate(1:row), ...
    MeanExitWithoutDetection(1:row), CycleTime_s(1:row), SceneN(1:row), ...
    'VariableNames', {'Region', 'Trajectory', 'Scenario', 'Beam', ...
    'BeamLabel', 'TargetCount', 'wD_mrad', 'wd_mrad', 'Threshold_s', ...
    'PAll', 'MeanRecall', 'PRecall80', 'MeanCompleteTime_s', ...
    'MedianCompleteTime_s', 'CompletionConditionalN', ...
    'MeanRepeatedEncounters', 'MeanResolvedHits', ...
    'MeanMaxTargetsPerPulse', 'MeanMultiTargetPulses', ...
    'RangeAmbiguityRate', 'MeanExitWithoutDetection', ...
    'CycleTime_s', 'SceneN'});

comparisons = addSceneConfidenceIntervals(buildComparisons(summary), results);
end

function comparisons = buildComparisons(summary)
keyNames = {'Region', 'Trajectory', 'Scenario', 'TargetCount', 'Threshold_s'};
comparisonRows = cell(0, 15);
pairDefinitions = {
    "EqualWidth125", "line", 125, "ring", 125;
    "SingleTargetReference", "line", 195, "ring", 125};

for definitionIndex = 1:size(pairDefinitions, 1)
    definition = pairDefinitions(definitionIndex, :);
    lineRows = summary.Beam == definition{2} & ...
        abs(summary.wD_mrad - definition{3}) < 1e-9;
    ringRows = summary.Beam == definition{4} & ...
        abs(summary.wD_mrad - definition{5}) < 1e-9;
    lineTable = summary(lineRows, :);
    ringTable = summary(ringRows, :);
    for lineIndex = 1:height(lineTable)
        match = true(height(ringTable), 1);
        for keyIndex = 1:numel(keyNames)
            key = keyNames{keyIndex};
            match = match & ringTable.(key) == lineTable.(key)(lineIndex);
        end
        ringIndex = find(match, 1);
        if isempty(ringIndex)
            continue;
        end
        L = lineTable(lineIndex, :);
        R = ringTable(ringIndex, :);
        comparisonRows(end + 1, :) = {definition{1}, L.Region, ... %#ok<AGROW>
            L.Trajectory, L.Scenario, L.TargetCount, L.Threshold_s, ...
            L.wD_mrad, R.wD_mrad, 100 * (R.PAll - L.PAll), ...
            100 * (R.MeanRecall - L.MeanRecall), ...
            R.MeanCompleteTime_s - L.MeanCompleteTime_s, ...
            R.MeanRepeatedEncounters - L.MeanRepeatedEncounters, ...
            R.MeanMaxTargetsPerPulse - L.MeanMaxTargetsPerPulse, ...
            R.MeanMultiTargetPulses - L.MeanMultiTargetPulses, ...
            R.RangeAmbiguityRate - L.RangeAmbiguityRate};
    end
end

variableNames = {'Comparison', 'Region', 'Trajectory', 'Scenario', ...
    'TargetCount', 'Threshold_s', 'LineWD_mrad', 'RingWD_mrad', ...
    'DeltaPAll_pp', 'DeltaMeanRecall_pp', 'DeltaMeanCompleteTime_s', ...
    'DeltaRepeatedEncounters', 'DeltaMaxTargetsPerPulse', ...
    'DeltaMultiTargetPulses', 'DeltaRangeAmbiguityRate'};
if isempty(comparisonRows)
    comparisons = cell2table(cell(0, numel(variableNames)), ...
        'VariableNames', variableNames);
else
    comparisons = cell2table(comparisonRows, 'VariableNames', variableNames);
end
end

function comparisons = addSceneConfidenceIntervals(comparisons, results)
rowCount = height(comparisons);
comparisons.DeltaPAll_CILow_pp = nan(rowCount, 1);
comparisons.DeltaPAll_CIHigh_pp = nan(rowCount, 1);
comparisons.DeltaMeanRecall_CILow_pp = nan(rowCount, 1);
comparisons.DeltaMeanRecall_CIHigh_pp = nan(rowCount, 1);
comparisons.DeltaCompleteTime_CILow_s = nan(rowCount, 1);
comparisons.DeltaCompleteTime_CIHigh_s = nan(rowCount, 1);
comparisons.PairedCompleteSceneN = zeros(rowCount, 1);

for rowIndex = 1:rowCount
    row = comparisons(rowIndex, :);
    lineResult = findResult(results, row, "line", row.LineWD_mrad);
    ringResult = findResult(results, row, "ring", row.RingWD_mrad);
    threshold = row.Threshold_s;
    lineTimely = lineResult.detect & lineResult.firstTime_s <= threshold;
    ringTimely = ringResult.detect & ringResult.firstTime_s <= threshold;
    allDifference = double(all(ringTimely, 2)) - ...
        double(all(lineTimely, 2));
    recallDifference = mean(double(ringTimely), 2) - ...
        mean(double(lineTimely), 2);
    comparisons{rowIndex, {'DeltaPAll_CILow_pp', ...
        'DeltaPAll_CIHigh_pp'}} = 100 .* meanConfidence(allDifference);
    comparisons{rowIndex, {'DeltaMeanRecall_CILow_pp', ...
        'DeltaMeanRecall_CIHigh_pp'}} = ...
        100 .* meanConfidence(recallDifference);

    bothComplete = all(lineTimely, 2) & all(ringTimely, 2);
    completeDifference = max(ringResult.firstTime_s(bothComplete, :), ...
        [], 2) - max(lineResult.firstTime_s(bothComplete, :), [], 2);
    comparisons.PairedCompleteSceneN(rowIndex) = numel(completeDifference);
    if ~isempty(completeDifference)
        comparisons{rowIndex, {'DeltaCompleteTime_CILow_s', ...
            'DeltaCompleteTime_CIHigh_s'}} = ...
            meanConfidence(completeDifference);
    end
end
end

function result = findResult(results, row, beamName, width_mrad)
matchIndex = 0;
for resultIndex = 1:numel(results)
    candidate = results{resultIndex};
    match = candidate.region == row.Region & ...
        candidate.trajectory == row.Trajectory & ...
        candidate.scenario == row.Scenario & ...
        candidate.targetCount == row.TargetCount & ...
        candidate.beam == beamName & ...
        abs(candidate.wD_mrad - width_mrad) < 1e-9;
    if match
        matchIndex = resultIndex;
        break;
    end
end
if matchIndex == 0
    error('multiuav_summarize:MissingPairedResult', ...
        'Could not find a paired scene-level result.');
end
result = results{matchIndex};
end

function interval = meanConfidence(values)
values = values(isfinite(values));
if isempty(values)
    interval = [NaN, NaN];
    return;
end
meanValue = mean(values);
halfWidth = 1.96 * std(values, 0) / sqrt(numel(values));
interval = [meanValue - halfWidth, meanValue + halfWidth];
end

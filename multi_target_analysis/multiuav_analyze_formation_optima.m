function [optima, comparisons] = multiuav_analyze_formation_optima(results, cfg)
%MULTIUAV_ANALYZE_FORMATION_OPTIMA Compare re-optimized line/ring beams.
% The discrete optimum is selected over the manuscript wD grid. A paired
% K-fold estimate re-selects wD on each training fold and evaluates it on
% held-out scenes, avoiding an optimistic same-sample maximum comparison.

arguments
    results (1,:) cell
    cfg (1,1) struct
end

threshold = cfg.optimization.threshold_s;
foldCount = cfg.optimization.cvFolds;
metadata = resultMetadata(results);
groupKeys = unique(metadata(:, {'Region', 'Trajectory', 'Scenario', ...
    'TargetCount'}), 'rows', 'stable');
optimumRows = cell(0, 18);
comparisonRows = cell(0, 23);

for groupIndex = 1:height(groupKeys)
    key = groupKeys(groupIndex, :);
    line = extractBeamData(results, key, "line", threshold);
    ring = extractBeamData(results, key, "ring", threshold);
    assertPaired(line, ring);

    lineBest = selectOptimalIndex(line.pAll, line.meanRecall, line.wD_mrad);
    ringBest = selectOptimalIndex(ring.pAll, ring.meanRecall, ring.wD_mrad);
    [lineCV, lineSelectedWD] = crossValidatedPerformance(line, foldCount);
    [ringCV, ringSelectedWD] = crossValidatedPerformance(ring, foldCount);

    lineInterval = wilsonInterval(sum(line.allCaptured(:, lineBest)), line.N);
    ringInterval = wilsonInterval(sum(ring.allCaptured(:, ringBest)), ring.N);
    optimumRows(end + 1, :) = makeOptimumRow( ... %#ok<AGROW>
        key, line, lineBest, lineInterval, lineCV, lineSelectedWD, threshold);
    optimumRows(end + 1, :) = makeOptimumRow( ... %#ok<AGROW>
        key, ring, ringBest, ringInterval, ringCV, ringSelectedWD, threshold);

    deltaCV = ringCV - lineCV;
    deltaInterval = pairedMeanInterval(deltaCV);
    apparentDelta = double(ring.allCaptured(:, ringBest)) - ...
        double(line.allCaptured(:, lineBest));
    apparentInterval = pairedMeanInterval(apparentDelta);
    apparentWinner = classifyWinner(apparentInterval);
    winner = "Indistinguishable";
    if deltaInterval(1) > 0
        winner = "Ring";
    elseif deltaInterval(2) < 0
        winner = "Line";
    end
    comparisonRows(end + 1, :) = {key.Region, key.Trajectory, ... %#ok<AGROW>
        key.Scenario, key.TargetCount, threshold, ...
        line.wD_mrad(lineBest), ring.wD_mrad(ringBest), ...
        line.pAll(lineBest), ring.pAll(ringBest), ...
        100 * (ring.pAll(ringBest) - line.pAll(lineBest)), ...
        100 * apparentInterval(1), 100 * apparentInterval(2), ...
        apparentWinner, ...
        mean(lineCV), mean(ringCV), 100 * mean(deltaCV), ...
        100 * deltaInterval(1), 100 * deltaInterval(2), winner, ...
        ring.wD_mrad(ringBest) / line.wD_mrad(lineBest), ...
        median(lineSelectedWD), median(ringSelectedWD), line.N};
end

optima = cell2table(optimumRows, 'VariableNames', { ...
    'Region', 'Trajectory', 'Scenario', 'TargetCount', 'Beam', ...
    'Threshold_s', 'BestWD_mrad', 'MaxPAll', 'PAllCILow', ...
    'PAllCIHigh', 'MeanRecallAtBest', 'MeanCompleteTimeAtBest_s', ...
    'RangeAmbiguityAtBest', 'CVPAll', 'CVSelectedWDMedian_mrad', ...
    'CVSelectedWDMin_mrad', 'CVSelectedWDMax_mrad', 'SceneN'});
comparisons = cell2table(comparisonRows, 'VariableNames', { ...
    'Region', 'Trajectory', 'Scenario', 'TargetCount', 'Threshold_s', ...
    'LineBestWD_mrad', 'RingBestWD_mrad', 'LineMaxPAll', ...
    'RingMaxPAll', 'ApparentDeltaPAll_pp', ...
    'ApparentDeltaPAllCILow_pp', 'ApparentDeltaPAllCIHigh_pp', ...
    'ApparentWinner95', 'LineCVPAll', 'RingCVPAll', ...
    'DeltaCVPAll_pp', 'DeltaCVPAllCILow_pp', ...
    'DeltaCVPAllCIHigh_pp', 'Winner95', 'BestWDRatioRingToLine', ...
    'LineCVSelectedWDMedian_mrad', 'RingCVSelectedWDMedian_mrad', ...
    'SceneN'});
end

function winner = classifyWinner(interval)
winner = "Indistinguishable";
if interval(1) > 0
    winner = "Ring";
elseif interval(2) < 0
    winner = "Line";
end
end

function metadata = resultMetadata(results)
count = numel(results);
Region = strings(count, 1);
Trajectory = strings(count, 1);
Scenario = strings(count, 1);
TargetCount = zeros(count, 1);
Beam = strings(count, 1);
wD_mrad = zeros(count, 1);
for index = 1:count
    item = results{index};
    Region(index) = item.region;
    Trajectory(index) = item.trajectory;
    Scenario(index) = item.scenario;
    TargetCount(index) = item.targetCount;
    Beam(index) = item.beam;
    wD_mrad(index) = item.wD_mrad;
end
metadata = table(Region, Trajectory, Scenario, TargetCount, Beam, wD_mrad);
end

function data = extractBeamData(results, key, beamName, threshold)
matching = false(1, numel(results));
for index = 1:numel(results)
    item = results{index};
    matching(index) = item.region == key.Region & ...
        item.trajectory == key.Trajectory & ...
        item.scenario == key.Scenario & ...
        item.targetCount == key.TargetCount & item.beam == beamName;
end
beamResults = results(matching);
if isempty(beamResults)
    error('multiuav_analyze_formation_optima:MissingBeam', ...
        'Each group must contain both line and ring results.');
end
widths = cellfun(@(x) x.wD_mrad, beamResults);
[widths, order] = sort(widths);
beamResults = beamResults(order);
N = size(beamResults{1}.detect, 1);
widthCount = numel(widths);
allCaptured = false(N, widthCount);
recall = zeros(N, widthCount);
meanCompleteTime = nan(1, widthCount);
rangeAmbiguity = zeros(1, widthCount);
scenarioId = beamResults{1}.scenarioId;
for widthIndex = 1:widthCount
    item = beamResults{widthIndex};
    if ~isequal(item.scenarioId, scenarioId)
        error('multiuav_analyze_formation_optima:UnpairedScenes', ...
            'Scenario IDs must match at every wD.');
    end
    timely = item.detect & item.firstTime_s <= threshold;
    allCaptured(:, widthIndex) = all(timely, 2);
    recall(:, widthIndex) = mean(double(timely), 2);
    completeTime = max(item.firstTime_s, [], 2);
    completeTime(~allCaptured(:, widthIndex)) = NaN;
    meanCompleteTime(widthIndex) = mean(completeTime, 'omitnan');
    ambiguous = sum(item.ambiguousHitCount, 'all');
    resolved = sum(item.resolvedHitCount, 'all');
    rangeAmbiguity(widthIndex) = ambiguous / max(ambiguous + resolved, 1);
end
data = struct('beam', beamName, 'wD_mrad', widths, ...
    'allCaptured', allCaptured, 'recall', recall, ...
    'pAll', mean(allCaptured, 1), 'meanRecall', mean(recall, 1), ...
    'meanCompleteTime_s', meanCompleteTime, ...
    'rangeAmbiguity', rangeAmbiguity, 'scenarioId', scenarioId, 'N', N);
end

function assertPaired(line, ring)
if line.N ~= ring.N || ~isequal(line.scenarioId, ring.scenarioId)
    error('multiuav_analyze_formation_optima:UnpairedBeams', ...
        'Line and ring curves must use the same scene IDs.');
end
end

function index = selectOptimalIndex(pAll, meanRecall, widths)
bestProbability = max(pAll);
candidates = find(abs(pAll - bestProbability) <= 10 * eps(max(1, bestProbability)));
bestRecall = max(meanRecall(candidates));
candidates = candidates(abs(meanRecall(candidates) - bestRecall) <= ...
    10 * eps(max(1, bestRecall)));
[~, localIndex] = min(widths(candidates));
index = candidates(localIndex);
end

function [heldOut, selectedWD] = crossValidatedPerformance(data, foldCount)
foldCount = min(foldCount, data.N);
fold = mod((1:data.N).' - 1, foldCount) + 1;
heldOut = nan(data.N, 1);
selectedWD = nan(foldCount, 1);
for foldIndex = 1:foldCount
    test = fold == foldIndex;
    train = ~test;
    pAll = mean(data.allCaptured(train, :), 1);
    recall = mean(data.recall(train, :), 1);
    optimumIndex = selectOptimalIndex(pAll, recall, data.wD_mrad);
    heldOut(test) = double(data.allCaptured(test, optimumIndex));
    selectedWD(foldIndex) = data.wD_mrad(optimumIndex);
end
end

function row = makeOptimumRow(key, data, bestIndex, interval, ...
    cvOutcome, selectedWD, threshold)
row = {key.Region, key.Trajectory, key.Scenario, key.TargetCount, ...
    data.beam, threshold, data.wD_mrad(bestIndex), data.pAll(bestIndex), ...
    interval(1), interval(2), data.meanRecall(bestIndex), ...
    data.meanCompleteTime_s(bestIndex), data.rangeAmbiguity(bestIndex), ...
    mean(cvOutcome), median(selectedWD), min(selectedWD), ...
    max(selectedWD), data.N};
end

function interval = pairedMeanInterval(values)
meanValue = mean(values);
halfWidth = 1.96 * std(values, 0) / sqrt(numel(values));
interval = [max(-1, meanValue - halfWidth), ...
    min(1, meanValue + halfWidth)];
end

function interval = wilsonInterval(successes, total)
z = 1.96;
probability = successes / total;
denominator = 1 + z^2 / total;
centre = (probability + z^2 / (2 * total)) / denominator;
halfWidth = z / denominator * sqrt(probability * (1 - probability) / total + ...
    z^2 / (4 * total^2));
interval = [max(0, centre - halfWidth), min(1, centre + halfWidth)];
end

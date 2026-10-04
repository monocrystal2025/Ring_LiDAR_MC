function [summary, comparisons, effects] = sens_analyze_results(rawResults, cfg)
%SENS_ANALYZE_RESULTS Summarize paired sensitivity simulations.
%
% rawResults is a struct array produced by run_sensitivity_analysis. Each
% entry represents one factor level, beam, and trajectory and contains the
% paired per-scenario vectors detect and firstTime_s.

arguments
    rawResults (1,:) struct
    cfg (1,1) struct
end

required = {'caseId', 'factor', 'factorLabel', 'level', 'levelValue', ...
    'unit', 'beam', 'trajectory', 'detect', 'firstTime_s'};
if isempty(rawResults) || ~all(isfield(rawResults, required))
    error('sens_analyze_results:InvalidRawResults', ...
        'rawResults is empty or missing required fields.');
end

summary = buildSummary(rawResults, cfg);
comparisons = buildComparisons(rawResults, cfg);
effects = buildEffects(comparisons, rawResults, cfg);
end

function summary = buildSummary(rawResults, cfg)
nRows = numel(rawResults);
CaseId = strings(nRows, 1);
Factor = strings(nRows, 1);
FactorLabel = strings(nRows, 1);
Level = strings(nRows, 1);
LevelValue = nan(nRows, 1);
Unit = strings(nRows, 1);
Beam = strings(nRows, 1);
Trajectory = strings(nRows, 1);
N = zeros(nRows, 1);
NDetected = zeros(nRows, 1);
PDetect = nan(nRows, 1);
PDetect_CI_Low = nan(nRows, 1);
PDetect_CI_High = nan(nRows, 1);
MeanWarning_s = nan(nRows, 1);
MeanWarning_CI_Low_s = nan(nRows, 1);
MeanWarning_CI_High_s = nan(nRows, 1);
PTimely5s = nan(nRows, 1);
PTimely5s_CI_Low = nan(nRows, 1);
PTimely5s_CI_High = nan(nRows, 1);

for rowIdx = 1:nRows
    item = rawResults(rowIdx);
    detected = logical(item.detect(:));
    firstTime = item.firstTime_s(:);
    validWarning = detected & isfinite(firstTime);
    timely = validWarning & firstTime <= cfg.timelyThreshold_s;
    sampleCount = numel(detected);
    detectedCount = nnz(detected);

    CaseId(rowIdx) = string(item.caseId);
    Factor(rowIdx) = string(item.factor);
    FactorLabel(rowIdx) = string(item.factorLabel);
    Level(rowIdx) = string(item.level);
    LevelValue(rowIdx) = item.levelValue;
    Unit(rowIdx) = string(item.unit);
    Beam(rowIdx) = string(item.beam);
    Trajectory(rowIdx) = string(item.trajectory);
    N(rowIdx) = sampleCount;
    NDetected(rowIdx) = detectedCount;
    PDetect(rowIdx) = detectedCount / sampleCount;
    [PDetect_CI_Low(rowIdx), PDetect_CI_High(rowIdx)] = ...
        wilsonInterval(detectedCount, sampleCount);
    if any(validWarning)
        warningTimes = firstTime(validWarning);
        MeanWarning_s(rowIdx) = mean(warningTimes);
        [MeanWarning_CI_Low_s(rowIdx), MeanWarning_CI_High_s(rowIdx)] = ...
            bootstrapMeanInterval(warningTimes, cfg.bootstrapN, cfg.seed + rowIdx);
    end
    PTimely5s(rowIdx) = nnz(timely) / sampleCount;
    [PTimely5s_CI_Low(rowIdx), PTimely5s_CI_High(rowIdx)] = ...
        wilsonInterval(nnz(timely), sampleCount);
end

summary = table(CaseId, Factor, FactorLabel, Level, LevelValue, Unit, Beam, ...
    Trajectory, N, NDetected, PDetect, PDetect_CI_Low, PDetect_CI_High, ...
    MeanWarning_s, MeanWarning_CI_Low_s, MeanWarning_CI_High_s, ...
    PTimely5s, PTimely5s_CI_Low, PTimely5s_CI_High);
end

function comparisons = buildComparisons(rawResults, cfg)
caseIds = string({rawResults.caseId}).';
trajectories = string({rawResults.trajectory}).';
[groupKeys, ~, groupIndex] = unique(caseIds + "|" + trajectories, 'stable');
nGroups = numel(groupKeys);

CaseId = strings(nGroups, 1);
Factor = strings(nGroups, 1);
FactorLabel = strings(nGroups, 1);
Level = strings(nGroups, 1);
LevelValue = nan(nGroups, 1);
Unit = strings(nGroups, 1);
Trajectory = strings(nGroups, 1);
DeltaP_RingLine = nan(nGroups, 1);
DeltaP_RingLine_CI_Low = nan(nGroups, 1);
DeltaP_RingLine_CI_High = nan(nGroups, 1);
DeltaWarning_LineMinusRing_s = nan(nGroups, 1);
DeltaWarning_LineMinusRing_CI_Low_s = nan(nGroups, 1);
DeltaWarning_LineMinusRing_CI_High_s = nan(nGroups, 1);
DeltaP_RingPoint = nan(nGroups, 1);
DeltaWarning_PointMinusRing_s = nan(nGroups, 1);

for groupIdx = 1:nGroups
    members = find(groupIndex == groupIdx);
    beams = lower(string({rawResults(members).beam}));
    ringIdx = members(beams == "ring");
    lineIdx = members(beams == "line");
    pointIdx = members(beams == "point");
    if numel(ringIdx) ~= 1 || numel(lineIdx) ~= 1 || numel(pointIdx) ~= 1
        error('sens_analyze_results:IncompletePairedGroup', ...
            'Each case and trajectory must contain one point, line, and ring result.');
    end

    ring = rawResults(ringIdx);
    line = rawResults(lineIdx);
    point = rawResults(pointIdx);
    verifyPairing(ring, line, point);

    CaseId(groupIdx) = string(ring.caseId);
    Factor(groupIdx) = string(ring.factor);
    FactorLabel(groupIdx) = string(ring.factorLabel);
    Level(groupIdx) = string(ring.level);
    LevelValue(groupIdx) = ring.levelValue;
    Unit(groupIdx) = string(ring.unit);
    Trajectory(groupIdx) = string(ring.trajectory);

    ringDetect = logical(ring.detect(:));
    lineDetect = logical(line.detect(:));
    pointDetect = logical(point.detect(:));
    DeltaP_RingLine(groupIdx) = mean(ringDetect) - mean(lineDetect);
    DeltaP_RingPoint(groupIdx) = mean(ringDetect) - mean(pointDetect);
    DeltaWarning_LineMinusRing_s(groupIdx) = conditionalMean(line) - conditionalMean(ring);
    DeltaWarning_PointMinusRing_s(groupIdx) = conditionalMean(point) - conditionalMean(ring);

    [probabilityCI, warningCI] = pairedBootstrap(ring, line, ...
        cfg.bootstrapN, cfg.seed + 10000 + groupIdx);
    DeltaP_RingLine_CI_Low(groupIdx) = probabilityCI(1);
    DeltaP_RingLine_CI_High(groupIdx) = probabilityCI(2);
    DeltaWarning_LineMinusRing_CI_Low_s(groupIdx) = warningCI(1);
    DeltaWarning_LineMinusRing_CI_High_s(groupIdx) = warningCI(2);
end

comparisons = table(CaseId, Factor, FactorLabel, Level, LevelValue, Unit, ...
    Trajectory, DeltaP_RingLine, DeltaP_RingLine_CI_Low, ...
    DeltaP_RingLine_CI_High, DeltaWarning_LineMinusRing_s, ...
    DeltaWarning_LineMinusRing_CI_Low_s, ...
    DeltaWarning_LineMinusRing_CI_High_s, DeltaP_RingPoint, ...
    DeltaWarning_PointMinusRing_s);
end

function effects = buildEffects(comparisons, rawResults, cfg)
if ~any(comparisons.Level == "Low") || ~any(comparisons.Level == "High")
    effects = table( ...
        strings(0,1), strings(0,1), strings(0,1), zeros(0,1), ...
        zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
        zeros(0,1), zeros(0,1), zeros(0,1), ...
        'VariableNames', {'Factor', 'FactorLabel', 'Trajectory', ...
        'ProbabilityEffect_pp', 'ProbabilityEffect_CI_Low_pp', ...
        'ProbabilityEffect_CI_High_pp', 'WarningEffect_s', ...
        'WarningEffect_CI_Low_s', 'WarningEffect_CI_High_s', ...
        'RingPointProbabilityEffect_pp', 'RingPointWarningEffect_s'});
    return;
end
factorNames = unique(comparisons.Factor(comparisons.Factor ~= "baseline"), 'stable');
trajectories = unique(comparisons.Trajectory, 'stable');
nRows = numel(factorNames) * numel(trajectories);

Factor = strings(nRows, 1);
FactorLabel = strings(nRows, 1);
Trajectory = strings(nRows, 1);
ProbabilityEffect_pp = nan(nRows, 1);
ProbabilityEffect_CI_Low_pp = nan(nRows, 1);
ProbabilityEffect_CI_High_pp = nan(nRows, 1);
WarningEffect_s = nan(nRows, 1);
WarningEffect_CI_Low_s = nan(nRows, 1);
WarningEffect_CI_High_s = nan(nRows, 1);
RingPointProbabilityEffect_pp = nan(nRows, 1);
RingPointWarningEffect_s = nan(nRows, 1);

rowIdx = 0;
for factorIdx = 1:numel(factorNames)
    for trajectoryIdx = 1:numel(trajectories)
        rowIdx = rowIdx + 1;
        factorName = factorNames(factorIdx);
        trajectory = trajectories(trajectoryIdx);
        low = comparisons(comparisons.Factor == factorName & ...
            comparisons.Trajectory == trajectory & comparisons.Level == "Low", :);
        high = comparisons(comparisons.Factor == factorName & ...
            comparisons.Trajectory == trajectory & comparisons.Level == "High", :);
        if height(low) ~= 1 || height(high) ~= 1
            error('sens_analyze_results:MissingFactorEndpoint', ...
                'Factor %s trajectory %s requires one Low and one High row.', ...
                factorName, trajectory);
        end
        Factor(rowIdx) = factorName;
        FactorLabel(rowIdx) = low.FactorLabel;
        Trajectory(rowIdx) = trajectory;
        ProbabilityEffect_pp(rowIdx) = 100 * ...
            (high.DeltaP_RingLine - low.DeltaP_RingLine);
        WarningEffect_s(rowIdx) = high.DeltaWarning_LineMinusRing_s - ...
            low.DeltaWarning_LineMinusRing_s;
        RingPointProbabilityEffect_pp(rowIdx) = 100 * ...
            (high.DeltaP_RingPoint - low.DeltaP_RingPoint);
        RingPointWarningEffect_s(rowIdx) = high.DeltaWarning_PointMinusRing_s - ...
            low.DeltaWarning_PointMinusRing_s;
        [probabilityCI, warningCI] = pairedEndpointEffectBootstrap( ...
            rawResults, factorName, trajectory, cfg.bootstrapN, ...
            cfg.seed + 20000 + rowIdx);
        ProbabilityEffect_CI_Low_pp(rowIdx) = 100 * probabilityCI(1);
        ProbabilityEffect_CI_High_pp(rowIdx) = 100 * probabilityCI(2);
        WarningEffect_CI_Low_s(rowIdx) = warningCI(1);
        WarningEffect_CI_High_s(rowIdx) = warningCI(2);
    end
end

effects = table(Factor, FactorLabel, Trajectory, ProbabilityEffect_pp, ...
    ProbabilityEffect_CI_Low_pp, ProbabilityEffect_CI_High_pp, ...
    WarningEffect_s, WarningEffect_CI_Low_s, WarningEffect_CI_High_s, ...
    RingPointProbabilityEffect_pp, RingPointWarningEffect_s);
end

function [probabilityCI, warningCI] = pairedEndpointEffectBootstrap( ...
        rawResults, factorName, trajectory, B, seed)
if B <= 0
    probabilityCI = [nan, nan];
    warningCI = [nan, nan];
    return;
end

factors = string({rawResults.factor});
paths = string({rawResults.trajectory});
levels = string({rawResults.level});
beams = lower(string({rawResults.beam}));
baseMask = factors == factorName & paths == trajectory;
lowRing = rawResults(find(baseMask & levels == "Low" & beams == "ring", 1));
lowLine = rawResults(find(baseMask & levels == "Low" & beams == "line", 1));
highRing = rawResults(find(baseMask & levels == "High" & beams == "ring", 1));
highLine = rawResults(find(baseMask & levels == "High" & beams == "line", 1));
if isempty(lowRing) || isempty(lowLine) || isempty(highRing) || isempty(highLine)
    error('sens_analyze_results:MissingEffectPair', ...
        'Factor %s trajectory %s lacks paired low/high ring-line results.', ...
        factorName, trajectory);
end

scenarioId = lowRing.scenarioId(:);
if ~isequal(scenarioId, lowLine.scenarioId(:), highRing.scenarioId(:), ...
        highLine.scenarioId(:))
    error('sens_analyze_results:UnpairedFactorEndpoints', ...
        'Low and high endpoints must share identical scenario IDs.');
end

N = numel(scenarioId);
stream = RandStream('mt19937ar', 'Seed', seed);
probabilitySamples = nan(B, 1);
warningSamples = nan(B, 1);
lowRingDetect = double(logical(lowRing.detect(:)));
lowLineDetect = double(logical(lowLine.detect(:)));
highRingDetect = double(logical(highRing.detect(:)));
highLineDetect = double(logical(highLine.detect(:)));
lowRingTime = maskedTimes(lowRing);
lowLineTime = maskedTimes(lowLine);
highRingTime = maskedTimes(highRing);
highLineTime = maskedTimes(highLine);

chunkSize = 100;
for chunkStart = 1:chunkSize:B
    chunkEnd = min(chunkStart + chunkSize - 1, B);
    idx = randi(stream, N, N, chunkEnd - chunkStart + 1);
    highProbability = mean(highRingDetect(idx) - highLineDetect(idx), 1);
    lowProbability = mean(lowRingDetect(idx) - lowLineDetect(idx), 1);
    probabilitySamples(chunkStart:chunkEnd) = ...
        (highProbability - lowProbability).';
    highWarning = mean(highLineTime(idx), 1, 'omitnan') - ...
        mean(highRingTime(idx), 1, 'omitnan');
    lowWarning = mean(lowLineTime(idx), 1, 'omitnan') - ...
        mean(lowRingTime(idx), 1, 'omitnan');
    warningSamples(chunkStart:chunkEnd) = (highWarning - lowWarning).';
end
probabilityCI = empiricalInterval(probabilitySamples);
warningCI = empiricalInterval(warningSamples);
end

function times = maskedTimes(item)
times = item.firstTime_s(:);
times(~logical(item.detect(:))) = nan;
end

function verifyPairing(ring, line, point)
ringId = ring.scenarioId(:);
if ~isequal(ringId, line.scenarioId(:)) || ~isequal(ringId, point.scenarioId(:))
    error('sens_analyze_results:UnpairedScenarios', ...
        'Point, line, and ring results do not share identical scenario IDs.');
end
end

function value = conditionalMean(item)
detected = logical(item.detect(:));
times = item.firstTime_s(:);
valid = detected & isfinite(times);
if any(valid)
    value = mean(times(valid));
else
    value = nan;
end
end

function [probabilityCI, warningCI] = pairedBootstrap(first, second, B, seed)
if B <= 0
    probabilityCI = [nan, nan];
    warningCI = [nan, nan];
    return;
end
stream = RandStream('mt19937ar', 'Seed', seed);
N = numel(first.detect);
probabilitySamples = nan(B, 1);
warningSamples = nan(B, 1);
firstDetect = double(logical(first.detect(:)));
secondDetect = double(logical(second.detect(:)));
firstTime = first.firstTime_s(:);
secondTime = second.firstTime_s(:);
firstTime(~logical(first.detect(:))) = nan;
secondTime(~logical(second.detect(:))) = nan;

chunkSize = 100;
for chunkStart = 1:chunkSize:B
    chunkEnd = min(chunkStart + chunkSize - 1, B);
    idx = randi(stream, N, N, chunkEnd - chunkStart + 1);
    probabilitySamples(chunkStart:chunkEnd) = ...
        mean(firstDetect(idx) - secondDetect(idx), 1).';
    warningSamples(chunkStart:chunkEnd) = ...
        (mean(secondTime(idx), 1, 'omitnan') - ...
        mean(firstTime(idx), 1, 'omitnan')).';
end
probabilityCI = empiricalInterval(probabilitySamples);
warningCI = empiricalInterval(warningSamples);
end

function [low, high] = bootstrapMeanInterval(values, B, seed)
if B <= 0
    low = nan;
    high = nan;
    return;
end
values = values(:);
N = numel(values);
stream = RandStream('mt19937ar', 'Seed', seed);
samples = nan(B, 1);
chunkSize = 100;
for chunkStart = 1:chunkSize:B
    chunkEnd = min(chunkStart + chunkSize - 1, B);
    idx = randi(stream, N, N, chunkEnd - chunkStart + 1);
    samples(chunkStart:chunkEnd) = mean(values(idx), 1).';
end
interval = empiricalInterval(samples);
low = interval(1);
high = interval(2);
end

function interval = empiricalInterval(samples)
samples = sort(samples(isfinite(samples)));
if isempty(samples)
    interval = [nan, nan];
    return;
end
n = numel(samples);
lowIdx = max(1, ceil(0.025 * n));
highIdx = max(1, min(n, floor(0.975 * n)));
interval = [samples(lowIdx), samples(highIdx)];
end

function [low, high] = wilsonInterval(successes, trials)
z = 1.95996398454005;
p = successes / trials;
denominator = 1 + z^2 / trials;
center = (p + z^2 / (2 * trials)) / denominator;
halfWidth = z / denominator * sqrt(p * (1 - p) / trials + ...
    z^2 / (4 * trials^2));
low = max(0, center - halfWidth);
high = min(1, center + halfWidth);
end

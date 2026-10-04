function [summary, timelySummary, isoPerformance, effects] = ...
        ea_sens_analyze_results(rawResults, cfg)
%EA_SENS_ANALYZE_RESULTS Summarize paired ROE sensitivity simulations.
%
% The primary robustness response is
%   kappa = minimum annular wD / minimum line wD
% for a common time-limited detection-probability requirement. kappa < 1
% means that the annular beam reaches the same performance with a smaller
% maximum off-axis field requirement.

arguments
    rawResults (1,:) struct
    cfg (1,1) struct
end

required = {"caseId", "factor", "factorLabel", "level", "levelValue", ...
    "unit", "beam", "wD_rad", "scenarioId", "detect", "firstTime_s"};
if isempty(rawResults) || ~all(isfield(rawResults, cellstr(required)))
    error("ea_sens_analyze_results:InvalidRawResults", ...
        "rawResults is empty or missing required fields.");
end

summary = buildSummary(rawResults, cfg);
timelySummary = buildTimelySummary(rawResults, cfg);
isoPerformance = buildIsoPerformance(rawResults, cfg);
effects = buildEffects(isoPerformance, cfg);
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
WD_mrad = nan(nRows, 1);
Wd_mrad = nan(nRows, 1);
N = zeros(nRows, 1);
NDetected = zeros(nRows, 1);
PDetect = nan(nRows, 1);
PDetect_CI_Low = nan(nRows, 1);
PDetect_CI_High = nan(nRows, 1);
MeanWarning_s = nan(nRows, 1);
MedianWarning_s = nan(nRows, 1);
MeanEncounter_s = nan(nRows, 1);
MeanEffectivePulses = nan(nRows, 1);
PWithinPrimaryTime = nan(nRows, 1);

for rowIdx = 1:nRows
    item = rawResults(rowIdx);
    detected = logical(item.detect(:));
    firstTime = item.firstTime_s(:);
    validWarning = detected & isfinite(firstTime);
    sampleCount = numel(detected);
    detectedCount = nnz(detected);

    CaseId(rowIdx) = string(item.caseId);
    Factor(rowIdx) = string(item.factor);
    FactorLabel(rowIdx) = string(item.factorLabel);
    Level(rowIdx) = string(item.level);
    LevelValue(rowIdx) = item.levelValue;
    Unit(rowIdx) = string(item.unit);
    Beam(rowIdx) = string(item.beam);
    WD_mrad(rowIdx) = 1e3 * item.wD_rad;
    Wd_mrad(rowIdx) = 1e3 * item.wd_rad;
    N(rowIdx) = sampleCount;
    NDetected(rowIdx) = detectedCount;
    PDetect(rowIdx) = detectedCount / sampleCount;
    [PDetect_CI_Low(rowIdx), PDetect_CI_High(rowIdx)] = ...
        wilsonInterval(detectedCount, sampleCount);
    if any(validWarning)
        MeanWarning_s(rowIdx) = mean(firstTime(validWarning));
        MedianWarning_s(rowIdx) = median(firstTime(validWarning));
    end
    if isfield(item, "firstEncounter_s")
        encounter = item.firstEncounter_s(:);
        MeanEncounter_s(rowIdx) = mean(encounter(isfinite(encounter)));
    end
    if isfield(item, "effectivePulses")
        pulses = item.effectivePulses(:);
        MeanEffectivePulses(rowIdx) = mean(pulses(isfinite(pulses)));
    end
    PWithinPrimaryTime(rowIdx) = mean( ...
        validWarning & firstTime <= cfg.primaryTimelyThreshold_s);
end

summary = table(CaseId, Factor, FactorLabel, Level, LevelValue, Unit, ...
    Beam, WD_mrad, Wd_mrad, N, NDetected, PDetect, ...
    PDetect_CI_Low, PDetect_CI_High, MeanWarning_s, MedianWarning_s, ...
    MeanEncounter_s, MeanEffectivePulses, PWithinPrimaryTime);
end

function timely = buildTimelySummary(rawResults, cfg)
thresholds = cfg.timelyThresholds_s(:);
nRows = numel(rawResults) * numel(thresholds);
CaseId = strings(nRows, 1);
Factor = strings(nRows, 1);
FactorLabel = strings(nRows, 1);
Level = strings(nRows, 1);
LevelValue = nan(nRows, 1);
Unit = strings(nRows, 1);
Beam = strings(nRows, 1);
WD_mrad = nan(nRows, 1);
TimeThreshold_s = nan(nRows, 1);
N = zeros(nRows, 1);
NTimely = zeros(nRows, 1);
PTimely = nan(nRows, 1);
PTimely_CI_Low = nan(nRows, 1);
PTimely_CI_High = nan(nRows, 1);

rowIdx = 0;
for itemIdx = 1:numel(rawResults)
    item = rawResults(itemIdx);
    detected = logical(item.detect(:));
    firstTime = item.firstTime_s(:);
    sampleCount = numel(detected);
    for thresholdIdx = 1:numel(thresholds)
        rowIdx = rowIdx + 1;
        timelyMask = detected & isfinite(firstTime) & ...
            firstTime <= thresholds(thresholdIdx);
        timelyCount = nnz(timelyMask);
        CaseId(rowIdx) = string(item.caseId);
        Factor(rowIdx) = string(item.factor);
        FactorLabel(rowIdx) = string(item.factorLabel);
        Level(rowIdx) = string(item.level);
        LevelValue(rowIdx) = item.levelValue;
        Unit(rowIdx) = string(item.unit);
        Beam(rowIdx) = string(item.beam);
        WD_mrad(rowIdx) = 1e3 * item.wD_rad;
        TimeThreshold_s(rowIdx) = thresholds(thresholdIdx);
        N(rowIdx) = sampleCount;
        NTimely(rowIdx) = timelyCount;
        PTimely(rowIdx) = timelyCount / sampleCount;
        [PTimely_CI_Low(rowIdx), PTimely_CI_High(rowIdx)] = ...
            wilsonInterval(timelyCount, sampleCount);
    end
end

timely = table(CaseId, Factor, FactorLabel, Level, LevelValue, Unit, ...
    Beam, WD_mrad, TimeThreshold_s, N, NTimely, PTimely, ...
    PTimely_CI_Low, PTimely_CI_High);
end

function iso = buildIsoPerformance(rawResults, cfg)
caseIds = unique(string({rawResults.caseId}), "stable");
thresholds = cfg.timelyThresholds_s(:);
probabilities = cfg.targetProbabilities(:);
nRows = numel(caseIds) * numel(thresholds) * numel(probabilities);

CaseId = strings(nRows, 1);
Factor = strings(nRows, 1);
FactorLabel = strings(nRows, 1);
Level = strings(nRows, 1);
LevelValue = nan(nRows, 1);
Unit = strings(nRows, 1);
TimeThreshold_s = nan(nRows, 1);
TargetProbability = nan(nRows, 1);
RingMinWD_mrad = nan(nRows, 1);
LineMinWD_mrad = nan(nRows, 1);
PointMinWD_mrad = nan(nRows, 1);
KappaRingLine = nan(nRows, 1);
Kappa_CI_Low = nan(nRows, 1);
Kappa_CI_High = nan(nRows, 1);
AnnularAdvantage = false(nRows, 1);

rowIdx = 0;
for caseIdx = 1:numel(caseIds)
    members = rawResults(string({rawResults.caseId}) == caseIds(caseIdx));
    verifyPairing(members);
    item0 = members(1);
    for thresholdIdx = 1:numel(thresholds)
        for probabilityIdx = 1:numel(probabilities)
            rowIdx = rowIdx + 1;
            threshold = thresholds(thresholdIdx);
            targetProbability = probabilities(probabilityIdx);
            minima = minimumDivergences(members, threshold, ...
                targetProbability, []);
            kappa = minima.ring / minima.line;
            [ciLow, ciHigh] = bootstrapKappa(members, threshold, ...
                targetProbability, cfg.bootstrapN, ...
                cfg.seed + 10000 * caseIdx + rowIdx);

            CaseId(rowIdx) = caseIds(caseIdx);
            Factor(rowIdx) = string(item0.factor);
            FactorLabel(rowIdx) = string(item0.factorLabel);
            Level(rowIdx) = string(item0.level);
            LevelValue(rowIdx) = item0.levelValue;
            Unit(rowIdx) = string(item0.unit);
            TimeThreshold_s(rowIdx) = threshold;
            TargetProbability(rowIdx) = targetProbability;
            RingMinWD_mrad(rowIdx) = 1e3 * minima.ring;
            LineMinWD_mrad(rowIdx) = 1e3 * minima.line;
            PointMinWD_mrad(rowIdx) = 1e3 * minima.point;
            KappaRingLine(rowIdx) = kappa;
            Kappa_CI_Low(rowIdx) = ciLow;
            Kappa_CI_High(rowIdx) = ciHigh;
            AnnularAdvantage(rowIdx) = isfinite(kappa) && kappa < 1;
        end
    end
end

iso = table(CaseId, Factor, FactorLabel, Level, LevelValue, Unit, ...
    TimeThreshold_s, TargetProbability, RingMinWD_mrad, ...
    LineMinWD_mrad, PointMinWD_mrad, KappaRingLine, ...
    Kappa_CI_Low, Kappa_CI_High, AnnularAdvantage);
end

function effects = buildEffects(iso, cfg)
primary = iso(abs(iso.TimeThreshold_s - ...
    cfg.primaryTimelyThreshold_s) < 1e-12 & ...
    abs(iso.TargetProbability - cfg.primaryProbability) < 1e-12, :);
factorNames = unique(primary.Factor(primary.Factor ~= "baseline"), "stable");

Factor = strings(numel(factorNames), 1);
FactorLabel = strings(numel(factorNames), 1);
LowKappa = nan(numel(factorNames), 1);
HighKappa = nan(numel(factorNames), 1);
Effect = nan(numel(factorNames), 1);
AbsEffect = nan(numel(factorNames), 1);
RelativeEffect_pct = nan(numel(factorNames), 1);
RobustAtBothEndpoints = false(numel(factorNames), 1);

writeIdx = 0;
for factorIdx = 1:numel(factorNames)
    low = primary(primary.Factor == factorNames(factorIdx) & ...
        primary.Level == "Low", :);
    high = primary(primary.Factor == factorNames(factorIdx) & ...
        primary.Level == "High", :);
    if height(low) ~= 1 || height(high) ~= 1
        continue;
    end
    writeIdx = writeIdx + 1;
    Factor(writeIdx) = factorNames(factorIdx);
    FactorLabel(writeIdx) = low.FactorLabel;
    LowKappa(writeIdx) = low.KappaRingLine;
    HighKappa(writeIdx) = high.KappaRingLine;
    Effect(writeIdx) = HighKappa(writeIdx) - LowKappa(writeIdx);
    AbsEffect(writeIdx) = abs(Effect(writeIdx));
    RelativeEffect_pct(writeIdx) = 100 * Effect(writeIdx) / ...
        LowKappa(writeIdx);
    RobustAtBothEndpoints(writeIdx) = ...
        isfinite(LowKappa(writeIdx)) && isfinite(HighKappa(writeIdx)) && ...
        LowKappa(writeIdx) < 1 && HighKappa(writeIdx) < 1;
end

Factor = Factor(1:writeIdx);
FactorLabel = FactorLabel(1:writeIdx);
LowKappa = LowKappa(1:writeIdx);
HighKappa = HighKappa(1:writeIdx);
Effect = Effect(1:writeIdx);
AbsEffect = AbsEffect(1:writeIdx);
RelativeEffect_pct = RelativeEffect_pct(1:writeIdx);
RobustAtBothEndpoints = RobustAtBothEndpoints(1:writeIdx);
effects = table(Factor, FactorLabel, LowKappa, HighKappa, Effect, ...
    AbsEffect, RelativeEffect_pct, RobustAtBothEndpoints);
if ~isempty(effects)
    effects = sortrows(effects, "AbsEffect", "descend");
end
end

function minima = minimumDivergences(items, threshold, targetProbability, idx)
beamNames = ["ring", "line", "point"];
values = nan(size(beamNames));
for beamIdx = 1:numel(beamNames)
    beamItems = items(strcmpi(string({items.beam}), beamNames(beamIdx)));
    if isempty(beamItems)
        continue;
    end
    wD = [beamItems.wD_rad].';
    probability = nan(size(wD));
    for itemIdx = 1:numel(beamItems)
        detected = logical(beamItems(itemIdx).detect(:));
        firstTime = beamItems(itemIdx).firstTime_s(:);
        if isempty(idx)
            selected = 1:numel(detected);
        else
            selected = idx;
        end
        probability(itemIdx) = mean(detected(selected) & ...
            isfinite(firstTime(selected)) & firstTime(selected) <= threshold);
    end
    [wD, order] = sort(wD);
    probability = probability(order);
    values(beamIdx) = firstMonotoneCrossing(wD, probability, ...
        targetProbability);
end
minima = struct("ring", values(1), "line", values(2), "point", values(3));
end

function wDMin = firstMonotoneCrossing(wD, probability, target)
probability = cummax(probability(:));
wD = wD(:);
crossingIdx = find(probability >= target, 1, "first");
if isempty(crossingIdx)
    wDMin = nan;
elseif crossingIdx == 1
    wDMin = wD(1);
else
    x = probability(crossingIdx-1:crossingIdx);
    y = wD(crossingIdx-1:crossingIdx);
    if abs(diff(x)) < eps(max(abs(x)))
        wDMin = y(2);
    else
        wDMin = interp1(x, y, target, "linear");
    end
end
end

function [low, high] = bootstrapKappa(items, threshold, target, B, seed)
if B <= 0
    low = nan;
    high = nan;
    return;
end
N = numel(items(1).scenarioId);
stream = RandStream("mt19937ar", "Seed", seed);
samples = nan(B, 1);
for bootstrapIdx = 1:B
    idx = randi(stream, N, N, 1);
    minima = minimumDivergences(items, threshold, target, idx);
    samples(bootstrapIdx) = minima.ring / minima.line;
end
interval = empiricalInterval(samples);
low = interval(1);
high = interval(2);
end

function verifyPairing(items)
reference = items(1).scenarioId(:);
for itemIdx = 2:numel(items)
    if ~isequal(reference, items(itemIdx).scenarioId(:))
        error("ea_sens_analyze_results:UnpairedScenarios", ...
            "All beams and divergence levels in a case must share IDs.");
    end
end
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

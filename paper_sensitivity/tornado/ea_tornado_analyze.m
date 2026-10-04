function analysis = ea_tornado_analyze(outputDirectory, replicates, seed)
%EA_TORNADO_ANALYZE Analyze completed paired EA tornado simulations.

arguments
    outputDirectory (1,1) string
    replicates (1,1) double {mustBeInteger,mustBePositive} = 1000
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 20260929
end

designFile = fullfile(outputDirectory, "design.mat");
assert(isfile(designFile), "EATornado:MissingDesign", ...
    "Run simulation before analysis: %s", designFile);
saved = load(designFile, "design");
design = saved.design;
cfg = design.config;
assert(all(design.complete(design.required)), "EATornado:Incomplete", ...
    "Every required configuration/width checkpoint must be complete.");

weights = bootstrapWeights(cfg.N, replicates, seed);
nCases = numel(design.cases);
caseResults = repmat(emptyCaseResult(), nCases, 1);
curveRecords = cell(nCases, 1);
for caseIndex = 1:nCases
    [caseResults(caseIndex), curveRecords{caseIndex}] = analyzeCase( ...
        outputDirectory, design, caseIndex, weights);
end
curveTable = vertcat(curveRecords{:});

baselineKappa = caseResults(1).kappa;
baselineBootstrap = caseResults(1).bootstrapKappa;
levelRecords = cell(numel(cfg.factorNames) * numel(cfg.factorScales), 1);
writeIndex = 0;
for factorIndex = 1:numel(cfg.factorNames)
    for scale = cfg.factorScales
        writeIndex = writeIndex + 1;
        if abs(scale - 1) < 1e-12
            caseIndex = 1;
        else
            caseIndex = find([design.cases.factorIndex] == factorIndex & ...
                abs([design.cases.scale] - scale) < 1e-12, 1);
        end
        item = design.cases(caseIndex);
        result = caseResults(caseIndex);
        differenceBootstrap = result.bootstrapKappa - baselineBootstrap;
        differenceFinite = isfinite(differenceBootstrap);
        differenceCI = [nan, nan];
        if mean(differenceFinite) >= cfg.feasibilityThreshold && ...
                isfinite(result.kappa) && isfinite(baselineKappa)
            differenceCI = empiricalInterval( ...
                differenceBootstrap(differenceFinite));
        end
        levelRecords{writeIndex} = struct( ...
            "Factor", cfg.factorNames(factorIndex), ...
            "FactorLabel", cfg.factorLabels(factorIndex), ...
            "Scale", scale, "ChangePercent", 100 * (scale - 1), ...
            "Value", item.values(factorIndex) * ...
                cfg.displayScale(factorIndex), ...
            "Unit", cfg.factorUnits(factorIndex), ...
            "CaseIndex", caseIndex, "CaseId", item.id, ...
            "LineMinWD_mrad", result.minimumWidth(1), ...
            "LineMinWD_CI_Low", result.widthCI(1, 1), ...
            "LineMinWD_CI_High", result.widthCI(1, 2), ...
            "RingMinWD_mrad", result.minimumWidth(2), ...
            "RingMinWD_CI_Low", result.widthCI(2, 1), ...
            "RingMinWD_CI_High", result.widthCI(2, 2), ...
            "LineGridWD_mrad", result.gridWidth(1), ...
            "RingGridWD_mrad", result.gridWidth(2), ...
            "LineBracketLow_mrad", result.bracket(1, 1), ...
            "LineBracketHigh_mrad", result.bracket(1, 2), ...
            "RingBracketLow_mrad", result.bracket(2, 1), ...
            "RingBracketHigh_mrad", result.bracket(2, 2), ...
            "Kappa", result.kappa, ...
            "KappaCI_Low", result.kappaCI(1), ...
            "KappaCI_High", result.kappaCI(2), ...
            "BothFeasibleBootstrap", result.bothFeasible, ...
            "DeltaKappa", result.kappa - baselineKappa, ...
            "DeltaKappaCI_Low", differenceCI(1), ...
            "DeltaKappaCI_High", differenceCI(2), ...
            "AnnularAdvantage", isfinite(result.kappa) && ...
                result.kappa < 1);
    end
end
levels = struct2table(vertcat(levelRecords{:}));
effects = ea_tornado_build_effects(levels);
physicalChecks = checkPhysicalDirections(design, caseResults);

writetable(levels, fullfile(outputDirectory, ...
    "ea_tornado_level_results.csv"));
writetable(effects, fullfile(outputDirectory, ...
    "ea_tornado_effects.csv"));
writetable(curveTable, fullfile(outputDirectory, ...
    "ea_tornado_probability_curves.csv"));
writetable(physicalChecks, fullfile(outputDirectory, ...
    "physical_direction_checks.csv"));

analysis = struct();
analysis.schemaVersion = 1;
analysis.createdAt = string(datetime("now", ...
    "Format", "yyyy-MM-dd HH:mm:ss Z"));
analysis.options = struct("Replicates", replicates, "Seed", seed, ...
    "FeasibilityThreshold", cfg.feasibilityThreshold);
analysis.design = design;
analysis.levels = levels;
analysis.effects = effects;
analysis.curves = curveTable;
analysis.physicalChecks = physicalChecks;
analysis.baselineKappa = baselineKappa;
analysis.baselineKappaCI = caseResults(1).kappaCI;
compactResults = rmfield(caseResults, 'events');
analysis.caseResults = rmfield(compactResults, 'bootstrapKappa');
save(fullfile(outputDirectory, "ea_tornado_analysis.mat"), ...
    "analysis", "-v7.3");
end

function [result, curveTable] = analyzeCase( ...
        outputDirectory, design, caseIndex, weights)
cfg = design.config;
widthMask = design.required(:, caseIndex);
widths = design.widthsMrad(widthMask).';
nWidths = numel(widths);
events = false(cfg.N, nWidths, 2);
counts = zeros(nWidths, 2);
for widthIndex = 1:nWidths
    rawFile = fullfile(outputDirectory, "raw", ...
        design.cases(caseIndex).id, sprintf( ...
        "ea_tornado_width_%03d.mat", round(widths(widthIndex))));
    raw = load(rawFile, "meta", "hit", "time_s");
    validateRaw(raw, design, caseIndex, widths(widthIndex));
    oneEvent = raw.hit & raw.time_s <= cfg.deadline_s;
    events(:, widthIndex, :) = reshape(oneEvent, cfg.N, 1, 2);
    counts(widthIndex, :) = sum(oneEvent, 1);
end
probability = counts / cfg.N;
[minimumWidth, gridWidth, bracket] = ea_tornado_min_width( ...
    widths, probability, cfg.requiredProbability);

replicates = size(weights, 2);
flatEvents = reshape(events, cfg.N, []);
bootstrapProbability = flatEvents.' * weights / cfg.N;
bootstrapProbability = reshape(bootstrapProbability, nWidths, 2, replicates);
bootstrapWidth = nan(replicates, 2);
for replicate = 1:replicates
    bootstrapWidth(replicate, :) = ea_tornado_min_width(widths, ...
        bootstrapProbability(:, :, replicate), cfg.requiredProbability);
end
bootstrapKappa = bootstrapWidth(:, 2) ./ bootstrapWidth(:, 1);
bothFinite = all(isfinite(bootstrapWidth), 2);
bothFeasible = mean(bothFinite);
kappa = minimumWidth(2) / minimumWidth(1);
kappaCI = [nan, nan];
widthCI = nan(2, 2);
if bothFeasible >= cfg.feasibilityThreshold && isfinite(kappa)
    kappaCI = empiricalInterval(bootstrapKappa(bothFinite));
end
for beamIndex = 1:2
    finiteWidth = isfinite(bootstrapWidth(:, beamIndex));
    if mean(finiteWidth) >= cfg.feasibilityThreshold && ...
            isfinite(minimumWidth(beamIndex))
        widthCI(beamIndex, :) = empiricalInterval( ...
            bootstrapWidth(finiteWidth, beamIndex));
    end
end

item = design.cases(caseIndex);
curveTable = table();
curveTable.CaseIndex = repmat(caseIndex, nWidths, 1);
curveTable.CaseId = repmat(item.id, nWidths, 1);
curveTable.FactorIndex = repmat(item.factorIndex, nWidths, 1);
curveTable.Scale = repmat(item.scale, nWidths, 1);
curveTable.WD_mrad = widths;
curveTable.LineTimelyN = counts(:, 1);
curveTable.RingTimelyN = counts(:, 2);
curveTable.LinePTimely = probability(:, 1);
curveTable.RingPTimely = probability(:, 2);

result = emptyCaseResult();
result.caseIndex = caseIndex;
result.caseId = item.id;
result.widthsMrad = widths;
result.probability = probability;
result.minimumWidth = minimumWidth;
result.gridWidth = gridWidth;
result.bracket = bracket;
result.widthCI = widthCI;
result.kappa = kappa;
result.kappaCI = kappaCI;
result.bothFeasible = bothFeasible;
result.bootstrapKappa = bootstrapKappa;
result.events = events;
end

function checks = checkPhysicalDirections(design, results)
cfg = design.config;
factorIndices = [1, 2, 4];
expectedDirection = ["nondecreasing", "nondecreasing", "nonincreasing"];
records = cell(numel(factorIndices), 1);
for checkIndex = 1:numel(factorIndices)
    factorIndex = factorIndices(checkIndex);
    caseIndices = zeros(size(cfg.factorScales));
    for levelIndex = 1:numel(cfg.factorScales)
        scale = cfg.factorScales(levelIndex);
        if abs(scale - 1) < 1e-12
            caseIndices(levelIndex) = 1;
        else
            caseIndices(levelIndex) = find( ...
                [design.cases.factorIndex] == factorIndex & ...
                abs([design.cases.scale] - scale) < 1e-12, 1);
        end
    end
    commonWidths = results(caseIndices(1)).widthsMrad;
    for levelIndex = 2:numel(caseIndices)
        commonWidths = intersect(commonWidths, ...
            results(caseIndices(levelIndex)).widthsMrad, "stable");
    end
    violations = 0;
    comparisons = 0;
    for levelIndex = 1:(numel(caseIndices) - 1)
        low = selectEvents(results(caseIndices(levelIndex)), commonWidths);
        high = selectEvents(results(caseIndices(levelIndex + 1)), commonWidths);
        if factorIndex == 4
            violations = violations + nnz(high & ~low);
        else
            violations = violations + nnz(low & ~high);
        end
        comparisons = comparisons + numel(low);
    end
    records{checkIndex} = struct( ...
        "Factor", cfg.factorNames(factorIndex), ...
        "FactorLabel", cfg.factorLabels(factorIndex), ...
        "ExpectedDirection", expectedDirection(checkIndex), ...
        "ComparedEvents", comparisons, "Violations", violations, ...
        "Passed", violations == 0);
end
checks = struct2table(vertcat(records{:}));
end

function selected = selectEvents(result, widths)
[found, locations] = ismember(widths, result.widthsMrad);
assert(all(found), "EATornado:DirectionGrid", ...
    "Direction-check widths are not shared.");
selected = result.events(:, locations, :);
end

function validateRaw(raw, design, caseIndex, widthMrad)
cfg = design.config;
item = design.cases(caseIndex);
assert(raw.meta.CaseIndex == caseIndex && raw.meta.CaseId == item.id && ...
    raw.meta.WidthMrad == widthMrad && raw.meta.N == cfg.N && ...
    raw.meta.Seed == cfg.seed && isequaln(raw.meta.Physics, item.physics) && ...
    raw.meta.Speed_m_s == item.speed_m_s && ...
    isequaln(raw.meta.Signature, design.signature) && ...
    isequal(size(raw.hit), [cfg.N, 2]) && ...
    isequal(size(raw.time_s), [cfg.N, 2]), ...
    "EATornado:RawSchema", "Invalid or incompatible raw result.");
end

function weights = bootstrapWeights(n, replicates, seed)
stream = RandStream("mt19937ar", "Seed", seed);
weights = zeros(n, replicates);
edges = 0.5:1:(n + 0.5);
for replicate = 1:replicates
    sample = randi(stream, n, n, 1);
    weights(:, replicate) = histcounts(sample, edges).';
end
end

function interval = empiricalInterval(values)
values = sort(values(isfinite(values)));
if isempty(values)
    interval = [nan, nan];
    return;
end
n = numel(values);
lowIndex = max(1, ceil(0.025 * n));
highIndex = min(n, max(1, floor(0.975 * n)));
interval = [values(lowIndex), values(highIndex)];
end

function result = emptyCaseResult()
result = struct("caseIndex", nan, "caseId", "", ...
    "widthsMrad", zeros(0, 1), "probability", zeros(0, 2), ...
    "minimumWidth", [nan, nan], "gridWidth", [nan, nan], ...
    "bracket", nan(2, 2), "widthCI", nan(2, 2), ...
    "kappa", nan, "kappaCI", [nan, nan], ...
    "bothFeasible", nan, "bootstrapKappa", zeros(0, 1), ...
    "events", false(0, 0, 2));
end

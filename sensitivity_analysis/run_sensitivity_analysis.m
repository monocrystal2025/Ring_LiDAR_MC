function outputs = run_sensitivity_analysis(varargin)
%RUN_SENSITIVITY_ANALYSIS Execute the standalone paired ROI sensitivity study.
%
% Examples
%   run_sensitivity_analysis
%   run_sensitivity_analysis('Mode','smoke','UseParallel',false)
%   run_sensitivity_analysis('Mode','screen')
%   run_sensitivity_analysis('Mode','plotonly')
%
% The default full mode performs screening, selects two influential factors,
% runs five-level refinement, writes CSV/MAT files, and creates the figure.

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'Mode', "full", @(x) any(strcmpi(string(x), ...
    ["full", "screen", "smoke", "plotonly"])));
addParameter(parser, 'OutputRoot', "", @(x) isstring(x) || ischar(x));
addParameter(parser, 'ScreenN', [], @(x) isempty(x) || (isscalar(x) && x > 0));
addParameter(parser, 'RefineN', [], @(x) isempty(x) || (isscalar(x) && x > 0));
addParameter(parser, 'BootstrapN', [], @(x) isempty(x) || (isscalar(x) && x >= 0));
addParameter(parser, 'UseParallel', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
addParameter(parser, 'Force', false, @(x) islogical(x) || isnumeric(x));
parse(parser, varargin{:});
options = parser.Results;
mode = lower(string(options.Mode));

cfg = sens_default_config();
cfg = applyOptions(cfg, options, mode);
ensureOutputFolders(cfg);

if mode == "plotonly"
    figureHandles = plot_sensitivity_results(cfg.outputRoot);
    outputs = struct('config', cfg, 'figures', figureHandles);
    return;
end

screenBank = sens_build_scenario_bank(cfg.screenN, cfg.seed);
screenCases = buildScreenCases(cfg);
screenRaw = runCases(cfg, screenCases, screenBank, "screen", logical(options.Force));
[screenSummary, screenComparisons, screenEffects] = ...
    sens_analyze_results(screenRaw, cfg);

writeAnalysisTables(cfg.resultsDir, "sensitivity_screen", ...
    screenSummary, screenComparisons, screenEffects);
save(fullfile(cfg.resultsDir, 'sensitivity_screen_results.mat'), ...
    'cfg', 'screenRaw', 'screenSummary', 'screenComparisons', ...
    'screenEffects', '-v7.3');

if mode == "smoke"
    validation = table("Skipped for reduced smoke configuration", true, ...
        'VariableNames', {'Message', 'Pass'});
else
    validation = validateAgainstExisting(screenSummary, cfg);
end
writetable(validation, fullfile(cfg.resultsDir, 'baseline_validation.csv'));

[probabilityFactor, warningFactor] = selectRefinementFactors(screenEffects);
selectedFactors = table(probabilityFactor, warningFactor, ...
    'VariableNames', {'ProbabilityFactor', 'WarningTimeFactor'});
writetable(selectedFactors, fullfile(cfg.resultsDir, 'selected_factors.csv'));

outputs = struct();
outputs.config = cfg;
outputs.screenSummary = screenSummary;
outputs.screenComparisons = screenComparisons;
outputs.screenEffects = screenEffects;
outputs.validation = validation;
outputs.selectedFactors = selectedFactors;

if mode == "screen"
    return;
end

refineBank = sens_build_scenario_bank(cfg.refineN, cfg.seed + 1);
refineCases = buildRefinementCases(cfg, probabilityFactor, warningFactor);
refineRaw = runCases(cfg, refineCases, refineBank, "refine", logical(options.Force));
[refineSummary, refineComparisons] = sens_analyze_results(refineRaw, cfg);
directionValidation = validateEffectDirections(screenEffects, ...
    refineComparisons, probabilityFactor, warningFactor);
writetable(directionValidation, fullfile(cfg.resultsDir, ...
    'effect_direction_validation.csv'));
if mode == "smoke"
    refinedValidation = table("Skipped for reduced smoke configuration", true, ...
        'VariableNames', {'Message', 'Pass'});
else
    refinedValidation = validateRefinedBaseline( ...
        refineSummary, validation, cfg, probabilityFactor);
end
writetable(refinedValidation, fullfile(cfg.resultsDir, ...
    'baseline_validation_refined_N10000.csv'));
writeAnalysisTables(cfg.resultsDir, "sensitivity_refined", ...
    refineSummary, refineComparisons, table());
save(fullfile(cfg.resultsDir, 'sensitivity_refined_results.mat'), ...
    'cfg', 'refineRaw', 'refineSummary', 'refineComparisons', ...
    'selectedFactors', 'directionValidation', 'refinedValidation', '-v7.3');

outputs.refineSummary = refineSummary;
outputs.refineComparisons = refineComparisons;
outputs.directionValidation = directionValidation;
outputs.refinedValidation = refinedValidation;
outputs.figures = plot_sensitivity_results(cfg.outputRoot);
end

function validation = validateEffectDirections(screenEffects, ...
        refineComparisons, probabilityFactor, warningFactor)
metrics = ["Detection probability", "Warning time"];
factors = [probabilityFactor, warningFactor];
screenVariables = ["ProbabilityEffect_pp", "WarningEffect_s"];
refineVariables = ["DeltaP_RingLine", "DeltaWarning_LineMinusRing_s"];
trajectories = ["RA", "SP"];
nRows = numel(metrics) * numel(trajectories);
Metric = strings(nRows, 1);
Factor = strings(nRows, 1);
Trajectory = strings(nRows, 1);
ScreenEffect = nan(nRows, 1);
RefinedEffect = nan(nRows, 1);
DirectionConsistent = false(nRows, 1);

rowIdx = 0;
for metricIdx = 1:numel(metrics)
    for trajectoryIdx = 1:numel(trajectories)
        rowIdx = rowIdx + 1;
        factor = factors(metricIdx);
        trajectory = trajectories(trajectoryIdx);
        screen = screenEffects(screenEffects.Factor == factor & ...
            screenEffects.Trajectory == trajectory, :);
        refined = refineComparisons(refineComparisons.Factor == factor & ...
            refineComparisons.Trajectory == trajectory, :);
        refined = sortrows(refined, 'LevelValue');
        if height(screen) ~= 1 || height(refined) < 2
            error('run_sensitivity_analysis:MissingDirectionData', ...
                'Cannot validate %s for factor %s trajectory %s.', ...
                metrics(metricIdx), factor, trajectory);
        end
        Metric(rowIdx) = metrics(metricIdx);
        Factor(rowIdx) = factor;
        Trajectory(rowIdx) = trajectory;
        ScreenEffect(rowIdx) = screen.(screenVariables(metricIdx));
        RefinedEffect(rowIdx) = refined.(refineVariables(metricIdx))(end) - ...
            refined.(refineVariables(metricIdx))(1);
        if metricIdx == 1
            RefinedEffect(rowIdx) = 100 * RefinedEffect(rowIdx);
        end
        DirectionConsistent(rowIdx) = sign(ScreenEffect(rowIdx)) == ...
            sign(RefinedEffect(rowIdx));
    end
end
validation = table(Metric, Factor, Trajectory, ScreenEffect, ...
    RefinedEffect, DirectionConsistent);
end

function validation = validateRefinedBaseline( ...
        refineSummary, screenValidation, cfg, probabilityFactor)
required = {'Beam', 'Trajectory', 'ReferenceProbability', 'ReferenceWarning_s'};
if ~all(ismember(required, screenValidation.Properties.VariableNames))
    validation = table("Reference validation data unavailable", false, ...
        'VariableNames', {'Message', 'Pass'});
    return;
end

factor = cfg.factors([cfg.factors.name] == probabilityFactor);
baselineValue = factorBaselineValue(cfg, factor);
candidate = refineSummary(refineSummary.Factor == probabilityFactor, :);
[~, nearestIdx] = min(abs(candidate.LevelValue - baselineValue));
baselineLevel = candidate.LevelValue(nearestIdx);
candidate = candidate(candidate.LevelValue == baselineLevel, :);
nRows = height(screenValidation);
Beam = screenValidation.Beam;
Trajectory = screenValidation.Trajectory;
RefinedProbability = nan(nRows, 1);
ReferenceProbability = screenValidation.ReferenceProbability;
ProbabilityDifference = nan(nRows, 1);
RefinedWarning_s = nan(nRows, 1);
ReferenceWarning_s = screenValidation.ReferenceWarning_s;
WarningRelativeDifference = nan(nRows, 1);
Pass = false(nRows, 1);
for rowIdx = 1:nRows
    row = candidate(candidate.Beam == Beam(rowIdx) & ...
        candidate.Trajectory == Trajectory(rowIdx), :);
    if height(row) ~= 1
        error('run_sensitivity_analysis:MissingRefinedBaseline', ...
            'Missing refined baseline for %s %s.', Beam(rowIdx), Trajectory(rowIdx));
    end
    RefinedProbability(rowIdx) = row.PDetect;
    ProbabilityDifference(rowIdx) = RefinedProbability(rowIdx) - ...
        ReferenceProbability(rowIdx);
    RefinedWarning_s(rowIdx) = row.MeanWarning_s;
    WarningRelativeDifference(rowIdx) = abs(RefinedWarning_s(rowIdx) - ...
        ReferenceWarning_s(rowIdx)) / ReferenceWarning_s(rowIdx);
    Pass(rowIdx) = abs(ProbabilityDifference(rowIdx)) <= ...
        cfg.validation.probabilityAbsTolerance && ...
        WarningRelativeDifference(rowIdx) <= ...
        cfg.validation.warningTimeRelativeTolerance;
end
validation = table(Beam, Trajectory, RefinedProbability, ...
    ReferenceProbability, ProbabilityDifference, RefinedWarning_s, ...
    ReferenceWarning_s, WarningRelativeDifference, Pass);
end

function value = factorBaselineValue(cfg, factor)
if factor.key == "beamWidth_rad"
    beam = cfg.beams([cfg.beams.name] == "ring");
    value = beam.wd_rad;
elseif factor.key == "environmentSeverity"
    value = cfg.environment.nominalSeverity;
else
    value = cfg.system.(factor.key);
end
end

function cfg = applyOptions(cfg, options, mode)
if strlength(string(options.OutputRoot)) > 0
    cfg.outputRoot = char(string(options.OutputRoot));
    cfg.resultsDir = fullfile(cfg.outputRoot, 'results');
    cfg.figuresDir = fullfile(cfg.outputRoot, 'figures');
end
if ~isempty(options.ScreenN)
    cfg.screenN = options.ScreenN;
end
if ~isempty(options.RefineN)
    cfg.refineN = options.RefineN;
end
if ~isempty(options.BootstrapN)
    cfg.bootstrapN = options.BootstrapN;
end
if ~isempty(options.UseParallel)
    cfg.execution.useParallel = logical(options.UseParallel);
end
if mode == "smoke"
    cfg.screenN = min(cfg.screenN, 8);
    cfg.refineN = min(cfg.refineN, 10);
    cfg.bootstrapN = min(cfg.bootstrapN, 20);
    cfg.execution.useParallel = false;
    cfg.execution.blockSize = 2000;
    cfg.system.maxRange_m = 250;
    cfg.system.prf_Hz = 250;
    cfg.system.noiseLutStep_m = 2;
    cfg.factors = cfg.factors([2, 9]);
end
end

function ensureOutputFolders(cfg)
if ~isfolder(cfg.outputRoot)
    mkdir(cfg.outputRoot);
end
if ~isfolder(cfg.resultsDir)
    mkdir(cfg.resultsDir);
end
if ~isfolder(cfg.figuresDir)
    mkdir(cfg.figuresDir);
end
end

function cases = buildScreenCases(cfg)
baseline = makeCase("baseline", "baseline", "Baseline", "Baseline", ...
    0, "", cfg);
cases = baseline;
for factorIdx = 1:numel(cfg.factors)
    factor = cfg.factors(factorIdx);
    lowCfg = applyFactor(cfg, factor, factor.lowValue);
    highCfg = applyFactor(cfg, factor, factor.highValue);
    lowCase = makeCase(factor.name + "_low", factor.name, factor.label, ...
        "Low", factor.lowValue, factor.unit, lowCfg);
    highCase = makeCase(factor.name + "_high", factor.name, factor.label, ...
        "High", factor.highValue, factor.unit, highCfg);
    cases = [cases, lowCase, highCase]; %#ok<AGROW>
end
end

function cases = buildRefinementCases(cfg, probabilityFactor, warningFactor)
factorNames = unique([probabilityFactor, warningFactor], 'stable');
cases = repmat(makeCase("", "", "", "", nan, "", cfg), 1, 0);
for selectedIdx = 1:numel(factorNames)
    factor = cfg.factors(string({cfg.factors.name}) == factorNames(selectedIdx));
    values = linspace(factor.lowValue, factor.highValue, 5);
    for valueIdx = 1:numel(values)
        caseCfg = applyFactor(cfg, factor, values(valueIdx));
        caseId = "refine_" + factor.name + "_" + valueIdx;
        level = "Refine" + valueIdx;
        cases(end+1) = makeCase(caseId, factor.name, factor.label, ...
            level, values(valueIdx), factor.unit, caseCfg); %#ok<AGROW>
    end
end
end

function item = makeCase(caseId, factor, factorLabel, level, levelValue, unit, cfg)
item = struct( ...
    'caseId', string(caseId), ...
    'factor', string(factor), ...
    'factorLabel', string(factorLabel), ...
    'level', string(level), ...
    'levelValue', levelValue, ...
    'unit', string(unit), ...
    'config', cfg);
end

function caseCfg = applyFactor(cfg, factor, value)
caseCfg = cfg;
switch factor.key
    case "beamWidth_rad"
        for beamIdx = 1:numel(caseCfg.beams)
            if any(strcmpi(caseCfg.beams(beamIdx).name, ["line", "ring"]))
                caseCfg.beams(beamIdx).wd_rad = value;
            end
        end
    case "environmentSeverity"
        caseCfg.system = applyEnvironmentSeverity(caseCfg.system, cfg, value);
    otherwise
        caseCfg.system.(char(factor.key)) = value;
end
end

function system = applyEnvironmentSeverity(system, cfg, severity)
severity = max(-1, min(1, severity));
if severity <= 0
    t = severity + 1;
    extinctionScale = exp(interp1([0, 1], ...
        log([cfg.environment.clearExtinctionScale, 1]), t));
    backscatterScale = exp(interp1([0, 1], ...
        log([cfg.environment.clearBackscatterScale, 1]), t));
    radianceScale = exp(interp1([0, 1], ...
        log([cfg.environment.clearRadianceScale, 1]), t));
else
    t = severity;
    extinctionScale = exp(interp1([0, 1], ...
        log([1, cfg.environment.adverseExtinctionScale]), t));
    backscatterScale = exp(interp1([0, 1], ...
        log([1, cfg.environment.adverseBackscatterScale]), t));
    radianceScale = exp(interp1([0, 1], ...
        log([1, cfg.environment.adverseRadianceScale]), t));
end
system.extinction_m_inv = cfg.system.extinction_m_inv * extinctionScale;
system.backscatter_m_inv_sr_inv = cfg.system.backscatter_m_inv_sr_inv * backscatterScale;
system.skyRadiance_W_m2_sr_nm = cfg.system.skyRadiance_W_m2_sr_nm * radianceScale;
end

function raw = runCases(cfg, cases, bank, stage, force)
template = struct('caseId', "", 'factor', "", 'factorLabel', "", ...
    'level', "", 'levelValue', nan, 'unit', "", 'beam', "", ...
    'beamLabel', "", 'trajectory', "", 'scenarioId', [], ...
    'detect', [], 'firstTime_s', [], 'elapsed_s', nan, ...
    'pathPointCount', nan, 'cycleLength', nan, 'windowPulses', nan);
raw = repmat(template, 1, numel(cases) * numel(cfg.trajectories) * numel(cfg.beams));
writeIdx = 0;

for caseIdx = 1:numel(cases)
    caseInfo = cases(caseIdx);
    checkpoint = fullfile(cfg.resultsDir, sprintf('raw_%s_%s.mat', ...
        stage, char(caseInfo.caseId)));
    if isfile(checkpoint) && ~force
        loaded = load(checkpoint, 'caseRaw', 'bankSeed', 'bankN');
        if isfield(loaded, 'caseRaw') && isfield(loaded, 'bankSeed') && ...
                isfield(loaded, 'bankN') && loaded.bankSeed == bank.seed && ...
                loaded.bankN == bank.N
            count = numel(loaded.caseRaw);
            raw(writeIdx + (1:count)) = loaded.caseRaw;
            writeIdx = writeIdx + count;
            fprintf('Loaded checkpoint %s\n', checkpoint);
            continue;
        end
    end

    caseRaw = repmat(template, 1, numel(cfg.trajectories) * numel(cfg.beams));
    localIdx = 0;
    for trajectoryIdx = 1:numel(cfg.trajectories)
        trajectory = cfg.trajectories(trajectoryIdx);
        for beamIdx = 1:numel(caseInfo.config.beams)
            beam = caseInfo.config.beams(beamIdx);
            if cfg.execution.showProgress
                fprintf('[%s %d/%d] %s %s %s\n', upper(stage), caseIdx, ...
                    numel(cases), caseInfo.caseId, trajectory, beam.name);
            end
            sim = sens_simulate_roi_case(caseInfo.config, beam, trajectory, bank);
            localIdx = localIdx + 1;
            caseRaw(localIdx) = packageRaw(sim, caseInfo);
        end
    end
    bankSeed = bank.seed;
    bankN = bank.N;
    save(checkpoint, 'caseRaw', 'caseInfo', 'bankSeed', 'bankN', '-v7.3');
    raw(writeIdx + (1:numel(caseRaw))) = caseRaw;
    writeIdx = writeIdx + numel(caseRaw);
end
raw = raw(1:writeIdx);
end

function item = packageRaw(sim, caseInfo)
item = struct( ...
    'caseId', caseInfo.caseId, ...
    'factor', caseInfo.factor, ...
    'factorLabel', caseInfo.factorLabel, ...
    'level', caseInfo.level, ...
    'levelValue', caseInfo.levelValue, ...
    'unit', caseInfo.unit, ...
    'beam', sim.beam, ...
    'beamLabel', sim.beamLabel, ...
    'trajectory', sim.trajectory, ...
    'scenarioId', sim.scenarioId, ...
    'detect', sim.detect, ...
    'firstTime_s', sim.firstTime_s, ...
    'elapsed_s', sim.elapsed_s, ...
    'pathPointCount', sim.pathPointCount, ...
    'cycleLength', sim.cycleLength, ...
    'windowPulses', sim.windowPulses);
end

function writeAnalysisTables(resultsDir, stem, summary, comparisons, effects)
writetable(summary, fullfile(resultsDir, stem + "_summary.csv"));
writetable(comparisons, fullfile(resultsDir, stem + "_comparisons.csv"));
if ~isempty(effects)
    writetable(effects, fullfile(resultsDir, stem + "_effects.csv"));
end
end

function [probabilityFactor, warningFactor] = selectRefinementFactors(effects)
factorNames = unique(effects.Factor, 'stable');
probabilityScores = zeros(size(factorNames));
warningScores = zeros(size(factorNames));
for factorIdx = 1:numel(factorNames)
    rows = effects.Factor == factorNames(factorIdx);
    probabilityScores(factorIdx) = max(abs(effects.ProbabilityEffect_pp(rows)));
    warningScores(factorIdx) = max(abs(effects.WarningEffect_s(rows)));
end
[~, probabilityOrder] = sort(probabilityScores, 'descend');
[~, warningOrder] = sort(warningScores, 'descend');
probabilityFactor = factorNames(probabilityOrder(1));
warningFactor = factorNames(warningOrder(1));
if warningFactor == probabilityFactor && numel(warningOrder) > 1
    warningFactor = factorNames(warningOrder(2));
end
end

function validation = validateAgainstExisting(summary, cfg)
referenceFile = cfg.validation.referenceSummaryFile;
if ~isfile(referenceFile)
    validation = table("Reference file unavailable", false, ...
        'VariableNames', {'Message', 'Pass'});
    return;
end
reference = readtable(referenceFile, TextType='string');
baseline = summary(summary.Factor == "baseline" & ...
    any(summary.Beam == ["line", "ring"], 2), :);

nRows = height(baseline);
Beam = strings(nRows, 1);
Trajectory = strings(nRows, 1);
TargetD_mrad = nan(nRows, 1);
NewProbability = nan(nRows, 1);
ReferenceProbability = nan(nRows, 1);
ProbabilityDifference = nan(nRows, 1);
NewWarning_s = nan(nRows, 1);
ReferenceWarning_s = nan(nRows, 1);
WarningRelativeDifference = nan(nRows, 1);
Pass = false(nRows, 1);

for rowIdx = 1:nRows
    beamName = baseline.Beam(rowIdx);
    pathName = baseline.Trajectory(rowIdx);
    beamCfg = cfg.beams(string({cfg.beams.name}) == beamName);
    targetD = beamCfg.wD_rad * 1e3;
    refRows = reference(lower(reference.Beam) == beamName & ...
        upper(reference.Path) == pathName, :);
    [dSorted, order] = sort(refRows.D_mrad);
    pRef = interp1(dSorted, refRows.Probability(order), targetD, 'linear');
    tRef = interp1(dSorted, refRows.MeanFirstWarningTime_s(order), targetD, 'linear');

    Beam(rowIdx) = beamName;
    Trajectory(rowIdx) = pathName;
    TargetD_mrad(rowIdx) = targetD;
    NewProbability(rowIdx) = baseline.PDetect(rowIdx);
    ReferenceProbability(rowIdx) = pRef;
    ProbabilityDifference(rowIdx) = NewProbability(rowIdx) - pRef;
    NewWarning_s(rowIdx) = baseline.MeanWarning_s(rowIdx);
    ReferenceWarning_s(rowIdx) = tRef;
    WarningRelativeDifference(rowIdx) = abs(NewWarning_s(rowIdx) - tRef) / tRef;
    Pass(rowIdx) = abs(ProbabilityDifference(rowIdx)) <= ...
        cfg.validation.probabilityAbsTolerance && ...
        WarningRelativeDifference(rowIdx) <= ...
        cfg.validation.warningTimeRelativeTolerance;
end
validation = table(Beam, Trajectory, TargetD_mrad, NewProbability, ...
    ReferenceProbability, ProbabilityDifference, NewWarning_s, ...
    ReferenceWarning_s, WarningRelativeDifference, Pass);
if any(~Pass)
    warning('run_sensitivity_analysis:BaselineValidation', ...
        'One or more nominal cases did not meet the reference tolerance.');
end
end

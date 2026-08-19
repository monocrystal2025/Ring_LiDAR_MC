function result = prf_omiga_data(varargin)
%PRF_OMIGA_DATA Generate the scan-sampling robustness data independently.
%
%   This file is the complete, standalone data-processing entry point for
%   the scan sampling-density study
%
%       eta_s = (f/f0) * (wd/wd0) * (Omega0/Omega).
%
%   It does not call any of the other three normalized-study files. Fresh
%   paired ROE scenarios are generated around the Table-1 reference and
%
%       kappa = min(wD_annular) / min(wD_line).
%
%   Public run examples:
%       prf_omiga_data('Mode','preview')
%       prf_omiga_data('Mode','paper')
%       prf_omiga_data('Mode','smoke','SaveData',false)

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end
moduleDir = fullfile(scriptDir, 'ea_sensitivity_analysis');
defaultOutput = fullfile(scriptDir, 'LUBANG_data', ...
    'prf_omiga_data.mat');

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'Mode', 'preview', @isTextScalar);
addParameter(parser, 'N', [], @validOptionalPositiveInteger);
addParameter(parser, 'WDGrid_mrad', [], @validOptionalPositiveVector);
addParameter(parser, 'EtaGrid', [], @validOptionalPositiveVector);
addParameter(parser, 'TimeLimit_s', 15, @validNonnegativeScalar);
addParameter(parser, 'TargetProbability', 0.80, @validProbability);
addParameter(parser, 'UseParallel', true, ...
    @(x) islogical(x) && isscalar(x));
addParameter(parser, 'Seed', 20260817, @validNonnegativeInteger);
addParameter(parser, 'OutputFile', defaultOutput, @isTextScalar);
addParameter(parser, 'SaveData', true, ...
    @(x) islogical(x) && isscalar(x));

parse(parser, varargin{:});

options = parser.Results;
profile = runProfile(lower(string(options.Mode)));
cfg = mergeOptions(profile, options);
cfg.scriptDir = scriptDir;
cfg.moduleDir = moduleDir;
cfg.smallWidth_rad = 1.0e-3;
cfg.smallWidth_mrad = 1.0;
cfg.referencePRF_Hz = 5000;
cfg.referenceScanRate_rad_s = 2 * pi;
cfg.wilsonZ = 1.95996398454005;
cfg.outputFile = char(options.OutputFile);
cfg.saveData = options.SaveData;
cfg.studyName = "scan_sampling";
cfg.panelTitle = "Scan-sampling density";
cfg.xLabel = "\eta_s = (f/f_0)(w_d/w_{d0})(\Omega_0/\Omega)";
cfg.xDescription = "Normalized scan-sampling density";

requiredFiles = { ...
    fullfile(moduleDir, 'ea_sens_default_config.m'), ...
    fullfile(moduleDir, 'ea_sens_build_scenario_bank.m'), ...
    fullfile(moduleDir, 'ea_sens_simulate_case.m')};
for fileIndex = 1:numel(requiredFiles)
    assert(isfile(requiredFiles{fileIndex}), ...
        'prf_omiga_data:MissingSimulator', ...
        'Required simulator file is missing: %s', requiredFiles{fileIndex});
end

pathEntries = string(strsplit(path, pathsep));
addedModulePath = ~any(strcmpi(pathEntries, string(moduleDir)));
if addedModulePath
    addpath(moduleDir);
end
pathCleanup = onCleanup(@() restoreModulePath(moduleDir, addedModulePath));

baseSimulationConfig = buildSimulationConfig(cfg);
caseDefinitions = buildSamplingCases(cfg);
scenarioBank = ea_sens_build_scenario_bank(cfg.N, cfg.Seed);
caseResults = repmat(emptyCaseResult(), numel(caseDefinitions), 1);

fprintf(['Generating %s study: N=%d, wd0=1.0 mrad, ', ...
    '%d wD values, %d normalized cases.\n'], ...
    cfg.studyName, cfg.N, numel(cfg.WDGrid_mrad), ...
    numel(caseDefinitions));
timerValue = tic;

for caseIndex = 1:numel(caseDefinitions)
    definition = caseDefinitions(caseIndex);
    caseConfig = applySystemOverrides( ...
        baseSimulationConfig, definition.overrides);
    fprintf('[%d/%d] %s, x=%.4g\n', caseIndex, ...
        numel(caseDefinitions), definition.caseLabel, definition.xValue);

    lineCurve = simulateBeamCurve( ...
        caseConfig, scenarioBank, 'line', cfg);
    ringCurve = simulateBeamCurve( ...
        caseConfig, scenarioBank, 'ring', cfg);
    equalPerformance = calculateCriterion( ...
        lineCurve, ringCurve, cfg.TargetProbability);

    item = emptyCaseResult();
    item.caseId = definition.caseId;
    item.caseLabel = definition.caseLabel;
    item.series = definition.series;
    item.xValue = definition.xValue;
    item.overrides = definition.overrides;
    item.prf_Hz = caseConfig.system.prf_Hz;
    item.scanRate_rad_s = caseConfig.system.scanRate_rad_s;
    item.wd_mrad = 1e3 * caseConfig.system.beamWidth_rad;
    item.line = lineCurve;
    item.annular = ringCurve;
    item.equalPerformance = equalPerformance;
    caseResults(caseIndex) = item;
end

plotData = buildPlotData(caseResults);
summary = buildSummaryTable(caseResults, plotData);
disp(summary);

result = struct();
result.metadata = struct( ...
    'description', "Fresh normalized-parameter paired ROE Monte Carlo", ...
    'studyName', cfg.studyName, ...
    'panelTitle', cfg.panelTitle, ...
    'xLabel', cfg.xLabel, ...
    'xDescription', cfg.xDescription, ...
    'generatedAt', string(datetime('now', 'Format', ...
    'yyyy-MM-dd HH:mm:ss Z')), ...
    'kappaDefinition', "min(wD_annular)/min(wD_line)", ...
    'historicalDataLoaded', false);
result.config = cfg;
result.caseDefinitions = caseDefinitions;
result.cases = caseResults;
result.plotData = plotData;
result.summary = summary;
result.elapsed_s = toc(timerValue);

if cfg.saveData
    outputFolder = fileparts(cfg.outputFile);
    if ~isfolder(outputFolder)
        mkdir(outputFolder);
    end
    save(cfg.outputFile, 'result', '-v7.3');
    fprintf('Saved MAT data: %s\n', cfg.outputFile);
else
    fprintf('SaveData=false; no MAT file was written.\n');
end
fprintf('Study finished in %.1f s.\n', result.elapsed_s);
end

function tf = isTextScalar(value)
tf = ischar(value) || (isstring(value) && isscalar(value));
end

function tf = validOptionalPositiveInteger(value)
tf = isempty(value) || (isnumeric(value) && isscalar(value) && ...
    isfinite(value) && value >= 1 && fix(value) == value);
end

function tf = validOptionalPositiveVector(value)
tf = isempty(value) || (isnumeric(value) && isvector(value) && ...
    all(isfinite(value)) && all(value > 0));
end

function tf = validNonnegativeScalar(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0;
end

function tf = validProbability(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) && ...
    value > 0 && value <= 1;
end

function tf = validNonnegativeInteger(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) && ...
    value >= 0 && fix(value) == value;
end

function profile = runProfile(mode)
switch mode
    case "smoke"
        profile.N = 4;
        profile.WDGrid_mrad = [80, 140, 200];
        profile.EtaGrid = [0.5, 1, 2];
    case "preview"
        profile.N = 200;
        profile.WDGrid_mrad = [50, 75, 100, 125, 150, 190, 230];
        profile.EtaGrid = [0.25, 0.5, 1, 2, 4];
    case "paper"
        profile.N = 3000;
        profile.WDGrid_mrad = 40:10:260;
        profile.EtaGrid = [0.20, 0.35, 0.50, 0.75, 1, 1.5, 2, 3, 5];
    otherwise
        error('prf_omiga_data:UnknownMode', ...
            'Mode must be smoke, preview or paper.');
end
profile.mode = mode;
end

function cfg = mergeOptions(profile, options)
cfg = profile;
fields = {'N', 'WDGrid_mrad', 'EtaGrid'};
for fieldIndex = 1:numel(fields)
    fieldName = fields{fieldIndex};
    if ~isempty(options.(fieldName))
        cfg.(fieldName) = options.(fieldName);
    end
end
cfg.N = double(cfg.N);
cfg.WDGrid_mrad = unique(double(cfg.WDGrid_mrad(:).'), 'sorted');
cfg.EtaGrid = unique(double(cfg.EtaGrid(:).'), 'sorted');
cfg.TimeLimit_s = double(options.TimeLimit_s);
cfg.TargetProbability = double(options.TargetProbability);
cfg.UseParallel = options.UseParallel;
cfg.Seed = double(options.Seed);
end

function definitions = buildSamplingCases(cfg)
etaGrid = cfg.EtaGrid;
caseCount = 3 * numel(etaGrid);
definitions = repmat(emptyDefinition(), caseCount, 1);
f0 = cfg.referencePRF_Hz;
wd0 = cfg.smallWidth_rad;
omega0 = cfg.referenceScanRate_rad_s;
for etaIndex = 1:numel(etaGrid)
    eta = etaGrid(etaIndex);
    baseIndex = 3 * (etaIndex - 1);
    f = f0 * eta;
    definitions(baseIndex + 1) = makeDefinition( ...
        sprintf('sampling_f_%d', etaIndex), ...
        sprintf('f=%.3g kHz', f / 1e3), 'Vary f', eta, ...
        struct('prf_Hz', f));
    wd = wd0 * eta;
    definitions(baseIndex + 2) = makeDefinition( ...
        sprintf('sampling_wd_%d', etaIndex), ...
        sprintf('wd=%.3g mrad', 1e3 * wd), 'Vary w_d', eta, ...
        struct('beamWidth_rad', wd));
    omega = omega0 / eta;
    definitions(baseIndex + 3) = makeDefinition( ...
        sprintf('sampling_omega_%d', etaIndex), ...
        sprintf('Omega=%.3g rad/s', omega), 'Vary Omega', eta, ...
        struct('scanRate_rad_s', omega));
end
end

function definition = emptyDefinition()
definition = struct('caseId', "", 'caseLabel', "", 'series', "", ...
    'xValue', nan, 'overrides', struct());
end

function definition = makeDefinition(id, label, series, xValue, overrides)
definition = struct('caseId', string(id), 'caseLabel', string(label), ...
    'series', string(series), 'xValue', double(xValue), ...
    'overrides', overrides);
end

function restoreModulePath(moduleDir, addedModulePath)
if addedModulePath
    rmpath(moduleDir);
end
end

function simulationConfig = buildSimulationConfig(cfg)
simulationConfig = ea_sens_default_config();
simulationConfig.seed = cfg.Seed;
simulationConfig.execution.useParallel = cfg.UseParallel;
simulationConfig.execution.showProgress = false;
simulationConfig.execution.blockSize = 25000;

simulationConfig.system.maxRange_m = 2000;
simulationConfig.system.prf_Hz = 5000;
simulationConfig.system.scanRate_rad_s = 2 * pi;
simulationConfig.system.targetSpeed_m_s = 30;
simulationConfig.system.beamWidth_rad = cfg.smallWidth_rad;
simulationConfig.system.pulseEnergy_J = 150e-6;
simulationConfig.system.wavelength_m = 1550e-9;
simulationConfig.system.receiverDiameter_m = 0.0508;
simulationConfig.system.opticalTransmission = 0.80;
simulationConfig.system.quantumEfficiency = 0.80;
simulationConfig.system.rangeGate_s = 50e-9;
simulationConfig.system.filterBandwidth_nm = 10;
simulationConfig.system.darkCountRate_Hz = 400;
simulationConfig.system.extinction_m_inv = 1.5e-5;
simulationConfig.system.backscatter_m_inv_sr_inv = 0.3e-6;
simulationConfig.system.skyRadiance_W_m2_sr_nm = 5e-9;
simulationConfig.system.targetWidth_m = 0.30;
simulationConfig.system.targetReflectivity = 0.50;
simulationConfig.system.snrThreshold = 2;
simulationConfig.system.gaussianCutoffWidths = 3;
simulationConfig.system.targetSampleGrid = 7;
simulationConfig.system.noiseLutStep_m = 0.5;
end

function config = applySystemOverrides(config, overrides)
names = fieldnames(overrides);
for index = 1:numel(names)
    name = names{index};
    assert(isfield(config.system, name), ...
        'prf_omiga_data:UnknownSystemField', ...
        'Unknown system override field: %s', name);
    config.system.(name) = overrides.(name);
end
end

function item = emptyCaseResult()
item = struct('caseId', "", 'caseLabel', "", 'series', "", ...
    'xValue', nan, 'overrides', struct(), 'prf_Hz', nan, ...
    'scanRate_rad_s', nan, 'wd_mrad', nan, ...
    'line', struct(), 'annular', struct(), ...
    'equalPerformance', struct());
end

function curve = simulateBeamCurve(simulationConfig, bank, beamName, cfg)
wDCount = numel(cfg.WDGrid_mrad);
timelyCounts = zeros(1, wDCount);
meanFirstWarning_s = nan(1, wDCount);
meanFirstEncounter_s = nan(1, wDCount);
meanEffectivePulses = nan(1, wDCount);

for wDIndex = 1:wDCount
    beam = struct('name', string(beamName), ...
        'label', string(beamName), ...
        'wD_rad', cfg.WDGrid_mrad(wDIndex) * 1e-3, ...
        'wd_rad', simulationConfig.system.beamWidth_rad);
    raw = ea_sens_simulate_case(simulationConfig, beam, bank);
    detected = logical(raw.detect(:));
    firstTime = raw.firstTime_s(:);
    timely = detected & isfinite(firstTime) & firstTime >= 0 & ...
        firstTime <= cfg.TimeLimit_s;
    timelyCounts(wDIndex) = nnz(timely);

    validWarning = detected & isfinite(firstTime) & firstTime >= 0;
    if any(validWarning)
        meanFirstWarning_s(wDIndex) = mean(firstTime(validWarning));
    end
    encounter = raw.firstEncounter_s(:);
    if any(isfinite(encounter))
        meanFirstEncounter_s(wDIndex) = mean(encounter(isfinite(encounter)));
    end
    pulses = raw.effectivePulses(:);
    if any(isfinite(pulses))
        meanEffectivePulses(wDIndex) = mean(pulses(isfinite(pulses)));
    end
end

probability = timelyCounts ./ bank.N;
[probabilityLow, probabilityHigh] = wilsonInterval( ...
    timelyCounts, bank.N, cfg.wilsonZ);
curve = struct('beam', string(beamName), ...
    'wD_mrad', cfg.WDGrid_mrad, ...
    'wd_mrad', 1e3 * simulationConfig.system.beamWidth_rad, ...
    'N', bank.N, 'timelyCounts', timelyCounts, ...
    'probability', probability, 'probabilityLow95', probabilityLow, ...
    'probabilityHigh95', probabilityHigh, ...
    'meanFirstWarning_s', meanFirstWarning_s, ...
    'meanFirstEncounter_s', meanFirstEncounter_s, ...
    'meanEffectivePulses', meanEffectivePulses);
end

function [low, high] = wilsonInterval(successes, trials, zValue)
phat = successes ./ trials;
denominator = 1 + zValue^2 / trials;
center = (phat + zValue^2 / (2 * trials)) ./ denominator;
halfWidth = zValue .* sqrt(phat .* (1 - phat) ./ trials + ...
    zValue^2 / (4 * trials^2)) ./ denominator;
low = max(0, center - halfWidth);
high = min(1, center + halfWidth);
end

function equal = calculateCriterion(lineCurve, ringCurve, target)
lineMinimum = firstCrossing( ...
    lineCurve.wD_mrad, lineCurve.probability, target);
ringMinimum = firstCrossing( ...
    ringCurve.wD_mrad, ringCurve.probability, target);
lineLower = firstCrossing( ...
    lineCurve.wD_mrad, lineCurve.probabilityHigh95, target);
lineUpper = firstCrossing( ...
    lineCurve.wD_mrad, lineCurve.probabilityLow95, target);
ringLower = firstCrossing( ...
    ringCurve.wD_mrad, ringCurve.probabilityHigh95, target);
ringUpper = firstCrossing( ...
    ringCurve.wD_mrad, ringCurve.probabilityLow95, target);

equal = struct('lineMinWD_mrad', lineMinimum, ...
    'ringMinWD_mrad', ringMinimum, ...
    'kappa', finiteRatio(ringMinimum, lineMinimum), ...
    'kappaLow95', finiteRatio(ringLower, lineUpper), ...
    'kappaHigh95', finiteRatio(ringUpper, lineLower));
equal.robustAdvantage = isfinite(equal.kappaHigh95) && ...
    equal.kappaHigh95 < 1;
equal.advantageDisappears = isfinite(equal.kappa) && equal.kappa >= 1;
end

function minimumWD = firstCrossing(wD, probability, target)
wD = wD(:);
probability = probability(:);
valid = isfinite(wD) & isfinite(probability);
wD = wD(valid);
probability = probability(valid);
if isempty(wD)
    minimumWD = nan;
    return;
end
[wD, order] = sort(wD);
probability = cummax(probability(order));
index = find(probability >= target, 1, 'first');
if isempty(index)
    minimumWD = nan;
elseif index == 1
    minimumWD = wD(1);
else
    p1 = probability(index - 1);
    p2 = probability(index);
    w1 = wD(index - 1);
    w2 = wD(index);
    if p2 <= p1 + eps(max(1, abs(p2)))
        minimumWD = w2;
    else
        minimumWD = w1 + (target - p1) * (w2 - w1) / (p2 - p1);
    end
end
end

function ratio = finiteRatio(numerator, denominator)
if isfinite(numerator) && isfinite(denominator) && denominator > 0
    ratio = numerator / denominator;
else
    ratio = nan;
end
end

function plotData = buildPlotData(caseResults)
count = numel(caseResults);
plotData = struct();
plotData.caseId = strings(count, 1);
plotData.caseLabel = strings(count, 1);
plotData.series = strings(count, 1);
plotData.x = nan(count, 1);
plotData.kappa = nan(count, 1);
plotData.kappaLow95 = nan(count, 1);
plotData.kappaHigh95 = nan(count, 1);
plotData.lineMinWD_mrad = nan(count, 1);
plotData.ringMinWD_mrad = nan(count, 1);
plotData.robustAdvantage = false(count, 1);
plotData.advantageDisappears = false(count, 1);
for index = 1:count
    equal = caseResults(index).equalPerformance;
    plotData.caseId(index) = caseResults(index).caseId;
    plotData.caseLabel(index) = caseResults(index).caseLabel;
    plotData.series(index) = caseResults(index).series;
    plotData.x(index) = caseResults(index).xValue;
    plotData.kappa(index) = equal.kappa;
    plotData.kappaLow95(index) = equal.kappaLow95;
    plotData.kappaHigh95(index) = equal.kappaHigh95;
    plotData.lineMinWD_mrad(index) = equal.lineMinWD_mrad;
    plotData.ringMinWD_mrad(index) = equal.ringMinWD_mrad;
    plotData.robustAdvantage(index) = equal.robustAdvantage;
    plotData.advantageDisappears(index) = equal.advantageDisappears;
end
end

function summary = buildSummaryTable(caseResults, plotData)
CaseId = plotData.caseId;
Series = plotData.series;
NormalizedParameter = plotData.x;
PRF_kHz = reshape([caseResults.prf_Hz], [], 1) ./ 1e3;
ScanRate_rad_s = reshape([caseResults.scanRate_rad_s], [], 1);
wd_mrad = reshape([caseResults.wd_mrad], [], 1);
LineMinWD_mrad = plotData.lineMinWD_mrad;
AnnularMinWD_mrad = plotData.ringMinWD_mrad;
Kappa = plotData.kappa;
KappaLow95 = plotData.kappaLow95;
KappaHigh95 = plotData.kappaHigh95;
RobustAdvantage = plotData.robustAdvantage;
AdvantageDisappears = plotData.advantageDisappears;
summary = table(CaseId, Series, NormalizedParameter, PRF_kHz, ...
    ScanRate_rad_s, wd_mrad, LineMinWD_mrad, AnnularMinWD_mrad, Kappa, ...
    KappaLow95, KappaHigh95, RobustAdvantage, AdvantageDisappears);
summary = sortrows(summary, {'Series', 'NormalizedParameter'});
end

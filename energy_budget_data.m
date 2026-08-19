function result = energy_budget_data(varargin)
%ENERGY_BUDGET_DATA Independent laser-return energy-budget study.
%
%   eta_E = (Ep/Ep0)(Dr/Dr0)^2(Wt/Wt0)^2(SNR0/SNRth).
%
%   This file performs its own case construction, Monte Carlo simulation,
%   confidence calculation and MAT-file output. It does not call any of the
%   other three normalized-study data files. wd is fixed at 1.0 mrad.

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
moduleDir = fullfile(scriptDir, 'ea_sensitivity_analysis');
defaultOutput = fullfile(scriptDir, 'LUBANG_data', ...
    'energy_budget_data.mat');

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Mode', 'preview', @isTextScalar);
addParameter(p, 'N', [], @optionalInteger);
addParameter(p, 'WDGrid_mrad', [], @optionalVector);
addParameter(p, 'TimeLimit_s', 15, @nonnegativeScalar);
addParameter(p, 'TargetProbability', 0.80, @probabilityScalar);
addParameter(p, 'UseParallel', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Seed', 20260817, @nonnegativeInteger);
addParameter(p, 'OutputFile', defaultOutput, @isTextScalar);
addParameter(p, 'SaveData', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
cfg = makeRunConfig(p.Results, scriptDir, moduleDir);
verifySimulator(moduleDir);
addedPath = addSimulatorPath(moduleDir);
pathCleanup = onCleanup(@() restorePath(moduleDir, addedPath));

definitions = buildEnergyCases();
baseConfig = makeSimulationConfig(cfg);
scenarioBank = ea_sens_build_scenario_bank(cfg.N, cfg.Seed);
cases = repmat(emptyCase(), numel(definitions), 1);
fprintf(['Generating independent energy-budget study: N=%d, ', ...
    'wd=1.0 mrad, %d wD values, %d cases.\n'], ...
    cfg.N, numel(cfg.WDGrid_mrad), numel(definitions));
timerValue = tic;
for k = 1:numel(definitions)
    fprintf('[%d/%d] %s, eta_E=%.4g\n', k, numel(definitions), ...
        definitions(k).caseLabel, definitions(k).xValue);
    simConfig = applyOverrides(baseConfig, definitions(k).overrides);
    line = simulateCurve(simConfig, scenarioBank, 'line', cfg);
    annular = simulateCurve(simConfig, scenarioBank, 'ring', cfg);
    cases(k) = packCase(definitions(k), line, annular, ...
        equalPerformance(line, annular, cfg.TargetProbability));
end

plotData = makePlotData(cases);
summary = makeSummary(plotData);
disp(summary);
result = packResult(cfg, definitions, cases, plotData, summary, ...
    toc(timerValue));
if cfg.SaveData
    outputFolder = fileparts(cfg.OutputFile);
    if ~isfolder(outputFolder), mkdir(outputFolder); end
    save(cfg.OutputFile, 'result', '-v7.3');
    fprintf('Saved MAT data: %s\n', cfg.OutputFile);
end
fprintf('Energy-budget study finished in %.1f s.\n', result.elapsed_s);
end

function definitions = buildEnergyCases()
% Log-like spacing covers two decades without an excessive case count.
ratios = [0.10, 0.20, 0.50, 1.00, 2.00, 5.00, 10.0];
definitions = repmat(emptyDefinition(), 4 * numel(ratios), 1);
k = 0;
for ratio = ratios
    k = k + 1;
    definitions(k) = definition(k, sprintf('Ep ratio %.3g', ratio), ...
        'Pulse energy', ratio, struct('pulseEnergy_J', 150e-6 * ratio));
    k = k + 1;
    definitions(k) = definition(k, ...
        sprintf('Aperture-area ratio %.3g', ratio), ...
        'Receiver aperture', ratio, ...
        struct('receiverDiameter_m', 0.0508 * sqrt(ratio)));
    k = k + 1;
    definitions(k) = definition(k, ...
        sprintf('Sensitivity ratio %.3g', ratio), ...
        'Detection threshold', ratio, ...
        struct('snrThreshold', 2 / ratio));
    k = k + 1;
    definitions(k) = definition(k, ...
        sprintf('Target-area ratio %.3g', ratio), 'Target size', ...
        ratio, struct('targetWidth_m', 0.30 * sqrt(ratio)));
end
end

function cfg = makeRunConfig(options, scriptDir, moduleDir)
switch lower(string(options.Mode))
    case "smoke"
        defaultN = 4; defaultWD = [80, 140, 200];
    case "preview"
        defaultN = 200; defaultWD = [50, 75, 100, 125, 150, 190, 230];
    case "paper"
        defaultN = 3000; defaultWD = 40:10:260;
    otherwise
        error('energy_budget_data:UnknownMode', ...
            'Mode must be smoke, preview or paper.');
end
cfg = struct();
cfg.mode = lower(string(options.Mode));
if isempty(options.N), cfg.N = defaultN; else, cfg.N = options.N; end
if isempty(options.WDGrid_mrad)
    cfg.WDGrid_mrad = defaultWD;
else
    cfg.WDGrid_mrad = unique(double(options.WDGrid_mrad(:).'), 'sorted');
end
cfg.TimeLimit_s = double(options.TimeLimit_s);
cfg.TargetProbability = double(options.TargetProbability);
cfg.UseParallel = options.UseParallel;
cfg.Seed = double(options.Seed);
cfg.OutputFile = char(options.OutputFile);
cfg.SaveData = options.SaveData;
cfg.smallWidth_mrad = 1.0;
cfg.scriptDir = scriptDir;
cfg.moduleDir = moduleDir;
cfg.wilsonZ = 1.95996398454005;
end

function config = makeSimulationConfig(cfg)
config = ea_sens_default_config();
config.seed = cfg.Seed;
config.execution.useParallel = cfg.UseParallel;
config.execution.showProgress = false;
config.execution.blockSize = 25000;
config.system.maxRange_m = 2000;
config.system.prf_Hz = 5000;
config.system.scanRate_rad_s = 2*pi;
config.system.targetSpeed_m_s = 30;
config.system.beamWidth_rad = 1e-3;
config.system.pulseEnergy_J = 150e-6;
config.system.wavelength_m = 1550e-9;
config.system.receiverDiameter_m = 0.0508;
config.system.opticalTransmission = 0.80;
config.system.quantumEfficiency = 0.80;
config.system.rangeGate_s = 50e-9;
config.system.filterBandwidth_nm = 10;
config.system.darkCountRate_Hz = 400;
config.system.extinction_m_inv = 1.5e-5;
config.system.backscatter_m_inv_sr_inv = 0.3e-6;
config.system.skyRadiance_W_m2_sr_nm = 5e-9;
config.system.targetWidth_m = 0.30;
config.system.targetReflectivity = 0.50;
config.system.snrThreshold = 2;
config.system.gaussianCutoffWidths = 3;
config.system.targetSampleGrid = 7;
config.system.noiseLutStep_m = 0.5;
end

function config = applyOverrides(config, overrides)
names = fieldnames(overrides);
for k = 1:numel(names), config.system.(names{k}) = overrides.(names{k}); end
config.system.beamWidth_rad = 1e-3;
end

function curve = simulateCurve(config, bank, beamName, cfg)
nWD = numel(cfg.WDGrid_mrad);
counts = zeros(1, nWD);
for j = 1:nWD
    beam = struct('name', string(beamName), 'label', string(beamName), ...
        'wD_rad', cfg.WDGrid_mrad(j)*1e-3, 'wd_rad', 1e-3);
    raw = ea_sens_simulate_case(config, beam, bank);
    firstTime = raw.firstTime_s(:);
    timely = logical(raw.detect(:)) & isfinite(firstTime) & ...
        firstTime >= 0 & firstTime <= cfg.TimeLimit_s;
    counts(j) = nnz(timely);
end
prob = counts ./ bank.N;
[low, high] = wilson(counts, bank.N, cfg.wilsonZ);
curve = struct('beam', string(beamName), 'wD_mrad', cfg.WDGrid_mrad, ...
    'wd_mrad', 1.0, 'N', bank.N, 'timelyCounts', counts, ...
    'probability', prob, 'probabilityLow95', low, ...
    'probabilityHigh95', high);
end

function equal = equalPerformance(line, ring, target)
lineMin = crossing(line.wD_mrad, line.probability, target);
ringMin = crossing(ring.wD_mrad, ring.probability, target);
lineLow = crossing(line.wD_mrad, line.probabilityHigh95, target);
lineHigh = crossing(line.wD_mrad, line.probabilityLow95, target);
ringLow = crossing(ring.wD_mrad, ring.probabilityHigh95, target);
ringHigh = crossing(ring.wD_mrad, ring.probabilityLow95, target);
equal = struct('lineMinWD_mrad', lineMin, 'ringMinWD_mrad', ringMin, ...
    'kappa', ratio(ringMin, lineMin), ...
    'kappaLow95', ratio(ringLow, lineHigh), ...
    'kappaHigh95', ratio(ringHigh, lineLow));
equal.robustAdvantage = isfinite(equal.kappaHigh95) && equal.kappaHigh95 < 1;
equal.advantageDisappears = isfinite(equal.kappa) && equal.kappa >= 1;
end

function value = crossing(wD, probability, target)
wD = wD(:); probability = probability(:);
[wD, order] = sort(wD); probability = probability(order);
valid = isfinite(wD) & isfinite(probability);
wD = wD(valid); probability = cummax(probability(valid));
index = find(probability >= target, 1);
if isempty(index), value = nan; return; end
if index == 1, value = wD(1); return; end
p1 = probability(index-1); p2 = probability(index);
if p2 <= p1 + eps(max(1, abs(p2)))
    value = wD(index);
else
    value = wD(index-1) + (target-p1) * ...
        (wD(index)-wD(index-1)) / (p2-p1);
end
end

function [low, high] = wilson(successes, n, z)
p = successes ./ n; d = 1 + z^2/n;
c = (p + z^2/(2*n)) ./ d;
h = z .* sqrt(p.*(1-p)./n + z^2/(4*n^2)) ./ d;
low = max(0, c-h); high = min(1, c+h);
end

function value = ratio(a, b)
if isfinite(a) && isfinite(b) && b > 0, value = a/b; else, value = nan; end
end

function item = packCase(def, line, annular, equal)
item = emptyCase();
item.caseId = def.caseId; item.caseLabel = def.caseLabel;
item.series = def.series; item.xValue = def.xValue;
item.overrides = def.overrides; item.line = line;
item.annular = annular; item.equalPerformance = equal;
end

function data = makePlotData(cases)
n = numel(cases);
data = struct('caseId', strings(n,1), 'caseLabel', strings(n,1), ...
    'series', strings(n,1), 'x', nan(n,1), 'kappa', nan(n,1), ...
    'kappaLow95', nan(n,1), 'kappaHigh95', nan(n,1), ...
    'lineMinWD_mrad', nan(n,1), 'ringMinWD_mrad', nan(n,1), ...
    'robustAdvantage', false(n,1), 'advantageDisappears', false(n,1));
for k = 1:n
    e = cases(k).equalPerformance;
    data.caseId(k)=cases(k).caseId; data.caseLabel(k)=cases(k).caseLabel;
    data.series(k)=cases(k).series; data.x(k)=cases(k).xValue;
    data.kappa(k)=e.kappa; data.kappaLow95(k)=e.kappaLow95;
    data.kappaHigh95(k)=e.kappaHigh95;
    data.lineMinWD_mrad(k)=e.lineMinWD_mrad;
    data.ringMinWD_mrad(k)=e.ringMinWD_mrad;
    data.robustAdvantage(k)=e.robustAdvantage;
    data.advantageDisappears(k)=e.advantageDisappears;
end
end

function summary = makeSummary(d)
summary = table(d.caseId,d.series,d.x,d.lineMinWD_mrad,d.ringMinWD_mrad, ...
    d.kappa,d.kappaLow95,d.kappaHigh95,d.robustAdvantage, ...
    d.advantageDisappears, 'VariableNames', {'CaseId','Series', ...
    'NormalizedParameter','LineMinWD_mrad','AnnularMinWD_mrad','Kappa', ...
    'KappaLow95','KappaHigh95','RobustAdvantage','AdvantageDisappears'});
summary = sortrows(summary, {'Series','NormalizedParameter'});
end

function result = packResult(cfg, definitions, cases, plotData, summary, elapsed)
metadata = struct('description', "Independent energy-budget paired ROE Monte Carlo", ...
    'studyName', "energy_budget", 'panelTitle', "Laser-return energy budget", ...
    'xLabel', "\eta_E (normalized photon-return budget)", ...
    'xDescription', "Ep * receiver area * target area / detection threshold", ...
    'generatedAt', string(datetime('now','Format','yyyy-MM-dd HH:mm:ss Z')), ...
    'kappaDefinition', "min(wD_annular)/min(wD_line)", ...
    'historicalDataLoaded', false, 'independentEntryPoint', true);
result = struct('metadata',metadata,'config',cfg, ...
    'caseDefinitions',definitions,'cases',cases,'plotData',plotData, ...
    'summary',summary,'elapsed_s',elapsed);
end

function d = definition(index, label, series, x, overrides)
d = struct('caseId',"energy_"+index,'caseLabel',string(label), ...
    'series',string(series),'xValue',double(x),'overrides',overrides);
end
function d = emptyDefinition()
d = struct('caseId',"",'caseLabel',"",'series',"", ...
    'xValue',nan,'overrides',struct());
end
function c = emptyCase()
c = struct('caseId',"",'caseLabel',"",'series',"",'xValue',nan, ...
    'overrides',struct(),'line',struct(),'annular',struct(), ...
    'equalPerformance',struct());
end

function verifySimulator(moduleDir)
names = {'ea_sens_default_config.m','ea_sens_build_scenario_bank.m', ...
    'ea_sens_simulate_case.m'};
for k=1:numel(names)
    assert(isfile(fullfile(moduleDir,names{k})), ...
        'energy_budget_data:MissingSimulator', ...
        'Required simulator file is missing: %s', fullfile(moduleDir,names{k}));
end
end
function added = addSimulatorPath(moduleDir)
entries = string(strsplit(path,pathsep)); added = ~any(strcmpi(entries,moduleDir));
if added, addpath(moduleDir); end
end
function restorePath(moduleDir, added)
if added, rmpath(moduleDir); end
end
function tf=isTextScalar(x), tf=ischar(x)||(isstring(x)&&isscalar(x)); end
function tf=optionalInteger(x), tf=isempty(x)||(isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&fix(x)==x); end
function tf=optionalVector(x), tf=isempty(x)||(isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>0)); end
function tf=nonnegativeScalar(x), tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0; end
function tf=probabilityScalar(x), tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<=1; end
function tf=nonnegativeInteger(x), tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&fix(x)==x; end

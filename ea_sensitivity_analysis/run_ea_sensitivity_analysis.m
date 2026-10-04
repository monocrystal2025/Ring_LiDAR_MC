function outputs = run_ea_sensitivity_analysis(varargin)
%RUN_EA_SENSITIVITY_ANALYSIS Run the standalone ROE sensitivity workflow.
%
% Examples
%   outputs = run_ea_sensitivity_analysis;
%   outputs = run_ea_sensitivity_analysis("Mode", "screen");
%   outputs = run_ea_sensitivity_analysis("Mode", "smoke", ...
%       "UseParallel", false, "Force", true);
%   run_ea_sensitivity_analysis("Mode", "plotonly");
%
% The default full run uses N=3000 for both endpoint screening and refined
% cases. It never starts a ten-thousand-level Monte Carlo run.

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, "Mode", "full", @(x) any(strcmpi(string(x), ...
    ["full", "screen", "smoke", "plotonly"])));
addParameter(parser, "OutputRoot", "", @(x) isstring(x) || ischar(x));
addParameter(parser, "ScreenN", [], @(x) isempty(x) || ...
    (isscalar(x) && x > 0 && fix(x) == x));
addParameter(parser, "RefineN", [], @(x) isempty(x) || ...
    (isscalar(x) && x > 0 && fix(x) == x));
addParameter(parser, "BootstrapN", [], @(x) isempty(x) || ...
    (isscalar(x) && x >= 0 && fix(x) == x));
addParameter(parser, "UseParallel", [], @(x) isempty(x) || ...
    islogical(x) || isnumeric(x));
addParameter(parser, "Force", false, @(x) islogical(x) || isnumeric(x));
parse(parser, varargin{:});
options = parser.Results;
mode = lower(string(options.Mode));

cfg = ea_sens_default_config();
cfg = applyOptions(cfg, options, mode);
ensureOutputFolders(cfg);

if mode == "plotonly"
    figures = ea_sens_plot_results(cfg.outputRoot);
    outputs = struct("config", cfg, "figures", figures);
    return;
end

screenBank = ea_sens_build_scenario_bank(cfg.screenN, cfg.seed);
screenCases = buildScreenCases(cfg);
screenRaw = runCases(cfg, screenCases, screenBank, ...
    cfg.wDGrid_rad, "screen", logical(options.Force));
[screenSummary, screenTimely, screenIso, screenEffects] = ...
    ea_sens_analyze_results(screenRaw, cfg);
writeAnalysisTables(cfg.resultsDir, "ea_sensitivity_screen", ...
    screenSummary, screenTimely, screenIso, screenEffects);
save(fullfile(cfg.resultsDir, "ea_sensitivity_screen_results.mat"), ...
    "cfg", "screenRaw", "screenSummary", "screenTimely", ...
    "screenIso", "screenEffects", "-v7.3");

[firstFactor, secondFactor] = selectRefinementFactors(screenEffects, cfg);
selectedFactors = table(firstFactor, secondFactor, ...
    "VariableNames", {"PrimaryFactor", "SecondaryFactor"});
writetable(selectedFactors, ...
    fullfile(cfg.resultsDir, "selected_factors.csv"));

outputs = struct();
outputs.config = cfg;
outputs.screenSummary = screenSummary;
outputs.screenTimely = screenTimely;
outputs.screenIso = screenIso;
outputs.screenEffects = screenEffects;
outputs.selectedFactors = selectedFactors;

if mode == "screen"
    outputs.figures = ea_sens_plot_results(cfg.outputRoot);
    return;
end

refineBank = ea_sens_build_scenario_bank(cfg.refineN, cfg.seed + 1);
refineCases = buildRefinementCases(cfg, firstFactor, secondFactor);
refineRaw = runCases(cfg, refineCases, refineBank, ...
    cfg.refineWDGrid_rad, "refine", logical(options.Force));
[refineSummary, refineTimely, refineIso] = ...
    ea_sens_analyze_results(refineRaw, cfg);
writeAnalysisTables(cfg.resultsDir, "ea_sensitivity_refined", ...
    refineSummary, refineTimely, refineIso, table());
save(fullfile(cfg.resultsDir, "ea_sensitivity_refined_results.mat"), ...
    "cfg", "refineRaw", "refineSummary", "refineTimely", ...
    "refineIso", "selectedFactors", "-v7.3");

interactionCases = buildInteractionCases(cfg, firstFactor, secondFactor);
interactionRaw = runCases(cfg, interactionCases, refineBank, ...
    cfg.refineWDGrid_rad, "interaction", logical(options.Force));
[interactionSummary, interactionTimely, interactionIso] = ...
    ea_sens_analyze_results(interactionRaw, cfg);
interactionIso = attachInteractionMetadata(interactionIso, interactionCases);
writeAnalysisTables(cfg.resultsDir, "ea_sensitivity_interaction", ...
    interactionSummary, interactionTimely, interactionIso, table());
save(fullfile(cfg.resultsDir, "ea_sensitivity_interaction_results.mat"), ...
    "cfg", "interactionRaw", "interactionSummary", ...
    "interactionTimely", "interactionIso", "selectedFactors", "-v7.3");

outputs.refineSummary = refineSummary;
outputs.refineTimely = refineTimely;
outputs.refineIso = refineIso;
outputs.interactionSummary = interactionSummary;
outputs.interactionTimely = interactionTimely;
outputs.interactionIso = interactionIso;
outputs.figures = ea_sens_plot_results(cfg.outputRoot);
end

function cfg = applyOptions(cfg, options, mode)
if strlength(string(options.OutputRoot)) > 0
    cfg.outputRoot = char(string(options.OutputRoot));
    cfg.resultsDir = fullfile(cfg.outputRoot, "results");
    cfg.figuresDir = fullfile(cfg.outputRoot, "figures");
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
    cfg.screenN = min(cfg.screenN, 6);
    cfg.refineN = min(cfg.refineN, 6);
    cfg.bootstrapN = min(cfg.bootstrapN, 10);
    cfg.execution.useParallel = false;
    cfg.execution.blockSize = 500;
    cfg.execution.showProgress = false;
    cfg.system.maxRange_m = 100;
    cfg.system.prf_Hz = 100;
    cfg.system.scanRate_rad_s = 2 * pi;
    cfg.system.targetSpeed_m_s = 20;
    cfg.system.noiseLutStep_m = 2;
    cfg.wDGrid_rad = [80, 140, 220] * 1e-3;
    cfg.refineWDGrid_rad = cfg.wDGrid_rad;
    cfg.timelyThresholds_s = 5;
    cfg.targetProbabilities = 0.05;
    cfg.primaryTimelyThreshold_s = 5;
    cfg.primaryProbability = 0.05;
    cfg.factors = cfg.factors(ismember( ...
        string({cfg.factors.name}), ["beamWidth", "targetSpeed"]));
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
baseline = makeCase("baseline", "baseline", "Baseline", ...
    "Baseline", 0, "", cfg);
cases = baseline;
for factorIdx = 1:numel(cfg.factors)
    factor = cfg.factors(factorIdx);
    lowCfg = applyFactor(cfg, factor, factor.lowValue);
    highCfg = applyFactor(cfg, factor, factor.highValue);
    cases(end+1) = makeCase(factor.name + "_low", factor.name, ...
        factor.label, "Low", factor.lowValue, factor.unit, lowCfg); %#ok<AGROW>
    cases(end+1) = makeCase(factor.name + "_high", factor.name, ...
        factor.label, "High", factor.highValue, factor.unit, highCfg); %#ok<AGROW>
end
end

function cases = buildRefinementCases(cfg, firstFactorName, secondFactorName)
factorNames = [firstFactorName, secondFactorName];
cases = repmat(emptyCase(), 1, 0);
for selectedIdx = 1:numel(factorNames)
    factor = findFactor(cfg, factorNames(selectedIdx));
    values = linspace(factor.lowValue, factor.highValue, 5);
    for levelIdx = 1:numel(values)
        caseCfg = applyFactor(cfg, factor, values(levelIdx));
        caseId = factor.name + "_refine_" + levelIdx;
        level = "Level" + levelIdx;
        cases(end+1) = makeCase(caseId, factor.name, factor.label, ...
            level, values(levelIdx), factor.unit, caseCfg); %#ok<AGROW>
    end
end
end

function cases = buildInteractionCases(cfg, firstFactorName, secondFactorName)
first = findFactor(cfg, firstFactorName);
second = findFactor(cfg, secondFactorName);
firstValues = [first.lowValue, factorBaseline(cfg, first), first.highValue];
secondValues = [second.lowValue, factorBaseline(cfg, second), second.highValue];
cases = repmat(emptyCase(), 1, 0);

for firstIdx = 1:3
    for secondIdx = 1:3
        caseCfg = applyFactor(cfg, first, firstValues(firstIdx));
        caseCfg = applyFactor(caseCfg, second, secondValues(secondIdx));
        caseId = "interaction_" + first.name + "_" + firstIdx + ...
            "_" + second.name + "_" + secondIdx;
        item = makeCase(caseId, "interaction", ...
            first.label + " x " + second.label, ...
            "I" + firstIdx + "J" + secondIdx, nan, "", caseCfg);
        item.factor1 = first.name;
        item.factor1Label = first.label;
        item.factor1Value = firstValues(firstIdx);
        item.factor1Unit = first.unit;
        item.factor2 = second.name;
        item.factor2Label = second.label;
        item.factor2Value = secondValues(secondIdx);
        item.factor2Unit = second.unit;
        cases(end+1) = item; %#ok<AGROW>
    end
end
end

function factor = findFactor(cfg, name)
factorNames = string({cfg.factors.name});
idx = find(factorNames == name, 1);
if isempty(idx)
    error("run_ea_sensitivity_analysis:UnknownFactor", ...
        "Unknown factor: %s.", name);
end
factor = cfg.factors(idx);
end

function value = factorBaseline(cfg, factor)
if factor.key == "environmentSeverity"
    value = 0;
else
    value = cfg.system.(factor.key);
end
end

function cfgOut = applyFactor(cfgIn, factor, value)
cfgOut = cfgIn;
if factor.key == "environmentSeverity"
    severity = max(-1, min(1, value));
    if severity < 0
        extinctionScale = 1 + (-severity) * ...
            (cfgIn.environment.clearExtinctionScale - 1);
        backscatterScale = 1 + (-severity) * ...
            (cfgIn.environment.clearBackscatterScale - 1);
        radianceScale = 1 + (-severity) * ...
            (cfgIn.environment.clearRadianceScale - 1);
    else
        extinctionScale = 1 + severity * ...
            (cfgIn.environment.adverseExtinctionScale - 1);
        backscatterScale = 1 + severity * ...
            (cfgIn.environment.adverseBackscatterScale - 1);
        radianceScale = 1 + severity * ...
            (cfgIn.environment.adverseRadianceScale - 1);
    end
    cfgOut.system.extinction_m_inv = ...
        cfgIn.environment.nominalExtinction_m_inv * extinctionScale;
    cfgOut.system.backscatter_m_inv_sr_inv = ...
        cfgIn.environment.nominalBackscatter_m_inv_sr_inv * backscatterScale;
    cfgOut.system.skyRadiance_W_m2_sr_nm = ...
        cfgIn.environment.nominalSkyRadiance_W_m2_sr_nm * radianceScale;
else
    cfgOut.system.(factor.key) = value;
end
end

function item = makeCase(caseId, factor, factorLabel, level, ...
        levelValue, unit, config)
item = emptyCase();
item.caseId = string(caseId);
item.factor = string(factor);
item.factorLabel = string(factorLabel);
item.level = string(level);
item.levelValue = levelValue;
item.unit = string(unit);
item.config = config;
end

function item = emptyCase()
item = struct( ...
    "caseId", "", ...
    "factor", "", ...
    "factorLabel", "", ...
    "level", "", ...
    "levelValue", nan, ...
    "unit", "", ...
    "config", struct(), ...
    "factor1", "", ...
    "factor1Label", "", ...
    "factor1Value", nan, ...
    "factor1Unit", "", ...
    "factor2", "", ...
    "factor2Label", "", ...
    "factor2Value", nan, ...
    "factor2Unit", "");
end

function raw = runCases(cfg, cases, bank, wDGrid, stage, force)
caseItemCount = numel(wDGrid) * numel(cfg.beamNames);
raw = repmat(emptyRaw(), 1, numel(cases) * caseItemCount);
writeIdx = 0;

for caseIdx = 1:numel(cases)
    caseInfo = cases(caseIdx);
    checkpoint = fullfile(cfg.resultsDir, ...
        "raw_" + stage + "_" + safeFilePart(caseInfo.caseId) + ".mat");
    loaded = false;
    if isfile(checkpoint) && ~force
        try
            saved = load(checkpoint, "caseRaw", "bankSeed", "bankN", ...
                "savedWDGrid", "savedSystem", "savedVersion");
            requiredSaved = {'caseRaw', 'bankSeed', 'bankN', ...
                'savedWDGrid', 'savedSystem', 'savedVersion'};
            loaded = all(isfield(saved, requiredSaved)) && ...
                saved.bankSeed == bank.seed && saved.bankN == bank.N && ...
                isequal(saved.savedWDGrid, wDGrid) && ...
                isequaln(saved.savedSystem, caseInfo.config.system) && ...
                string(saved.savedVersion) == cfg.version && ...
                numel(saved.caseRaw) == caseItemCount;
            if loaded
                caseRaw = saved.caseRaw;
            end
        catch checkpointError
            warning("run_ea_sensitivity_analysis:InvalidCheckpoint", ...
                "Ignoring unreadable checkpoint %s (%s).", ...
                checkpoint, checkpointError.message);
        end
    end

    if ~loaded
        caseRaw = repmat(emptyRaw(), 1, caseItemCount);
        localIdx = 0;
        for wDIdx = 1:numel(wDGrid)
            beams = makeBeams(caseInfo.config, wDGrid(wDIdx), cfg.beamNames);
            for beamIdx = 1:numel(beams)
                if cfg.execution.showProgress
                    fprintf("[%s %d/%d] %s, %s, wD=%.1f mrad\n", ...
                        upper(stage), caseIdx, numel(cases), ...
                        caseInfo.caseId, beams(beamIdx).name, ...
                        1e3 * wDGrid(wDIdx));
                end
                sim = ea_sens_simulate_case( ...
                    caseInfo.config, beams(beamIdx), bank);
                localIdx = localIdx + 1;
                caseRaw(localIdx) = packageRaw(sim, caseInfo);
            end
        end
        bankSeed = bank.seed;
        bankN = bank.N;
        savedWDGrid = wDGrid;
        savedSystem = caseInfo.config.system;
        savedVersion = cfg.version;
        save(checkpoint, "caseRaw", "bankSeed", "bankN", ...
            "savedWDGrid", "savedSystem", "savedVersion", "-v7.3");
    end

    raw(writeIdx + (1:numel(caseRaw))) = caseRaw;
    writeIdx = writeIdx + numel(caseRaw);
end
raw = raw(1:writeIdx);
end

function beams = makeBeams(cfg, wD, beamNames)
beams = repmat(struct("name", "", "label", "", ...
    "wD_rad", nan, "wd_rad", nan), 1, numel(beamNames));
for beamIdx = 1:numel(beamNames)
    name = lower(string(beamNames(beamIdx)));
    beams(beamIdx).name = name;
    beams(beamIdx).wD_rad = wD;
    if name == "point"
        beams(beamIdx).label = "Spot";
        beams(beamIdx).wd_rad = nan;
    elseif name == "line"
        beams(beamIdx).label = "Line";
        beams(beamIdx).wd_rad = cfg.system.beamWidth_rad;
    elseif name == "ring"
        beams(beamIdx).label = "Annular";
        beams(beamIdx).wd_rad = cfg.system.beamWidth_rad;
    else
        error("run_ea_sensitivity_analysis:UnknownBeam", ...
            "Unknown beam name: %s.", name);
    end
end
end

function item = packageRaw(sim, caseInfo)
item = emptyRaw();
item.caseId = caseInfo.caseId;
item.factor = caseInfo.factor;
item.factorLabel = caseInfo.factorLabel;
item.level = caseInfo.level;
item.levelValue = caseInfo.levelValue;
item.unit = caseInfo.unit;
item.beam = sim.beam;
item.beamLabel = sim.beamLabel;
item.wD_rad = sim.wD_rad;
item.wd_rad = sim.wd_rad;
item.scenarioId = sim.scenarioId;
item.detect = sim.detect;
item.firstTime_s = sim.firstTime_s;
item.firstEncounter_s = sim.firstEncounter_s;
item.effectivePulses = sim.effectivePulses;
item.pathPointCount = sim.pathPointCount;
item.cycleLength = sim.cycleLength;
item.windowPulses = sim.windowPulses;
item.elapsed_s = sim.elapsed_s;
end

function item = emptyRaw()
item = struct( ...
    "caseId", "", ...
    "factor", "", ...
    "factorLabel", "", ...
    "level", "", ...
    "levelValue", nan, ...
    "unit", "", ...
    "beam", "", ...
    "beamLabel", "", ...
    "wD_rad", nan, ...
    "wd_rad", nan, ...
    "scenarioId", zeros(0, 1), ...
    "detect", false(0, 1), ...
    "firstTime_s", zeros(0, 1), ...
    "firstEncounter_s", zeros(0, 1), ...
    "effectivePulses", zeros(0, 1), ...
    "pathPointCount", nan, ...
    "cycleLength", nan, ...
    "windowPulses", nan, ...
    "elapsed_s", nan);
end

function [firstFactor, secondFactor] = selectRefinementFactors(effects, cfg)
valid = effects(isfinite(effects.AbsEffect), :);
if height(valid) >= 2
    firstFactor = valid.Factor(1);
    secondFactor = valid.Factor(2);
elseif height(valid) == 1
    firstFactor = valid.Factor(1);
    fallback = string({cfg.factors.name});
    secondFactor = fallback(find(fallback ~= firstFactor, 1));
else
    fallback = string({cfg.factors.name});
    firstFactor = fallback(1);
    secondFactor = fallback(min(2, numel(fallback)));
end
end

function iso = attachInteractionMetadata(iso, cases)
caseIds = string({cases.caseId});
Factor1 = strings(height(iso), 1);
Factor1Label = strings(height(iso), 1);
Factor1Value = nan(height(iso), 1);
Factor1Unit = strings(height(iso), 1);
Factor2 = strings(height(iso), 1);
Factor2Label = strings(height(iso), 1);
Factor2Value = nan(height(iso), 1);
Factor2Unit = strings(height(iso), 1);
for rowIdx = 1:height(iso)
    caseIdx = find(caseIds == iso.CaseId(rowIdx), 1);
    Factor1(rowIdx) = cases(caseIdx).factor1;
    Factor1Label(rowIdx) = cases(caseIdx).factor1Label;
    Factor1Value(rowIdx) = cases(caseIdx).factor1Value;
    Factor1Unit(rowIdx) = cases(caseIdx).factor1Unit;
    Factor2(rowIdx) = cases(caseIdx).factor2;
    Factor2Label(rowIdx) = cases(caseIdx).factor2Label;
    Factor2Value(rowIdx) = cases(caseIdx).factor2Value;
    Factor2Unit(rowIdx) = cases(caseIdx).factor2Unit;
end
iso = addvars(iso, Factor1, Factor1Label, Factor1Value, Factor1Unit, ...
    Factor2, Factor2Label, Factor2Value, Factor2Unit, ...
    'After', 'Unit');
end

function writeAnalysisTables(resultsDir, stem, summary, timely, iso, effects)
writetable(summary, fullfile(resultsDir, stem + "_summary.csv"));
writetable(timely, fullfile(resultsDir, stem + "_timely.csv"));
writetable(iso, fullfile(resultsDir, stem + "_iso_performance.csv"));
if ~isempty(effects)
    writetable(effects, fullfile(resultsDir, stem + "_effects.csv"));
end
end

function value = safeFilePart(value)
value = regexprep(char(string(value)), "[^A-Za-z0-9_-]", "_");
end

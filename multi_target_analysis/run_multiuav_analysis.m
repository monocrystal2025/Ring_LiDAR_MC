function output = run_multiuav_analysis(varargin)
%RUN_MULTIUAV_ANALYSIS Run multi-UAV studies, including formation wD search.
%
% output = run_multiuav_analysis('Mode','smoke','N',24)

parser = inputParser;
addParameter(parser, 'Mode', "smoke", @(x) any(lower(string(x)) == ...
    ["smoke", "preview", "screen", "formal", "formationopt", "plotonly"]));
addParameter(parser, 'N', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(parser, 'OutputRoot', "", @(x) ischar(x) || isstring(x));
addParameter(parser, 'UseParallel', false, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'Force', false, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'MaxObservation', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(parser, 'FormationSpacing', NaN, ...
    @(x) isnumeric(x) && isscalar(x) && (isnan(x) || x > 0));
addParameter(parser, 'TargetCounts', [], ...
    @(x) isnumeric(x) && isvector(x) && all(x >= 1));
addParameter(parser, 'ParallelBatchSize', NaN, ...
    @(x) isnumeric(x) && isscalar(x) && (isnan(x) || x >= 1));
addParameter(parser, 'WDGrid_mrad', [], ...
    @(x) isnumeric(x) && isvector(x) && all(x > 0) && all(x < pi * 1e3));
parse(parser, varargin{:});
options = parser.Results;
mode = lower(string(options.Mode));

moduleDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(moduleDir);
addpath(moduleDir, projectRoot);
cfg = configureMode(multiuav_default_config(), mode, options);
createOutputFolders(cfg);

if mode == "plotonly"
    output = plotExisting(cfg);
    return;
end

caseList = buildCaseList(cfg, mode);
results = cell(1, numel(caseList));
missing = zeros(1, numel(caseList));
missingCount = 0;
for caseIndex = 1:numel(caseList)
    checkpointFile = checkpointPath(cfg, caseList(caseIndex));
    if isfile(checkpointFile) && ~cfg.execution.force
        saved = load(checkpointFile, 'result');
        results{caseIndex} = saved.result;
    else
        missingCount = missingCount + 1;
        missing(missingCount) = caseIndex;
    end
end
missing = missing(1:missingCount);

if cfg.execution.useParallel && ~isempty(missing)
    batchSize = cfg.execution.parallelBatchSize;
    for batchStart = 1:batchSize:numel(missing)
        batchStop = min(batchStart + batchSize - 1, numel(missing));
        batchCases = missing(batchStart:batchStop);
        pendingResults = cell(size(batchCases));
        parfor pendingIndex = 1:numel(batchCases)
            caseIndex = batchCases(pendingIndex);
            pendingResults{pendingIndex} = runOneCase(cfg, caseList(caseIndex));
        end
        for pendingIndex = 1:numel(batchCases)
            caseIndex = batchCases(pendingIndex);
            result = pendingResults{pendingIndex}; %#ok<NASGU>
            results{caseIndex} = result;
            save(checkpointPath(cfg, caseList(caseIndex)), 'result', '-v7.3');
        end
        if cfg.execution.showProgress
            fprintf('[parallel %d/%d] checkpoints saved\n', ...
                batchStop, numel(missing));
        end
    end
else
    for pendingIndex = 1:numel(missing)
        caseIndex = missing(pendingIndex);
        if cfg.execution.showProgress
            fprintf('[%d/%d] %s\n', caseIndex, numel(caseList), ...
                describeCase(caseList(caseIndex)));
        end
        result = runOneCase(cfg, caseList(caseIndex)); %#ok<NASGU>
        results{caseIndex} = result;
        save(checkpointPath(cfg, caseList(caseIndex)), 'result', '-v7.3');
    end
end

[summary, comparisons] = multiuav_summarize(results, cfg);
optima = table();
optimumComparisons = table();
if mode == "formationopt"
    [optima, optimumComparisons] = ...
        multiuav_analyze_formation_optima(results, cfg);
    figureFiles = multiuav_plot_formation_optima( ...
        summary, optima, optimumComparisons, string(cfg.figuresDir), ...
        cfg.optimization.threshold_s);
else
    figureFiles = multiuav_plot_results(summary, comparisons, ...
        string(cfg.figuresDir));
end
resultFile = fullfile(cfg.resultsDir, sprintf('%s_results.mat', mode));
save(resultFile, 'results', 'summary', 'comparisons', 'cfg', ...
    'caseList', 'optima', 'optimumComparisons', '-v7.3');
writetable(summary, fullfile(cfg.resultsDir, ...
    sprintf('%s_summary.csv', mode)));
writetable(comparisons, fullfile(cfg.resultsDir, ...
    sprintf('%s_comparisons.csv', mode)));
if mode == "formationopt"
    writetable(optima, fullfile(cfg.resultsDir, ...
        'formationopt_optima.csv'));
    writetable(optimumComparisons, fullfile(cfg.resultsDir, ...
        'formationopt_optimum_comparisons.csv'));
end

output = struct('mode', mode, 'config', cfg, 'caseList', caseList, ...
    'results', {results}, 'summary', summary, ...
    'comparisons', comparisons, 'optima', optima, ...
    'optimumComparisons', optimumComparisons, ...
    'figureFiles', figureFiles, ...
    'resultFile', string(resultFile));
end

function cfg = configureMode(cfg, mode, options)
switch mode
    case 'smoke'
        cfg.runN = cfg.smokeN;
        cfg.scene.targetCounts = [1, 2, 4];
        cfg.scene.scenarioNames = ["async_maneuver", "sync_formation"];
        cfg.trajectories = "RA";
        cfg.system.maxRange_m = 300;
        cfg.system.minimumEntryRange_m = 100;
        cfg.system.prf_Hz = 250;
        cfg.system.noiseLutStep_m = 1;
        cfg.system.maxObservation_s = 3;
        cfg.scene.entrySpan_s = 0.5;
        cfg.scene.timelyThresholds_s = [1, 2, 3];
    case 'preview'
        cfg.runN = cfg.previewN;
        cfg.scene.targetCounts = [1, 2, 4, 8];
        cfg.scene.scenarioNames = ["async_maneuver", "sync_formation"];
        cfg.trajectories = ["RA", "SP"];
        cfg.system.maxObservation_s = 5;
        cfg.scene.entrySpan_s = 1;
        cfg.scene.timelyThresholds_s = [1, 2, 5];
    case 'screen'
        cfg.runN = cfg.screenN;
        cfg.scene.targetCounts = [1, 2, 4, 8, 16];
    case 'formal'
        cfg.runN = cfg.formalN;
    case 'formationopt'
        cfg.runN = cfg.previewN;
        cfg.scene.targetCounts = [4, 8];
        cfg.scene.scenarioNames = "sync_formation";
        cfg.trajectories = ["RA", "SP"];
        cfg.system.maxObservation_s = cfg.optimization.threshold_s;
        cfg.scene.entrySpan_s = 0;
        cfg.scene.timelyThresholds_s = cfg.optimization.threshold_s;
    case 'plotonly'
        cfg.runN = 0;
end
if isfinite(options.N)
    cfg.runN = round(options.N);
end
if strlength(string(options.OutputRoot)) > 0
    cfg.outputRoot = char(options.OutputRoot);
    cfg.resultsDir = fullfile(cfg.outputRoot, 'results');
    cfg.figuresDir = fullfile(cfg.outputRoot, 'figures');
    cfg.checkpointDir = fullfile(cfg.outputRoot, 'checkpoints');
end
if isfinite(options.MaxObservation)
    cfg.system.maxObservation_s = options.MaxObservation;
    cfg.scene.timelyThresholds_s = ...
        cfg.scene.timelyThresholds_s( ...
        cfg.scene.timelyThresholds_s <= options.MaxObservation);
end
if isfinite(options.FormationSpacing)
    cfg.scene.formationSpacing_m = options.FormationSpacing;
end
if ~isempty(options.TargetCounts)
    cfg.scene.targetCounts = unique(round(options.TargetCounts), 'stable');
end
if isfinite(options.ParallelBatchSize)
    cfg.execution.parallelBatchSize = round(options.ParallelBatchSize);
end
if ~isempty(options.WDGrid_mrad)
    cfg.manuscriptWD_mrad = unique(options.WDGrid_mrad, 'sorted');
end
cfg.execution.useParallel = options.UseParallel;
cfg.execution.force = options.Force;
cfg.constants.rangeResolution_m = cfg.constants.speedOfLight_m_s * ...
    cfg.system.rangeGate_s / 2;
end

function caseList = buildCaseList(cfg, mode)
caseList = struct('region', {}, 'trajectory', {}, 'scenario', {}, ...
    'beam', {}, 'targetCount', {});
if mode == "smoke"
    beamList = [makeBeam("point", 55), makeBeam("line", 125), ...
        makeBeam("ring", 125)];
elseif mode == "preview"
    beamList = [makeBeam("line", 125), makeBeam("line", 195), ...
        makeBeam("ring", 125)];
elseif mode == "formationopt"
    beamList = struct('name', {}, 'label', {}, 'wD_rad', {}, 'wd_rad', {});
    for width = cfg.manuscriptWD_mrad
        beamList(end + 1) = makeBeam("line", width); %#ok<AGROW>
        beamList(end + 1) = makeBeam("ring", width); %#ok<AGROW>
    end
else
    beamList = makeBeam("point", 55);
    widths = cfg.screenWD_mrad;
    if mode == "formal"
        widths = cfg.formalWD_mrad;
    end
    for width = widths
        beamList(end + 1) = makeBeam("line", width); %#ok<AGROW>
        beamList(end + 1) = makeBeam("ring", width); %#ok<AGROW>
    end
end

for targetCount = cfg.scene.targetCounts
    for scenario = cfg.scene.scenarioNames
        for trajectory = cfg.trajectories
            for beamIndex = 1:numel(beamList)
                caseList(end + 1) = struct('region', "ROI", ... %#ok<AGROW>
                    'trajectory', trajectory, 'scenario', scenario, ...
                    'beam', beamList(beamIndex), ...
                    'targetCount', targetCount);
            end
        end
    end
end

% ROE is deliberately restricted to the two reference beam widths and
% key scenarios; its role is verification rather than the main sweep.
if mode == "screen" || mode == "formal"
    roeBeams = [makeBeam("line", 195), makeBeam("ring", 125)];
    for targetCount = cfg.scene.targetCounts
        for scenario = ["async_maneuver", "sync_formation"]
            for beamIndex = 1:numel(roeBeams)
                caseList(end + 1) = struct('region', "ROE", ... %#ok<AGROW>
                    'trajectory', "SP", 'scenario', scenario, ...
                    'beam', roeBeams(beamIndex), ...
                    'targetCount', targetCount);
            end
        end
    end
end
end

function beam = makeBeam(name, width_mrad)
labels = struct('point', "Spot", 'line', "Line", 'ring', "Annular");
beam = struct('name', name, 'label', labels.(char(name)), ...
    'wD_rad', width_mrad * 1e-3, 'wd_rad', 1e-3);
if name == "point"
    beam.wd_rad = NaN;
end
end

function result = runOneCase(cfg, caseDefinition)
bank = multiuav_build_scenario_bank(cfg.runN, ...
    caseDefinition.targetCount, cfg.seed + caseDefinition.targetCount);
result = multiuav_simulate_case(cfg, caseDefinition.beam, ...
    caseDefinition.region, caseDefinition.trajectory, ...
    caseDefinition.scenario, bank);
end

function path = checkpointPath(cfg, caseDefinition)
name = sprintf('%s_%s_%s_%s_%gmrad_M%d_N%d_seed%d.mat', ...
    caseDefinition.region, caseDefinition.trajectory, ...
    caseDefinition.scenario, caseDefinition.beam.name, ...
    caseDefinition.beam.wD_rad * 1e3, caseDefinition.targetCount, ...
    cfg.runN, cfg.seed);
if caseDefinition.scenario == "sync_formation" && ...
        abs(cfg.scene.formationSpacing_m - 15) > 1e-9
    name = replace(name, '.mat', sprintf('_spacing%gm.mat', ...
        cfg.scene.formationSpacing_m));
end
path = fullfile(cfg.checkpointDir, regexprep(name, '[^A-Za-z0-9_.-]', '_'));
end

function text = describeCase(caseDefinition)
text = sprintf('%s/%s/%s/%s/%g mrad/M=%d', ...
    caseDefinition.region, caseDefinition.trajectory, ...
    caseDefinition.scenario, caseDefinition.beam.name, ...
    caseDefinition.beam.wD_rad * 1e3, caseDefinition.targetCount);
end

function createOutputFolders(cfg)
folders = {cfg.outputRoot, cfg.resultsDir, cfg.figuresDir, cfg.checkpointDir};
for index = 1:numel(folders)
    if ~isfolder(folders{index})
        mkdir(folders{index});
    end
end
end

function output = plotExisting(cfg)
files = dir(fullfile(cfg.resultsDir, '*_results.mat'));
if isempty(files)
    error('run_multiuav_analysis:NoSavedResults', ...
        'No saved result file was found for plotonly mode.');
end
[~, latestIndex] = max([files.datenum]);
saved = load(fullfile(files(latestIndex).folder, files(latestIndex).name), ...
    'summary', 'comparisons', 'optima', 'optimumComparisons', 'cfg');
if isfield(saved, 'optima') && ~isempty(saved.optima)
    figureFiles = multiuav_plot_formation_optima(saved.summary, ...
        saved.optima, saved.optimumComparisons, string(cfg.figuresDir), ...
        saved.cfg.optimization.threshold_s);
else
    figureFiles = multiuav_plot_results(saved.summary, saved.comparisons, ...
        string(cfg.figuresDir));
end
output = struct('mode', "plotonly", 'summary', saved.summary, ...
    'comparisons', saved.comparisons, 'figureFiles', figureFiles);
end

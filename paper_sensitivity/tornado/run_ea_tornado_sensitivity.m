function outputs = run_ea_tornado_sensitivity(options)
%RUN_EA_TORNADO_SENSITIVITY Run the paired EA tornado workflow.
%
%   outputs = run_ea_tornado_sensitivity
%   outputs = run_ea_tornado_sensitivity(Mode="simulate")
%   outputs = run_ea_tornado_sensitivity(Mode="analyze")
%   outputs = run_ea_tornado_sensitivity(Mode="plot")
%   outputs = run_ea_tornado_sensitivity(Mode="smoke", ...
%       OutputDirectory=tempname, UseParallel=false)

arguments
    options.Mode (1,1) string {mustBeMember(options.Mode, ...
        ["full", "simulate", "analyze", "plot", "smoke"])} = "full"
    options.N (1,1) double {mustBeInteger,mustBePositive} = 3000
    options.BootstrapReplicates (1,1) double ...
        {mustBeInteger,mustBePositive} = 1000
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 20260928
    options.BootstrapSeed (1,1) double ...
        {mustBeInteger,mustBeNonnegative} = 20260929
    options.UseParallel (1,1) logical = true
    options.Workers (1,1) double {mustBeInteger,mustBeNonnegative} = 4
    options.OutputDirectory (1,1) string = ""
    options.BeamPathDirectory (1,1) string = "G:/BeamVEC_NEW"
    options.Visible (1,1) string ...
        {mustBeMember(options.Visible, ["on", "off"])} = "off"
    options.InitialWidthsMrad (1,:) double ...
        {mustBeFinite,mustBePositive} = 5:5:150
    options.MaxWidthMrad (1,1) double ...
        {mustBeFinite,mustBePositive} = 400
    options.TaskLimit (1,1) double {mustBeNonnegative} = inf
end

moduleDirectory = fileparts(mfilename("fullpath"));
projectRoot = fileparts(fileparts(moduleDirectory));
addpath(projectRoot, fileparts(moduleDirectory), moduleDirectory);

isSmoke = options.Mode == "smoke";
if isSmoke
    options.N = min(options.N, 2);
    options.BootstrapReplicates = min(options.BootstrapReplicates, 20);
    options.UseParallel = false;
    options.Workers = 0;
    options.InitialWidthsMrad = [50, 100, 150];
    options.MaxWidthMrad = 150;
end
if options.OutputDirectory == ""
    if isSmoke
        options.OutputDirectory = fullfile(moduleDirectory, "results", "smoke");
    else
        options.OutputDirectory = fullfile(moduleDirectory, "results", "formal");
    end
end

cfg = makeConfiguration(options, moduleDirectory, projectRoot, isSmoke);
ensureFolders(cfg);
outputs = struct("config", cfg);

if options.Mode == "plot"
    outputs.figures = ea_tornado_plot(cfg.outputDirectory, options.Visible);
    return;
end

if ismember(options.Mode, ["full", "simulate", "smoke"])
    design = simulateAdaptive(cfg, options.TaskLimit);
    outputs.design = design;
    if options.Mode == "simulate" || ~all(design.complete(design.required))
        return;
    end
end

if ismember(options.Mode, ["full", "analyze", "smoke"])
    analysis = ea_tornado_analyze(cfg.outputDirectory, ...
        cfg.bootstrapReplicates, cfg.bootstrapSeed);
    outputs.analysis = analysis;
end

if ismember(options.Mode, ["full", "smoke"])
    outputs.figures = ea_tornado_plot(cfg.outputDirectory, options.Visible);
    outputs.report = ea_tornado_write_report(cfg.outputDirectory);
    outputs.verification = ea_tornado_verify(cfg.outputDirectory);
end
end

function cfg = makeConfiguration(options, moduleDirectory, projectRoot, isSmoke)
assert(issorted(options.InitialWidthsMrad, "strictascend") && ...
    numel(unique(options.InitialWidthsMrad)) == ...
    numel(options.InitialWidthsMrad), ...
    "EATornado:InitialGrid", ...
    "Initial widths must be unique and strictly increasing.");
assert(options.MaxWidthMrad >= max(options.InitialWidthsMrad), ...
    "EATornado:MaximumWidth", ...
    "MaxWidthMrad must cover InitialWidthsMrad.");
assert(all(abs(5 * round(options.InitialWidthsMrad / 5) - ...
    options.InitialWidthsMrad) < 1e-10), ...
    "EATornado:ScanGrid", "EA trajectories use a 5 mrad grid.");

cfg = struct();
cfg.schemaVersion = 1;
cfg.moduleDirectory = string(moduleDirectory);
cfg.projectRoot = string(projectRoot);
cfg.outputDirectory = string(options.OutputDirectory);
cfg.rawDirectory = fullfile(cfg.outputDirectory, "raw");
cfg.figuresDirectory = fullfile(cfg.outputDirectory, "figures");
cfg.N = options.N;
cfg.bootstrapReplicates = options.BootstrapReplicates;
cfg.seed = options.Seed;
cfg.bootstrapSeed = options.BootstrapSeed;
cfg.useParallel = options.UseParallel;
cfg.workers = options.Workers;
cfg.beamPathDirectory = string(options.BeamPathDirectory);
cfg.initialWidthsMrad = options.InitialWidthsMrad;
cfg.maxWidthMrad = options.MaxWidthMrad;
cfg.factorScales = [0.8, 0.9, 1.0, 1.1, 1.2];
cfg.factorNames = ["energy", "reflectivity", "speed", "alpha"];
cfg.factorLabels = ["Pulse energy", "Target reflectivity", ...
    "Target speed", "Atmospheric extinction"];
cfg.factorUnits = ["uJ", "1", "m/s", "1e-5 /m"];
cfg.baseline = [150e-6, 0.5, 30, 1.5e-5];
cfg.displayScale = [1e6, 1, 1, 1e5];
cfg.alphaBetaRatio_sr = 50;
cfg.deadline_s = 15;
cfg.requiredProbability = 0.60;
cfg.feasibilityThreshold = 0.975;
cfg.f_Hz = 5000;
cfg.omega_rad_s = 2 * pi;
cfg.maxRange_m = 2000;
cfg.smallWidth_rad = 1e-3;
cfg.isSmoke = isSmoke;
if isSmoke
    cfg.deadline_s = 0.2;
    cfg.requiredProbability = 0.50;
end

cfg.extensionBreaksMrad = unique([max(cfg.initialWidthsMrad), ...
    min(250, cfg.maxWidthMrad), cfg.maxWidthMrad], "stable");
cfg.extensionBreaksMrad = cfg.extensionBreaksMrad( ...
    cfg.extensionBreaksMrad >= max(cfg.initialWidthsMrad));
extraWidths = (max(cfg.initialWidthsMrad) + 5):5:cfg.maxWidthMrad;
cfg.allWidthsMrad = unique([cfg.initialWidthsMrad, extraWidths]);
cfg.matlabRelease = string(version("-release"));
end

function ensureFolders(cfg)
folders = [cfg.outputDirectory, cfg.rawDirectory, cfg.figuresDirectory];
for folder = folders
    if ~isfolder(folder)
        mkdir(folder);
    end
end
end

function design = simulateAdaptive(cfg, taskLimit)
manifestFile = fullfile(cfg.outputDirectory, "design.mat");
sourceSummary = buildSourceSummary(cfg);
signature = simulationSignature(cfg, sourceSummary);

if isfile(manifestFile)
    saved = load(manifestFile, "design");
    design = saved.design;
    assert(isequaln(design.signature, signature), ...
        "EATornado:CacheSignature", ...
        "Configuration or source changed; choose a new OutputDirectory.");
else
    design = createDesign(cfg, signature, sourceSummary);
    saveDesign(manifestFile, design);
end

if taskLimit == 0
    return;
end
assert(exist("readNPY", "file") == 2, "EATornado:ReadNPY", ...
    "readNPY must be available on the MATLAB path.");

stageBreaks = cfg.extensionBreaksMrad;
completedTasks = 0;
candidateCases = true(1, numel(design.cases));
for stageIndex = 1:numel(stageBreaks)
    stageMaximum = stageBreaks(stageIndex);
    if stageIndex == 1
        stageWidths = cfg.initialWidthsMrad;
    else
        previousMaximum = stageBreaks(stageIndex - 1);
        stageWidths = (previousMaximum + 5):5:stageMaximum;
    end
    widthIndices = find(ismember(design.widthsMrad, stageWidths));
    design.required(widthIndices, candidateCases) = true;
    saveDesign(manifestFile, design);

    [design, completedTasks, stopped] = simulateRequiredTasks( ...
        cfg, design, manifestFile, completedTasks, taskLimit);
    if stopped || stageIndex == numel(stageBreaks)
        return;
    end

    candidateCases = extensionCases(cfg, design, stageMaximum, candidateCases);
    if ~any(candidateCases)
        return;
    end
end
end

function design = createDesign(cfg, signature, sourceSummary)
previousRng = rng;
cleanup = onCleanup(@() rng(previousRng)); %#ok<NASGU>
rng(cfg.seed, "twister");
[position, velocity] = init_UAV(cfg.maxRange_m, 1, 30, cfg.N);
phase = rand(cfg.N, 1);

cases = buildCases(cfg);
nWidths = numel(cfg.allWidthsMrad);
design = struct();
design.schemaVersion = cfg.schemaVersion;
design.signature = signature;
design.sourceSummary = sourceSummary;
design.config = cfg;
design.cases = cases;
design.widthsMrad = cfg.allWidthsMrad;
design.position_m = position;
design.velocityBaseline_m_s = velocity;
design.scanPhase_u = phase;
design.required = false(nWidths, numel(cases));
design.complete = false(nWidths, numel(cases));
initialIndices = ismember(design.widthsMrad, cfg.initialWidthsMrad);
design.required(initialIndices, :) = true;
design.createdAt = string(datetime("now", ...
    "Format", "yyyy-MM-dd HH:mm:ss Z"));
end

function cases = buildCases(cfg)
prototype = struct("id", "", "factorIndex", 0, "scale", 1, ...
    "values", cfg.baseline, "physics", zeros(1, 4), "speed_m_s", 30);
cases = prototype;
cases(1).id = "baseline";
cases(1).physics = [cfg.baseline(1:2), cfg.baseline(4), ...
    cfg.baseline(4) / cfg.alphaBetaRatio_sr];
writeIndex = 1;
for factorIndex = 1:numel(cfg.factorNames)
    for scale = cfg.factorScales(cfg.factorScales ~= 1)
        writeIndex = writeIndex + 1;
        item = prototype;
        item.factorIndex = factorIndex;
        item.scale = scale;
        item.values = cfg.baseline;
        item.values(factorIndex) = item.values(factorIndex) * scale;
        item.speed_m_s = item.values(3);
        item.physics = [item.values(1:2), item.values(4), ...
            item.values(4) / cfg.alphaBetaRatio_sr];
        item.id = cfg.factorNames(factorIndex) + "_" + scaleCode(scale);
        cases(writeIndex) = item; %#ok<AGROW>
    end
end
end

function code = scaleCode(scale)
percent = round(100 * (scale - 1));
if percent < 0
    code = "m" + abs(percent);
else
    code = "p" + percent;
end
end

function [design, completedTasks, stopped] = simulateRequiredTasks( ...
        cfg, design, manifestFile, completedTasks, taskLimit)
stopped = false;
for widthIndex = 1:numel(design.widthsMrad)
    caseIndices = find(design.required(widthIndex, :) & ...
        ~design.complete(widthIndex, :));
    if isempty(caseIndices)
        continue;
    end
    widthMrad = design.widthsMrad(widthIndex);
    scan = loadScan(cfg, widthMrad);
    initialIndex = floor(design.scanPhase_u * scan.phaseLength) + 1;
    for caseIndex = caseIndices
        if completedTasks >= taskLimit
            stopped = true;
            return;
        end
        rawFile = rawResultFile(cfg, design.cases(caseIndex).id, widthMrad);
        if isfile(rawFile)
            cached = load(rawFile, "meta");
            validateRawMetadata(cached.meta, design, caseIndex, widthMrad);
            design.complete(widthIndex, caseIndex) = true;
            continue;
        end

        item = design.cases(caseIndex);
        velocity = design.velocityBaseline_m_s * (item.speed_m_s / 30);
        hit = false(cfg.N, 2);
        time_s = nan(cfg.N, 2);
        for beamIndex = 1:2
            if beamIndex == 1
                kernel = @MC_line_snr;
                beamName = "line";
            else
                kernel = @MC_ring_snr;
                beamName = "ring";
            end
            blockSize = mc_ea_block_size(widthMrad * 1e-3, beamName);
            [hit(:, beamIndex), time_s(:, beamIndex)] = simulateBank( ...
                cfg, kernel, scan, design.position_m, velocity, ...
                initialIndex, widthMrad, blockSize, item.physics);
        end
        meta = struct("SchemaVersion", cfg.schemaVersion, ...
            "CaseIndex", caseIndex, "CaseId", item.id, ...
            "WidthMrad", widthMrad, "N", cfg.N, "Seed", cfg.seed, ...
            "Deadline_s", cfg.deadline_s, "Physics", item.physics, ...
            "Speed_m_s", item.speed_m_s, ...
            "Signature", design.signature);
        if ~isfolder(fileparts(rawFile))
            mkdir(fileparts(rawFile));
        end
        temporaryFile = rawFile + ".partial.mat";
        save(temporaryFile, "meta", "hit", "time_s", "-v7.3");
        movefile(temporaryFile, rawFile, "f");
        design.complete(widthIndex, caseIndex) = true;
        completedTasks = completedTasks + 1;
        saveDesign(manifestFile, design);
        fprintf("EA TORNADO %s, w_D=%g mrad: P_D=%s (%d tasks)\n", ...
            item.id, widthMrad, ...
            mat2str(mean(hit & time_s <= cfg.deadline_s, 1), 3), ...
            completedTasks);
    end
end
saveDesign(manifestFile, design);
end

function scan = loadScan(cfg, widthMrad)
scanFile = fullfile(cfg.beamPathDirectory, sprintf( ...
    "beam_vector_fast_5000Hz_h3000m_D%gu.npy", widthMrad * 1000));
assert(isfile(scanFile), "EATornado:MissingScan", ...
    "Missing EA scan trajectory: %s", scanFile);
scan.beam = readNPY(char(scanFile));
assert(isnumeric(scan.beam) && ismatrix(scan.beam) && ...
    size(scan.beam, 2) == 3 && ~isempty(scan.beam) && ...
    all(isfinite(scan.beam), "all"), ...
    "EATornado:InvalidScan", "EA scan trajectory must be finite N-by-3 data.");
scan.path1 = boundaryCircle(scan.beam(1, :), cfg.f_Hz, cfg.omega_rad_s);
scan.path2 = boundaryCircle(scan.beam(end, :), cfg.f_Hz, cfg.omega_rad_s);
scan.phaseLength = size(scan.path1, 1) + size(scan.beam, 1) + ...
    size(scan.path2, 1);
end

function [hit, time_s] = simulateBank(cfg, kernel, scan, position, ...
        velocity, initialIndex, widthMrad, blockSize, physics)
n = size(position, 1);
hit = false(n, 1);
time_s = nan(n, 1);
maximumStep = floor(cfg.deadline_s * cfg.f_Hz);
f = cfg.f_Hz;
omega = cfg.omega_rad_s;
smallWidth = cfg.smallWidth_rad;
range = cfg.maxRange_m;
beam = scan.beam;
path1 = scan.path1;
path2 = scan.path2;

if cfg.useParallel
    pool = gcp("nocreate");
    if isempty(pool)
        if cfg.workers > 0
            pool = parpool("threads", cfg.workers);
        else
            pool = parpool("threads");
        end
    end
    workerCount = pool.NumWorkers;
    parfor (targetIndex = 1:n, workerCount)
        [detected, detectionStep] = kernel(widthMrad * 1e-3, ...
            smallWidth, path1, beam, path2, f, omega, ...
            position(targetIndex, :), velocity(targetIndex, :), ...
            initialIndex(targetIndex), range, blockSize, maximumStep, ...
            "fixed", physics);
        hit(targetIndex) = logical(detected);
        time_s(targetIndex) = detectionStep / f;
    end
else
    for targetIndex = 1:n
        [detected, detectionStep] = kernel(widthMrad * 1e-3, ...
            smallWidth, path1, beam, path2, f, omega, ...
            position(targetIndex, :), velocity(targetIndex, :), ...
            initialIndex(targetIndex), range, blockSize, maximumStep, ...
            "fixed", physics);
        hit(targetIndex) = logical(detected);
        time_s(targetIndex) = detectionStep / f;
    end
end
end

function selected = extensionCases(cfg, design, stageMaximum, candidates)
selected = false(size(candidates));
weights = bootstrapWeights(cfg.N, cfg.bootstrapReplicates, ...
    cfg.bootstrapSeed);
for caseIndex = find(candidates)
    widthMask = design.required(:, caseIndex) & ...
        design.widthsMrad(:) <= stageMaximum;
    widths = design.widthsMrad(widthMask).';
    events = loadCaseEvents(cfg, design, caseIndex, widthMask);
    probability = squeeze(mean(events, 1));
    minimum = ea_tornado_min_width(widths, probability, ...
        cfg.requiredProbability);
    feasible = bootstrapFeasibility(widths, events, weights, ...
        cfg.requiredProbability);
    selected(caseIndex) = any(~isfinite(minimum)) || ...
        feasible < cfg.feasibilityThreshold;
end
end

function events = loadCaseEvents(cfg, design, caseIndex, widthMask)
widths = design.widthsMrad(widthMask);
events = false(cfg.N, numel(widths), 2);
for widthIndex = 1:numel(widths)
    file = rawResultFile(cfg, design.cases(caseIndex).id, ...
        widths(widthIndex));
    raw = load(file, "meta", "hit", "time_s");
    validateRawMetadata(raw.meta, design, caseIndex, widths(widthIndex));
    events(:, widthIndex, :) = reshape(raw.hit & ...
        raw.time_s <= cfg.deadline_s, cfg.N, 1, 2);
end
end

function feasibility = bootstrapFeasibility(widths, events, weights, target)
n = size(events, 1);
replicates = size(weights, 2);
feasible = false(replicates, 1);
flatEvents = reshape(events, n, []);
boot = flatEvents.' * weights / n;
boot = reshape(boot, numel(widths), 2, replicates);
for replicate = 1:replicates
    minimum = ea_tornado_min_width(widths, ...
        boot(:, :, replicate), target);
    feasible(replicate) = all(isfinite(minimum));
end
feasibility = mean(feasible);
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

function file = rawResultFile(cfg, caseId, widthMrad)
file = fullfile(cfg.rawDirectory, caseId, sprintf( ...
    "ea_tornado_width_%03d.mat", round(widthMrad)));
end

function validateRawMetadata(meta, design, caseIndex, widthMrad)
item = design.cases(caseIndex);
assert(meta.CaseIndex == caseIndex && meta.CaseId == item.id && ...
    meta.WidthMrad == widthMrad && meta.N == design.config.N && ...
    meta.Seed == design.config.seed && ...
    isequaln(meta.Physics, item.physics) && ...
    meta.Speed_m_s == item.speed_m_s && ...
    isequaln(meta.Signature, design.signature), ...
    "EATornado:RawCache", "Incompatible raw checkpoint.");
end

function circle = boundaryCircle(point, f, omega)
theta0 = atan2(point(2), point(1));
radius = hypot(point(1), point(2));
theta = (theta0:omega / f / radius:(theta0 + 2 * pi)).';
circle = [radius * cos(theta), radius * sin(theta), ...
    ones(numel(theta), 1) * point(3)];
end

function saveDesign(file, design)
temporaryFile = file + ".partial.mat";
save(temporaryFile, "design", "-v7.3");
movefile(temporaryFile, file, "f");
end

function sourceSummary = buildSourceSummary(cfg)
names = ["run_ea_tornado_sensitivity.m", ...
    "ea_tornado_analyze.m", "ea_tornado_min_width.m", ...
    "MC_line_snr.m", "MC_ring_snr.m", "init_UAV.m", ...
    "append_and_check_photon_window.m", ...
    "empty_detection_window_photons.m", "mc_ea_block_size.m"];
paths = strings(size(names));
paths(1:3) = fullfile(cfg.moduleDirectory, names(1:3));
paths(4:end) = fullfile(cfg.projectRoot, names(4:end));
bytes = zeros(size(names));
modified = strings(size(names));
sha256 = strings(size(names));
for sourceIndex = 1:numel(names)
    info = dir(paths(sourceIndex));
    assert(~isempty(info), "EATornado:MissingSource", ...
        "Missing source file: %s", paths(sourceIndex));
    bytes(sourceIndex) = info.bytes;
    modified(sourceIndex) = string(datetime(info.datenum, ...
        "ConvertFrom", "datenum", "Format", "yyyy-MM-dd HH:mm:ss"));
    sha256(sourceIndex) = fileSha256(paths(sourceIndex));
end
sourceSummary = table();
sourceSummary.Name = names(:);
sourceSummary.Path = paths(:);
sourceSummary.Bytes = bytes(:);
sourceSummary.Modified = modified(:);
sourceSummary.SHA256 = sha256(:);
end

function signature = simulationSignature(cfg, sourceSummary)
signature = struct();
signature.schemaVersion = cfg.schemaVersion;
signature.N = cfg.N;
signature.seed = cfg.seed;
signature.factorScales = cfg.factorScales;
signature.baseline = cfg.baseline;
signature.deadline_s = cfg.deadline_s;
signature.requiredProbability = cfg.requiredProbability;
signature.initialWidthsMrad = cfg.initialWidthsMrad;
signature.maxWidthMrad = cfg.maxWidthMrad;
signature.beamPathDirectory = cfg.beamPathDirectory;
signature.sourceSHA256 = sourceSummary.SHA256;
end

function digest = fileSha256(file)
bytes = uint8(fileread(file));
engine = java.security.MessageDigest.getInstance("SHA-256");
engine.update(bytes);
raw = typecast(engine.digest(), "uint8");
digest = lower(string(reshape(dec2hex(raw, 2).', 1, [])));
end

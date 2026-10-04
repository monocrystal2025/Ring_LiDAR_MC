function result = multiuav_simulate_case(cfg, beam, region, trajectory, ...
    scenario, bank)
%MULTIUAV_SIMULATE_CASE Simulate one paired scene-level multi-UAV case.

arguments
    cfg (1,1) struct
    beam (1,1) struct
    region (1,1) string
    trajectory (1,1) string
    scenario (1,1) string
    bank (1,1) struct
end

validateInputs(cfg, beam, region, trajectory, scenario, bank);
region = upper(region);
trajectory = upper(trajectory);
scenario = lower(scenario);
system = cfg.system;
N = bank.N;
M = bank.maxTargets;

scan = multiuav_build_scan(cfg, beam, region, trajectory);
photon = buildPhotonModel(cfg, beam);
[position, velocity, entryTime] = mapSceneBank( ...
    cfg, bank, region, scenario);
initialPosition = position;
initialPhase = floor(bank.scanPhase_u .* scan.cycleLength) + 1;
initialPhase = min(initialPhase, scan.cycleLength);

detected = false(N, M);
exited = false(N, M);
firstTime_s = nan(N, M);
encounterCount = zeros(N, M);
resolvedHitCount = zeros(N, M);
ambiguousHitCount = zeros(N, M);
previousCandidate = false(N, M);
maxResolvedPerPulse = zeros(N, 1);
multiTargetPulseCount = zeros(N, 1);
maxInterceptedEnergyFraction = zeros(N, 1);

windowPulses = photon.windowPulses;
scoreBuffer = zeros(N, M, windowPulses);
rollingScore = zeros(N, M);
bufferIndex = 0;
omegaState = zeros(N, M, 3);
turnStream = RandStream('mt19937ar', 'Seed', bank.seed + 7919);

dt = 1 / system.prf_Hz;
motionPulses = max(1, round(cfg.execution.motionUpdate_s * system.prf_Hz));
lastEntry_s = max(entryTime, [], 'all');
lastPulse = ceil((lastEntry_s + system.maxObservation_s) * system.prf_Hz);
timerValue = tic;

for pulseIndex = 1:lastPulse
    time_s = (pulseIndex - 1) * dt;
    entered = time_s + 10 * eps(time_s + 1) >= entryTime;
    inside = isInsideRegion(position, system, region);
    exited = exited | (entered & ~inside & ~detected);
    active = entered & ~exited & ~detected;

    scanIndex = mod(initialPhase + pulseIndex - 2, scan.cycleLength) + 1;
    B = reshape(scan.B(scanIndex, :), N, 1, 3);
    ranges = sqrt(sum(position.^2, 3));
    safeRanges = max(ranges, eps);
    directions = position ./ safeRanges;
    illumination = illuminationFactor( ...
        cfg, beam, directions, ranges, B, scan, scanIndex);
    illumination(~active) = 0;
    candidate = illumination > 0;

    [resolved, ambiguous] = multiuav_apply_range_resolution( ...
        ranges, candidate, cfg.constants.rangeResolution_m);
    [signal, interceptedFraction] = calcSignal( ...
        photon, ranges, illumination);
    noise = lookupNoise(photon, ranges);
    signal(~resolved) = 0;
    score = zeros(N, M);
    validScore = resolved & signal > 0 & signal + noise > 0;
    score(validScore) = signal(validScore).^2 ./ ...
        (signal(validScore) + noise(validScore));

    bufferIndex = mod(bufferIndex, windowPulses) + 1;
    rollingScore = rollingScore - scoreBuffer(:, :, bufferIndex) + score;
    scoreBuffer(:, :, bufferIndex) = score;
    newDetection = active & ~detected & ...
        sqrt(max(rollingScore, 0)) >= system.snrThreshold;
    firstTime_s(newDetection) = time_s - entryTime(newDetection);
    detected = detected | newDetection;

    encounterCount = encounterCount + double(candidate & ~previousCandidate);
    resolvedHitCount = resolvedHitCount + double(resolved);
    ambiguousHitCount = ambiguousHitCount + double(ambiguous);
    previousCandidate = candidate;

    simultaneousCount = sum(resolved, 2);
    maxResolvedPerPulse = max(maxResolvedPerPulse, simultaneousCount);
    multiTargetPulseCount = multiTargetPulseCount + ...
        double(simultaneousCount >= 2);
    maxInterceptedEnergyFraction = max(maxInterceptedEnergyFraction, ...
        sum(interceptedFraction, 2));

    allFinished = all(detected | exited, 2);
    if all(allFinished) && time_s >= lastEntry_s
        break;
    end

    moveMask = entered & ~exited & ~detected;
    moveMask3 = reshape(moveMask, N, M, 1);
    position = position + velocity .* dt .* moveMask3;
    if scenario == "async_maneuver" && mod(pulseIndex, motionPulses) == 0
        [velocity, omegaState] = updateRandomTurn( ...
            velocity, omegaState, moveMask, cfg, turnStream);
    end
end

if any(maxInterceptedEnergyFraction > 1 + 1e-9)
    error('multiuav_simulate_case:EnergyFractionExceeded', ...
        'The summed intercepted target-energy fraction exceeded unity.');
end

result = struct();
result.caseId = makeCaseId(region, trajectory, scenario, beam, M);
result.region = region;
result.trajectory = trajectory;
result.scenario = scenario;
result.beam = string(beam.name);
result.beamLabel = string(beam.label);
result.wD_mrad = beam.wD_rad * 1e3;
result.wd_mrad = beam.wd_rad * 1e3;
result.scenarioId = bank.scenarioId;
result.targetId = bank.targetId;
result.targetCount = M;
result.detect = detected;
result.exitWithoutDetection = exited & ~detected;
result.firstTime_s = firstTime_s;
result.entryTime_s = entryTime;
result.encounterCount = encounterCount;
result.resolvedHitCount = resolvedHitCount;
result.ambiguousHitCount = ambiguousHitCount;
result.maxResolvedPerPulse = maxResolvedPerPulse;
result.multiTargetPulseCount = multiTargetPulseCount;
result.maxInterceptedEnergyFraction = maxInterceptedEnergyFraction;
result.initialScanPhase_u = bank.scanPhase_u;
result.initialPosition_m = initialPosition;
result.finalPosition_m = position;
result.finalVelocity_m_s = velocity;
result.cycleLength = scan.cycleLength;
result.cycleTime_s = scan.cycleTime_s;
result.windowPulses = windowPulses;
result.rangeResolution_m = cfg.constants.rangeResolution_m;
result.simulatedPulses = pulseIndex;
result.elapsed_s = toc(timerValue);
end

function validateInputs(cfg, beam, region, trajectory, scenario, bank)
requiredCfg = {'system', 'scene', 'execution', 'constants'};
requiredBank = {'N', 'maxTargets', 'seed', 'scenarioId', 'targetId', ...
    'entryAzimuth_u', 'entryRange_u', 'velocityAzimuth_u', ...
    'velocityCosPolar_u', 'entryTime_u', 'scanPhase_u'};
if ~all(isfield(cfg, requiredCfg))
    error('multiuav_simulate_case:InvalidConfig', ...
        'cfg is missing a required section.');
end
if ~all(isfield(bank, requiredBank))
    error('multiuav_simulate_case:InvalidScenarioBank', ...
        'bank is missing required paired draws.');
end
if ~any(strcmpi(string(beam.name), ["point", "line", "ring"]))
    error('multiuav_simulate_case:UnknownBeam', ...
        'beam.name must be point, line, or ring.');
end
if ~any(upper(region) == ["ROI", "ROE"])
    error('multiuav_simulate_case:UnknownRegion', ...
        'region must be ROI or ROE.');
end
if upper(region) == "ROI" && ~any(upper(trajectory) == ["RA", "SP"])
    error('multiuav_simulate_case:UnknownTrajectory', ...
        'ROI trajectory must be RA or SP.');
end
if ~any(lower(scenario) == ["async_maneuver", ...
        "sync_dispersed", "sync_formation"])
    error('multiuav_simulate_case:UnknownScenario', ...
        'Unknown multi-UAV scenario.');
end
if bank.N ~= size(bank.entryAzimuth_u, 1) || ...
        bank.maxTargets ~= size(bank.entryAzimuth_u, 2)
    error('multiuav_simulate_case:ScenarioCountMismatch', ...
        'bank dimensions are inconsistent.');
end
end

function [P, V, entryTime] = mapSceneBank(cfg, bank, region, scenario)
system = cfg.system;
N = bank.N;
M = bank.maxTargets;
azimuth = 2 * pi .* bank.entryAzimuth_u;

if region == "ROI"
    alpha = deg2rad(system.roiHalfAngle_deg);
    slantRange = system.minimumEntryRange_m + ...
        (system.maxRange_m - system.minimumEntryRange_m) .* bank.entryRange_u;
    P = cat(3, slantRange .* sin(alpha) .* cos(azimuth), ...
        slantRange .* sin(alpha) .* sin(azimuth), ...
        slantRange .* cos(alpha));
else
    cosPolar = bank.entryRange_u;
    sinPolar = sqrt(max(1 - cosPolar.^2, 0));
    P = system.maxRange_m .* cat(3, sinPolar .* cos(azimuth), ...
        sinPolar .* sin(azimuth), cosPolar);
end

velocityAzimuth = 2 * pi .* bank.velocityAzimuth_u;
velocityCosPolar = 2 .* bank.velocityCosPolar_u - 1;
velocitySinPolar = sqrt(max(1 - velocityCosPolar.^2, 0));
direction = cat(3, velocitySinPolar .* cos(velocityAzimuth), ...
    velocitySinPolar .* sin(velocityAzimuth), velocityCosPolar);

if scenario == "sync_formation"
    [P, direction] = makeFormation(P, cfg, bank, region);
end

direction = forceInward(P, direction, cfg, region);
V = system.targetSpeed_m_s .* direction;
if scenario == "async_maneuver"
    entryTime = cfg.scene.entrySpan_s .* bank.entryTime_u;
else
    entryTime = zeros(N, M);
end
end

function [P, direction] = makeFormation(P, cfg, bank, region)
N = bank.N;
M = bank.maxTargets;
centroid = P(:, 1, :);
azimuth = atan2(centroid(:, 1, 2), centroid(:, 1, 1));
eAz = cat(3, -sin(azimuth), cos(azimuth), zeros(N, 1));
eGenerator = centroid ./ max(sqrt(sum(centroid.^2, 3)), eps);

% A centre-first phyllotaxis layout keeps every prefix meaningful: the
% M=1 subset is the formation centroid and larger target counts expand
% outwards without changing the positions already assigned.
memberIndex = 0:M-1;
goldenAngle = pi * (3 - sqrt(5));
gridRadius = sqrt(memberIndex);
gridX = gridRadius .* cos(memberIndex .* goldenAngle);
gridY = gridRadius .* sin(memberIndex .* goldenAngle);
rotation = 2 * pi .* bank.formationRotation_u;
offsetX = reshape(cos(rotation) .* gridX - sin(rotation) .* gridY, ...
    N, M, 1);
offsetY = reshape(sin(rotation) .* gridX + cos(rotation) .* gridY, ...
    N, M, 1);
spacing = cfg.scene.formationSpacing_m;
P = centroid + spacing .* offsetX .* eAz + ...
    spacing .* offsetY .* eGenerator;

% Put every formation member exactly back on the region entry surface.
% The tangent-plane construction above otherwise leaves small members just
% outside a conical ROI because of surface curvature.
if region == "ROI"
    targetAzimuth = atan2(P(:, :, 2), P(:, :, 1));
    targetRange = sqrt(sum(P.^2, 3));
    targetRange = min(targetRange, cfg.system.maxRange_m);
    alpha = deg2rad(cfg.system.roiHalfAngle_deg);
    P = cat(3, targetRange .* sin(alpha) .* cos(targetAzimuth), ...
        targetRange .* sin(alpha) .* sin(targetAzimuth), ...
        targetRange .* cos(alpha));
else
    P = cfg.system.maxRange_m .* P ./ ...
        max(sqrt(sum(P.^2, 3)), eps);
end

centroidDirection = -centroid ./ max(sqrt(sum(centroid.^2, 3)), eps);
if region == "ROI"
    alpha = deg2rad(cfg.system.roiHalfAngle_deg);
    normal = cat(3, cos(azimuth), sin(azimuth), ...
        -tan(alpha) .* ones(N, 1));
    normal = normal ./ sqrt(sum(normal.^2, 3));
    centroidDirection = centroidDirection - 0.2 .* normal;
end
centroidDirection = centroidDirection ./ ...
    max(sqrt(sum(centroidDirection.^2, 3)), eps);
direction = repmat(centroidDirection, 1, M, 1);
end

function direction = forceInward(P, direction, cfg, region)
direction = direction ./ max(sqrt(sum(direction.^2, 3)), eps);
if region == "ROI"
    azimuth = atan2(P(:, :, 2), P(:, :, 1));
    alpha = deg2rad(cfg.system.roiHalfAngle_deg);
    normal = cat(3, cos(azimuth), sin(azimuth), ...
        -tan(alpha) .* ones(size(azimuth)));
    normal = normal ./ sqrt(sum(normal.^2, 3));
    outwardComponent = sum(direction .* normal, 3);
else
    normal = P ./ max(sqrt(sum(P.^2, 3)), eps);
    outwardComponent = sum(direction .* normal, 3);
end
outward = reshape(outwardComponent >= 0, size(outwardComponent, 1), ...
    size(outwardComponent, 2), 1);
outwardComponent = reshape(outwardComponent, size(outwardComponent, 1), ...
    size(outwardComponent, 2), 1);
direction = direction - 2 .* outwardComponent .* normal .* outward;
direction = direction ./ max(sqrt(sum(direction.^2, 3)), eps);
end

function inside = isInsideRegion(P, system, region)
ranges = sqrt(sum(P.^2, 3));
if region == "ROE"
    inside = P(:, :, 3) >= 0 & ranges <= system.maxRange_m + 1e-8;
else
    rho = hypot(P(:, :, 1), P(:, :, 2));
    inside = P(:, :, 3) >= 0 & ranges <= system.maxRange_m + 1e-8 & ...
        rho <= P(:, :, 3) .* tand(system.roiHalfAngle_deg) + 1e-8;
end
end

function [V, omegaState] = updateRandomTurn(V, omegaState, active, cfg, stream)
dt = cfg.execution.motionUpdate_s;
rho = exp(-dt / cfg.scene.turnCorrelation_s);
sigma = deg2rad(cfg.scene.turnRateRms_deg_s);
omegaState = rho .* omegaState + sigma .* sqrt(1 - rho^2) .* ...
    randn(stream, size(omegaState));
omegaNorm = sqrt(sum(omegaState.^2, 3));
limit = deg2rad(cfg.scene.turnRateLimit_deg_s);
scale = min(1, limit ./ max(omegaNorm, eps));
omegaState = omegaState .* scale;

speed = sqrt(sum(V.^2, 3));
unitVelocity = V ./ max(speed, eps);
turnDerivative = cross(omegaState, unitVelocity, 3);
active3 = reshape(active, size(active, 1), size(active, 2), 1);
unitVelocity = unitVelocity + dt .* turnDerivative .* active3;
unitVelocity = unitVelocity ./ max(sqrt(sum(unitVelocity.^2, 3)), eps);
V = unitVelocity .* speed;
end

function illumination = illuminationFactor(cfg, beam, U, ranges, B, scan, index)
system = cfg.system;
cosTheta = sum(U .* B, 3);
switch lower(char(beam.name))
    case 'point'
        width = beam.wD_rad / 2;
        fovHalf = 3 * width;
        targetHalf = atan((sqrt(2) * system.targetWidth_m / 2) ./ max(ranges, eps));
        candidate = cosTheta >= cos(min(fovHalf + targetHalf, pi));
        theta = safeAcos(cosTheta);
        illumination = sampledFactor('point', theta, zeros(size(theta)), ...
            ranges, system, beam, candidate);
    case 'ring'
        radius = beam.wD_rad / 2;
        width = beam.wd_rad / 2;
        theta = safeAcos(cosTheta);
        targetHalf = atan((sqrt(2) * system.targetWidth_m / 2) ./ max(ranges, eps));
        candidate = abs(theta - radius) <= targetHalf + ...
            system.gaussianCutoffWidths * width;
        illumination = sampledFactor('ring', theta, zeros(size(theta)), ...
            ranges, system, beam, candidate);
    case 'line'
        longAxis = reshape(scan.longAxes(index, :), size(B));
        shortAxis = reshape(scan.shortAxes(index, :), size(B));
        forward = sum(U .* B, 3);
        longOffset = atan2(sum(U .* longAxis, 3), forward);
        shortOffset = atan2(sum(U .* shortAxis, 3), forward);
        targetHalf = atan((system.targetWidth_m / 2) ./ max(ranges, eps));
        candidate = forward > 0 & ...
            abs(longOffset) <= beam.wD_rad / 2 + targetHalf & ...
            abs(shortOffset) <= system.gaussianCutoffWidths * beam.wd_rad + targetHalf;
        illumination = sampledFactor('line', longOffset, shortOffset, ...
            ranges, system, beam, candidate);
end
end

function factor = sampledFactor(kind, firstOffset, secondOffset, ...
    ranges, system, beam, candidate)
factor = zeros(size(ranges));
if ~any(candidate, 'all')
    return;
end
[sampleX, sampleY] = sampleGrid(system.targetSampleGrid);
rangeVector = ranges(candidate);
angularSide = 2 .* atan((system.targetWidth_m / 2) ./ rangeVector);
x = firstOffset(candidate) + angularSide .* sampleX;
y = secondOffset(candidate) + angularSide .* sampleY;
switch kind
    case 'point'
        radial = sqrt(x.^2 + y.^2);
        width = beam.wD_rad / 2;
        intensity = (radial <= 3 * width) .* exp(-2 .* (radial ./ width).^2);
    case 'ring'
        radial = sqrt(x.^2 + y.^2);
        radius = beam.wD_rad / 2;
        width = beam.wd_rad / 2;
        intensity = exp(-2 .* ((radial - radius) ./ width).^2);
    case 'line'
        intensity = (abs(x) <= beam.wD_rad / 2) .* ...
            exp(-2 .* (y ./ (beam.wd_rad / 2)).^2);
end
factor(candidate) = mean(intensity, 2);
end

function [sampleX, sampleY] = sampleGrid(gridN)
persistent cachedN cachedX cachedY
if isempty(cachedN) || cachedN ~= gridN
    nodes = ((1:gridN) - 0.5) ./ gridN - 0.5;
    [xx, yy] = meshgrid(nodes, nodes);
    cachedN = gridN;
    cachedX = xx(:).';
    cachedY = yy(:).';
end
sampleX = cachedX;
sampleY = cachedY;
end

function model = buildPhotonModel(cfg, beam)
system = cfg.system;
h = 6.62607015e-34;
c = cfg.constants.speedOfLight_m_s;
photonEnergy = h * c / system.wavelength_m;
receiverArea = pi * (system.receiverDiameter_m / 2)^2;
eta = system.opticalTransmission * system.quantumEfficiency;
transmitted = system.pulseEnergy_J / photonEnergy;

switch lower(char(beam.name))
    case 'point'
        width = beam.wD_rad / 2;
        omegaBeam = pointSolidAngle(width);
        omegaFov = 2 * pi * (1 - cos(min(3 * width, pi)));
        windowPulses = ceil(system.prf_Hz * 3 * beam.wD_rad / ...
            system.scanRate_rad_s);
    case 'line'
        width = beam.wd_rad / 2;
        omegaBeam = beam.wD_rad * width * sqrt(pi / 2);
        omegaFov = beam.wD_rad * 6 * beam.wd_rad;
        windowPulses = ceil(system.prf_Hz * beam.wD_rad / ...
            system.scanRate_rad_s);
    case 'ring'
        width = beam.wd_rad / 2;
        radius = beam.wD_rad / 2;
        omegaBeam = annularSolidAngle(radius, width);
        thetaInner = max(radius - 3 * width, 0);
        thetaOuter = min(radius + 3 * width, pi);
        omegaFov = 2 * pi * (cos(thetaInner) - cos(thetaOuter));
        windowPulses = ceil(system.prf_Hz * beam.wD_rad / ...
            system.scanRate_rad_s);
end

background = system.skyRadiance_W_m2_sr_nm * ...
    system.filterBandwidth_nm * receiverArea * omegaFov * ...
    system.rangeGate_s * eta / photonEnergy;
dark = system.darkCountRate_Hz * system.rangeGate_s;
z = (0:system.noiseLutStep_m:system.maxRange_m).';
if z(end) < system.maxRange_m
    z = [z; system.maxRange_m];
end
backscatter = calcBackscatter(transmitted, eta, receiverArea, ...
    system.backscatter_m_inv_sr_inv, system.extinction_m_inv, ...
    z, system.rangeGate_s, c);

model = struct('transmittedPhotons', transmitted, ...
    'receiverArea', receiverArea, 'systemEfficiency', eta, ...
    'targetArea', system.targetWidth_m^2, ...
    'targetReflectivity', system.targetReflectivity, ...
    'extinction', system.extinction_m_inv, 'omegaBeam', omegaBeam, ...
    'backgroundPhotons', background, 'darkPhotons', dark, ...
    'backscatterPhotons', backscatter, ...
    'noiseLutStep_m', system.noiseLutStep_m, ...
    'windowPulses', max(1, windowPulses));
end

function [signal, interceptedFraction] = calcSignal( ...
    model, ranges, illumination)
signal = zeros(size(ranges));
valid = isfinite(ranges) & ranges > 0 & illumination > 0;
beamArea = model.omegaBeam .* ranges(valid).^2;
rawFraction = zeros(size(ranges));
rawFraction(valid) = model.targetArea .* illumination(valid) ./ beamArea;
individualFraction = min(rawFraction, 1);
sceneFraction = sum(individualFraction, 2);
sceneScale = min(1, 1 ./ max(sceneFraction, realmin));
interceptedFraction = individualFraction .* sceneScale;
transmission = exp(-2 * model.extinction .* ranges(valid));
signal(valid) = model.transmittedPhotons .* ...
    interceptedFraction(valid) .* ...
    model.targetReflectivity .* ...
    (model.receiverArea ./ (pi .* ranges(valid).^2)) .* ...
    transmission .* model.systemEfficiency;
end

function noise = lookupNoise(model, ranges)
index = round(ranges ./ model.noiseLutStep_m) + 1;
index(~isfinite(index)) = 1;
index = max(1, min(numel(model.backscatterPhotons), index));
noise = reshape(model.backscatterPhotons(index), size(ranges)) + ...
    model.backgroundPhotons + model.darkPhotons;
end

function value = pointSolidAngle(width)
integrand = @(theta) exp(-2 .* (theta ./ width).^2) .* sin(theta);
value = 2 * pi * integral(integrand, 0, pi, ...
    'RelTol', 1e-10, 'AbsTol', 1e-14);
end

function value = annularSolidAngle(radius, width)
u0 = -sqrt(2) * radius / width;
radialIntegral = radius * width * sqrt(pi) / (2 * sqrt(2)) .* ...
    erfc(u0) + width^2 / 4 .* exp(-2 * (radius / width)^2);
value = 2 * pi * radialIntegral;
end

function photons = calcBackscatter(Ntx, eta, receiverArea, beta, ...
    alpha, z, gate, c)
resolution = c * gate / 2;
offsets = linspace(-resolution / 2, resolution / 2, 101);
r = z(:) + offsets;
r(r <= 0) = 1e-6;
integrand = beta .* exp(-2 .* alpha .* r) ./ r.^2;
photons = Ntx .* eta .* receiverArea .* trapz(offsets, integrand, 2);
end

function value = safeAcos(value)
value = acos(max(min(value, 1), -1));
end

function id = makeCaseId(region, trajectory, scenario, beam, targetCount)
id = sprintf('%s_%s_%s_%s_%gmrad_M%d', char(region), char(trajectory), ...
    char(scenario), char(string(beam.name)), beam.wD_rad * 1e3, ...
    targetCount);
end

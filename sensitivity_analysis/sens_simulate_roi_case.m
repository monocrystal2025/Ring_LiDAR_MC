function result = sens_simulate_roi_case(cfg, beam, trajectory, bank)
%SENS_SIMULATE_ROI_CASE Run one paired ROI Monte Carlo beam/path case.
%
% RESULT contains per-scenario detection flags and first-warning times. The
% implementation is standalone but follows the geometry, photon budget,
% 7-by-7 finite-target sampling, and variable-length sliding-window SNR
% decision used by the current ROI_NEW model.

arguments
    cfg (1,1) struct
    beam (1,1) struct
    trajectory (1,1) string
    bank (1,1) struct
end

validateInputs(cfg, beam, trajectory, bank);
trajectory = upper(trajectory);
system = cfg.system;
N = bank.N;

[positions, velocities] = mapScenarioBank(bank, system);
preparedPath = preparePath(system, beam, trajectory);
photonModel = buildPhotonModel(system, beam);
initialIndex = floor(bank.scanPhase_u .* preparedPath.cycleLength) + 1;
initialIndex = min(initialIndex, preparedPath.cycleLength);

detect = false(N, 1);
firstTime_s = nan(N, 1);
blockSize = cfg.execution.blockSize;
useParallel = cfg.execution.useParallel && ...
    license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));

timerValue = tic;
if useParallel
    pool = gcp('nocreate');
    if isempty(pool)
        try
            parpool('threads');
        catch
            parpool('local');
        end
    end
    parfor scenarioIdx = 1:N
        [detect(scenarioIdx), firstTime_s(scenarioIdx)] = simulateTarget( ...
            positions(scenarioIdx, :), velocities(scenarioIdx, :), ...
            initialIndex(scenarioIdx), system, beam, preparedPath, ...
            photonModel, blockSize);
    end
else
    for scenarioIdx = 1:N
        [detect(scenarioIdx), firstTime_s(scenarioIdx)] = simulateTarget( ...
            positions(scenarioIdx, :), velocities(scenarioIdx, :), ...
            initialIndex(scenarioIdx), system, beam, preparedPath, ...
            photonModel, blockSize);
    end
end

result = struct();
result.beam = string(beam.name);
result.beamLabel = string(beam.label);
result.trajectory = trajectory;
result.scenarioId = bank.scenarioId;
result.detect = detect;
result.firstTime_s = firstTime_s;
result.initialPhase_u = bank.scanPhase_u;
result.pathPointCount = preparedPath.pathPointCount;
result.cycleLength = preparedPath.cycleLength;
result.windowPulses = photonModel.windowPulses;
result.elapsed_s = toc(timerValue);
end

function validateInputs(cfg, beam, trajectory, bank)
requiredCfg = {'system', 'execution'};
requiredBank = {'N', 'scenarioId', 'entryAzimuth_u', 'entryRange_u', ...
    'velocityAzimuth_u', 'velocityCosPolar_u', 'scanPhase_u'};
if ~all(isfield(cfg, requiredCfg))
    error('sens_simulate_roi_case:InvalidConfig', ...
        'cfg must contain system and execution fields.');
end
if ~all(isfield(bank, requiredBank))
    error('sens_simulate_roi_case:InvalidScenarioBank', ...
        'Scenario bank is missing one or more required fields.');
end
if ~isfield(beam, 'name') || ~isfield(beam, 'wD_rad') || ~isfield(beam, 'wd_rad')
    error('sens_simulate_roi_case:InvalidBeam', ...
        'beam must contain name, wD_rad, and wd_rad.');
end
if ~any(strcmpi(string(beam.name), ["point", "line", "ring"]))
    error('sens_simulate_roi_case:UnknownBeam', ...
        'Beam name must be point, line, or ring.');
end
if ~any(strcmpi(trajectory, ["RA", "SP"]))
    error('sens_simulate_roi_case:UnknownTrajectory', ...
        'Trajectory must be RA or SP.');
end
if bank.N ~= numel(bank.scenarioId)
    error('sens_simulate_roi_case:ScenarioCountMismatch', ...
        'bank.N must match the number of scenarios.');
end
end

function [P, V] = mapScenarioBank(bank, system)
alpha = deg2rad(system.roiHalfAngle_deg);
azimuth = 2 * pi .* bank.entryAzimuth_u;
slantRange = system.minimumEntryRange_m + ...
    (system.maxRange_m - system.minimumEntryRange_m) .* bank.entryRange_u;

P = [ ...
    slantRange .* sin(alpha) .* cos(azimuth), ...
    slantRange .* sin(alpha) .* sin(azimuth), ...
    slantRange .* cos(alpha)];

velocityAzimuth = 2 * pi .* bank.velocityAzimuth_u;
velocityCosPolar = 2 .* bank.velocityCosPolar_u - 1;
velocityPolar = acos(velocityCosPolar);
direction = [ ...
    sin(velocityPolar) .* cos(velocityAzimuth), ...
    sin(velocityPolar) .* sin(velocityAzimuth), ...
    cos(velocityPolar)];

normal = [cos(azimuth), sin(azimuth), -tan(alpha) .* ones(bank.N, 1)];
normal = normal ./ vecnorm(normal, 2, 2);
outward = sum(direction .* normal, 2) >= 0;
direction(outward, :) = direction(outward, :) - ...
    2 .* sum(direction(outward, :) .* normal(outward, :), 2) .* normal(outward, :);
V = direction .* system.targetSpeed_m_s;
end

function prepared = preparePath(system, beam, trajectory)
scanRadius = tan(deg2rad(system.roiHalfAngle_deg));
beamDiameter = 2 * tan(beam.wD_rad / 2);
stepSize = system.scanRate_rad_s / system.prf_Hz;

if strcmpi(trajectory, 'RA')
    path = generateRasterPath(beamDiameter, stepSize, scanRadius);
else
    path = generateSpiralPath(scanRadius, beamDiameter, stepSize);
end

pathCycle = [path; flipud(path)];
prepared = struct();
prepared.pathPointCount = size(path, 1);
prepared.cycleLength = size(pathCycle, 1);
prepared.B = tangentToUnit(pathCycle);

if strcmpi(beam.name, 'line')
    if strcmpi(trajectory, 'SP')
        [eLong, eShort] = localFrameSpiral(path, beamDiameter);
    else
        eLong = repmat([1, 0], size(path, 1), 1);
        eShort = repmat([0, 1], size(path, 1), 1);
    end
    eLong = [eLong; flipud(eLong)];
    eShort = [eShort; flipud(eShort)];
    prepared.EL = liftTangent(prepared.B, eLong);
    prepared.ES = liftTangent(prepared.B, eShort);
else
    prepared.EL = zeros(0, 3);
    prepared.ES = zeros(0, 3);
end
end

function path = generateSpiralPath(rEnd, pitch, stepSize)
if pitch <= 0 || stepSize <= 0
    error('sens_simulate_roi_case:InvalidSpiralStep', ...
        'Spiral pitch and step size must be positive.');
end
a = pitch / (2 * pi);
thetaMax = rEnd / a;
arcLength = 0.5 * a * (thetaMax * sqrt(thetaMax^2 + 1) + ...
    log(thetaMax + sqrt(thetaMax^2 + 1)));
path = zeros(ceil(arcLength / stepSize) + 2, 2);
path(1, :) = [0, 0];
count = 1;
theta = 0;
while true
    dtheta = stepSize / (a * sqrt(theta^2 + 1));
    nextTheta = theta + dtheta;
    nextRadius = a * nextTheta;
    count = count + 1;
    if nextRadius >= rEnd
        finalTheta = rEnd / a;
        path(count, :) = [rEnd * cos(finalTheta), rEnd * sin(finalTheta)];
        break;
    end
    path(count, :) = [nextRadius * cos(nextTheta), nextRadius * sin(nextTheta)];
    theta = nextTheta;
end
path = path(1:count, :);
end

function path = generateRasterPath(D, stepSize, radius)
if D <= 0 || stepSize <= 0 || radius <= 0
    error('sens_simulate_roi_case:InvalidRasterStep', ...
        'Raster spacing, pulse step, and ROI radius must be positive.');
end
if D >= 2 * radius
    y = (-radius-stepSize:stepSize:radius+stepSize).';
    path = [zeros(size(y)), y];
    return;
end

lineCount = max(1, ceil(2 * radius / D));
parts = cell(2 * lineCount - 1, 1);
partIdx = 0;
for lineIdx = 1:lineCount
    x0 = -radius + D / 2 + D * (lineIdx - 1);
    y = linspace(-radius, radius, max(2, ceil(2 * radius / stepSize))).';
    if mod(lineIdx, 2) == 0
        y = flipud(y);
    end
    partIdx = partIdx + 1;
    parts{partIdx} = [x0 .* ones(size(y)), y];
    if lineIdx ~= lineCount
        xNext = -radius + D / 2 + D * lineIdx;
        connectorCount = max(2, ceil(abs(xNext - x0) / stepSize));
        xConnector = linspace(x0, xNext, connectorCount).';
        partIdx = partIdx + 1;
        parts{partIdx} = [xConnector(2:end), y(end) .* ones(connectorCount - 1, 1)];
    end
end
path = vertcat(parts{1:partIdx});
path = path(hypot(path(:, 1), path(:, 2)) <= radius + eps(radius), :);
if isempty(path)
    path = [0, 0];
end
end

function [eLong, eShort] = localFrameSpiral(path, pitch)
x = path(:, 1);
y = path(:, 2);
r = hypot(x, y);
b = pitch / (2 * pi);
tx = zeros(size(x));
ty = zeros(size(y));
mask = r > eps;
tx(mask) = (b ./ r(mask)) .* x(mask) - y(mask);
ty(mask) = (b ./ r(mask)) .* y(mask) + x(mask);
tx(~mask) = 1;
ty(~mask) = 0;
normTangent = hypot(tx, ty);
eShort = [tx ./ normTangent, ty ./ normTangent];
eLong = [-eShort(:, 2), eShort(:, 1)];
end

function B = tangentToUnit(path)
B = [path(:, 1), path(:, 2), ones(size(path, 1), 1)];
B = B ./ vecnorm(B, 2, 2);
end

function E = liftTangent(B, e2)
E0 = [e2(:, 1), e2(:, 2), zeros(size(e2, 1), 1)];
E = E0 - B .* sum(B .* E0, 2);
E = E ./ vecnorm(E, 2, 2);
end

function model = buildPhotonModel(system, beam)
h = 6.62607015e-34;
c = 299792458;
photonEnergy = h * c / system.wavelength_m;
receiverArea = pi * (system.receiverDiameter_m / 2)^2;
systemEfficiency = system.opticalTransmission * system.quantumEfficiency;
transmittedPhotons = system.pulseEnergy_J / photonEnergy;
darkPhotons = system.darkCountRate_Hz * system.rangeGate_s;

switch lower(char(beam.name))
    case 'point'
        gaussianWidth = beam.wD_rad / 2;
        fovHalfAngle = min(3 * gaussianWidth, pi);
        omegaBeam = pointGaussianSolidAngle(gaussianWidth);
        omegaFov = 2 * pi * (1 - cos(fovHalfAngle));
        windowPulses = max(1, ceil(system.prf_Hz * (3 * beam.wD_rad) / ...
            system.scanRate_rad_s));
    case 'line'
        gaussianWidth = beam.wd_rad / 2;
        omegaBeam = beam.wD_rad * gaussianWidth * sqrt(pi / 2);
        omegaFov = beam.wD_rad * (6 * beam.wd_rad);
        windowPulses = max(1, ceil(system.prf_Hz * beam.wD_rad / ...
            system.scanRate_rad_s));
    case 'ring'
        gaussianWidth = beam.wd_rad / 2;
        ringRadius = beam.wD_rad / 2;
        omegaBeam = annularGaussianSolidAngle(ringRadius, gaussianWidth);
        fovHalfWidth = 3 * gaussianWidth;
        thetaInner = max(ringRadius - fovHalfWidth, 0);
        thetaOuter = min(ringRadius + fovHalfWidth, pi);
        omegaFov = 2 * pi * (cos(thetaInner) - cos(thetaOuter));
        windowPulses = max(1, ceil(system.prf_Hz * beam.wD_rad / ...
            system.scanRate_rad_s));
end

backgroundPhotons = system.skyRadiance_W_m2_sr_nm * ...
    system.filterBandwidth_nm * receiverArea * omegaFov * ...
    system.rangeGate_s * systemEfficiency / photonEnergy;

z = (0:system.noiseLutStep_m:system.maxRange_m).';
if z(end) < system.maxRange_m
    z = [z; system.maxRange_m];
end
backscatterPhotons = calcBackscatter(transmittedPhotons, systemEfficiency, ...
    receiverArea, system.backscatter_m_inv_sr_inv, ...
    system.extinction_m_inv, z, system.rangeGate_s, c);

model = struct();
model.transmittedPhotons = transmittedPhotons;
model.receiverArea = receiverArea;
model.systemEfficiency = systemEfficiency;
model.targetArea = system.targetWidth_m^2;
model.targetReflectivity = system.targetReflectivity;
model.extinction = system.extinction_m_inv;
model.omegaBeam = omegaBeam;
model.backgroundPhotons = backgroundPhotons;
model.darkPhotons = darkPhotons;
model.backscatterPhotons = backscatterPhotons;
model.noiseLutStep_m = system.noiseLutStep_m;
model.windowPulses = windowPulses;
end

function [detected, firstTime_s] = simulateTarget(UP, UV, initialIndex, ...
    system, beam, path, photonModel, blockSize)
detected = false;
firstTime_s = nan;
carrySignal = zeros(0, 1);
carryNoise = zeros(0, 1);
carryStep = zeros(0, 1);

if ~isInsideCone(UP, system)
    return;
end

[B0, EL0, ES0] = pathFrame(path, initialIndex, 1);
range0 = norm(UP);
illumination0 = illuminationFactor(beam, system, UP ./ range0, range0, ...
    B0, EL0, ES0);
noise0 = lookupNoise(photonModel, range0);
signal0 = calcSignal(photonModel, range0, illumination0);
[detected, firstStep, carrySignal, carryNoise, carryStep] = checkWindow( ...
    carrySignal, carryNoise, carryStep, signal0, noise0, 0, ...
    photonModel.windowPulses, system.snrThreshold);
if detected
    firstTime_s = firstStep / system.prf_Hz;
    return;
end

basePosition = UP;
beamStart = initialIndex + 1;
stepOffset = 0;
while true
    steps = (1:blockSize).';
    time = steps ./ system.prf_Hz;
    positions = basePosition + time .* UV;
    [B, EL, ES] = pathFrame(path, beamStart, blockSize);
    inside = isInsideCone(positions, system);
    firstOutside = find(~inside, 1, 'first');
    if isempty(firstOutside)
        validCount = blockSize;
    else
        validCount = firstOutside - 1;
    end

    if validCount > 0
        validPositions = positions(1:validCount, :);
        ranges = vecnorm(validPositions, 2, 2);
        unitDirections = validPositions ./ ranges;
        illumination = illuminationFactor(beam, system, unitDirections, ranges, ...
            B(1:validCount, :), EL(1:min(validCount, size(EL, 1)), :), ...
            ES(1:min(validCount, size(ES, 1)), :));
        noise = lookupNoise(photonModel, ranges);
        signal = calcSignal(photonModel, ranges, illumination);
        stepNumbers = stepOffset + (1:validCount).';
        [detected, firstStep, carrySignal, carryNoise, carryStep] = checkWindow( ...
            carrySignal, carryNoise, carryStep, signal, noise, stepNumbers, ...
            photonModel.windowPulses, system.snrThreshold);
        if detected
            firstTime_s = firstStep / system.prf_Hz;
            return;
        end
    end

    if ~isempty(firstOutside)
        return;
    end
    basePosition = basePosition + (blockSize / system.prf_Hz) .* UV;
    beamStart = beamStart + blockSize;
    stepOffset = stepOffset + blockSize;
end
end

function inside = isInsideCone(P, system)
if isvector(P)
    P = P(:).';
end
range = vecnorm(P, 2, 2);
rho = hypot(P(:, 1), P(:, 2));
inside = P(:, 3) >= 0 & range <= system.maxRange_m & ...
    rho <= P(:, 3) .* tand(system.roiHalfAngle_deg) + 1e-10;
end

function [B, EL, ES] = pathFrame(path, startIndex, count)
idx = mod((startIndex - 1) + (0:count-1), path.cycleLength) + 1;
B = path.B(idx, :);
if isempty(path.EL)
    EL = zeros(count, 0);
    ES = zeros(count, 0);
else
    EL = path.EL(idx, :);
    ES = path.ES(idx, :);
end
end

function illumination = illuminationFactor(beam, system, U, ranges, B, EL, ES)
ranges = ranges(:);
if isvector(U)
    U = U(:).';
end
switch lower(char(beam.name))
    case 'point'
        gaussianWidth = beam.wD_rad / 2;
        fovHalfAngle = 3 * gaussianWidth;
        cosTheta = sum(U .* B, 2);
        targetHalfDiagonal = atan((sqrt(2) * system.targetWidth_m / 2) ./ ranges);
        candidate = cosTheta >= cos(min(fovHalfAngle + targetHalfDiagonal, pi));
        theta = safeAcos(cosTheta(candidate));
        illumination = zeros(size(ranges));
        illumination(candidate) = pointTargetFactor(theta, ranges(candidate), ...
            system.targetWidth_m, gaussianWidth, fovHalfAngle, ...
            system.targetSampleGrid);
    case 'ring'
        ringRadius = beam.wD_rad / 2;
        gaussianWidth = beam.wd_rad / 2;
        theta = safeAcos(sum(U .* B, 2));
        targetHalfDiagonal = atan((sqrt(2) * system.targetWidth_m / 2) ./ ranges);
        candidate = abs(theta - ringRadius) <= targetHalfDiagonal + ...
            system.gaussianCutoffWidths * gaussianWidth;
        illumination = zeros(size(ranges));
        illumination(candidate) = ringTargetFactor(theta(candidate), ...
            ranges(candidate), system.targetWidth_m, ringRadius, ...
            gaussianWidth, system.targetSampleGrid);
    case 'line'
        c3 = sum(U .* B, 2);
        longOffset = atan2(sum(U .* EL, 2), c3);
        shortOffset = atan2(sum(U .* ES, 2), c3);
        targetHalf = atan((system.targetWidth_m / 2) ./ ranges);
        longHalfWidth = beam.wD_rad / 2;
        gaussianWidth = beam.wd_rad / 2;
        shortPrefilter = system.gaussianCutoffWidths * beam.wd_rad;
        candidate = c3 > 0 & abs(longOffset) <= longHalfWidth + targetHalf & ...
            abs(shortOffset) <= shortPrefilter + targetHalf;
        illumination = zeros(size(ranges));
        illumination(candidate) = lineTargetFactor(longOffset(candidate), ...
            shortOffset(candidate), ranges(candidate), system.targetWidth_m, ...
            longHalfWidth, gaussianWidth, system.targetSampleGrid);
end
end

function factor = pointTargetFactor(theta, z, targetWidth, gaussianWidth, fovHalf, gridN)
[sampleX, sampleY] = sampleGrid(gridN);
angularSide = 2 .* atan((targetWidth / 2) ./ z(:));
sampleTheta = sqrt((theta(:) + angularSide .* sampleX).^2 + ...
    (angularSide .* sampleY).^2);
intensity = (sampleTheta <= fovHalf) .* exp(-2 .* (sampleTheta ./ gaussianWidth).^2);
factor = mean(intensity, 2);
end

function factor = ringTargetFactor(theta, z, targetWidth, ringRadius, gaussianWidth, gridN)
[sampleX, sampleY] = sampleGrid(gridN);
angularSide = 2 .* atan((targetWidth / 2) ./ z(:));
sampleTheta = sqrt((theta(:) + angularSide .* sampleX).^2 + ...
    (angularSide .* sampleY).^2);
intensity = exp(-2 .* ((sampleTheta - ringRadius) ./ gaussianWidth).^2);
factor = mean(intensity, 2);
end

function factor = lineTargetFactor(longOffset, shortOffset, z, targetWidth, ...
    longHalfWidth, gaussianWidth, gridN)
[sampleLong, sampleShort] = sampleGrid(gridN);
angularSide = 2 .* atan((targetWidth / 2) ./ z(:));
sampleL = longOffset(:) + angularSide .* sampleLong;
sampleS = shortOffset(:) + angularSide .* sampleShort;
intensity = (abs(sampleL) <= longHalfWidth) .* ...
    exp(-2 .* (sampleS ./ gaussianWidth).^2);
factor = mean(intensity, 2);
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

function noise = lookupNoise(model, ranges)
idx = round(ranges(:) ./ model.noiseLutStep_m) + 1;
idx(~isfinite(idx)) = 1;
idx = max(1, min(numel(model.backscatterPhotons), idx));
noise = model.backscatterPhotons(idx) + model.backgroundPhotons + model.darkPhotons;
end

function signal = calcSignal(model, ranges, illumination)
ranges = ranges(:);
illumination = illumination(:);
signal = zeros(size(ranges));
valid = isfinite(ranges) & ranges > 0 & isfinite(illumination) & ...
    illumination > 0 & model.omegaBeam > 0;
beamArea = model.omegaBeam .* ranges(valid).^2;
twoWayTransmission = exp(-2 * model.extinction .* ranges(valid));
signal(valid) = model.transmittedPhotons .* ...
    (model.targetArea .* illumination(valid) ./ beamArea) .* ...
    model.targetReflectivity .* (model.receiverArea ./ (pi .* ranges(valid).^2)) .* ...
    twoWayTransmission .* model.systemEfficiency;
end

function [detected, firstStep, carrySignal, carryNoise, carryStep] = checkWindow( ...
    carrySignal, carryNoise, carryStep, signal, noise, steps, windowPulses, threshold)
detected = false;
firstStep = nan;
signal = signal(:);
noise = noise(:);
steps = steps(:);
signal(~isfinite(signal)) = 0;
noise(~isfinite(noise)) = 0;

carryCount = numel(carrySignal);
allSignal = [carrySignal; signal];
allNoise = [carryNoise; noise];
allSteps = [carryStep; steps];
signalPositions = carryCount + find(signal > 0);

if ~isempty(signalPositions)
    tail = signalPositions(:);
    lengths = 1:min(windowPulses, numel(allSignal));
    starts = tail - lengths + 1;
    validLength = starts >= 1;
    starts(~validLength) = 1;
    variance = allSignal + allNoise;
    score = zeros(size(allSignal));
    validScore = allSignal > 0 & variance > 0 & isfinite(variance);
    score(validScore) = allSignal(validScore).^2 ./ variance(validScore);
    cumulativeScore = [0; cumsum(score)];
    windowScore = cumulativeScore(tail + 1) - cumulativeScore(starts);
    windowSnr = sqrt(windowScore);
    windowSnr(~validLength | windowScore <= 0 | ~isfinite(windowSnr)) = -inf;
    hit = find(any(windowSnr >= threshold, 2), 1, 'first');
    if ~isempty(hit)
        detected = true;
        firstStep = allSteps(tail(hit));
        return;
    end
end

keepCount = min(windowPulses - 1, numel(allSignal));
if keepCount > 0
    carrySignal = allSignal(end-keepCount+1:end);
    carryNoise = allNoise(end-keepCount+1:end);
    carryStep = allSteps(end-keepCount+1:end);
else
    carrySignal = zeros(0, 1);
    carryNoise = zeros(0, 1);
    carryStep = zeros(0, 1);
end
end

function omega = pointGaussianSolidAngle(width)
integrand = @(theta) exp(-2 .* (theta ./ width).^2) .* sin(theta);
omega = 2 * pi * integral(integrand, 0, pi, 'RelTol', 1e-10, 'AbsTol', 1e-14);
end

function omega = annularGaussianSolidAngle(radius, width)
u0 = -sqrt(2) * radius / width;
radialIntegral = radius * width * sqrt(pi) / (2 * sqrt(2)) .* erfc(u0) + ...
    width^2 / 4 .* exp(-2 * (radius / width)^2);
omega = 2 * pi * radialIntegral;
end

function photons = calcBackscatter(Ntx, eta, receiverArea, beta, alpha, z, gate, c)
rangeResolution = c * gate / 2;
offsets = linspace(-rangeResolution / 2, rangeResolution / 2, 101);
r = z(:) + offsets;
r(r <= 0) = 1e-6;
integrand = beta .* exp(-2 .* alpha .* r) ./ r.^2;
photons = Ntx .* eta .* receiverArea .* trapz(offsets, integrand, 2);
end

function value = safeAcos(value)
value = acos(max(min(value, 1), -1));
end

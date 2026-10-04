function result = ea_sens_simulate_case(cfg, beam, bank)
%EA_SENS_SIMULATE_CASE Run one standalone ROE Monte Carlo beam case.
%
% This implementation follows the physical process used by MC_EA:
% uniformly distributed entries on the upper hemispherical boundary,
% inward constant-velocity targets, an equidistant spherical spiral with
% boundary circles, finite-target illumination, photon-counting noise, and
% optimal multi-pulse SNR accumulation. It does not call or modify MC_EA.m.

arguments
    cfg (1,1) struct
    beam (1,1) struct
    bank (1,1) struct
end

validateInputs(cfg, beam, bank);
system = cfg.system;
N = bank.N;

[positions, velocities] = mapScenarioBank(bank, system);
preparedPath = prepareRoePath(system, beam.wD_rad);
photonModel = buildPhotonModel(system, beam);
initialIndex = floor(bank.scanPhase_u .* preparedPath.cycleLength) + 1;
initialIndex = min(initialIndex, preparedPath.cycleLength);

detect = false(N, 1);
firstTime_s = nan(N, 1);
firstEncounter_s = nan(N, 1);
effectivePulses = nan(N, 1);
blockSize = cfg.execution.blockSize;
useParallel = cfg.execution.useParallel && ...
    license("test", "Distrib_Computing_Toolbox") && ~isempty(ver("parallel"));

timerValue = tic;
if useParallel
    pool = gcp("nocreate");
    if isempty(pool)
        try
            parpool("threads");
        catch
            parpool("local");
        end
    end
    parfor scenarioIdx = 1:N
        [detect(scenarioIdx), firstTime_s(scenarioIdx), ...
            firstEncounter_s(scenarioIdx), effectivePulses(scenarioIdx)] = ...
            simulateTarget(positions(scenarioIdx, :), ...
            velocities(scenarioIdx, :), initialIndex(scenarioIdx), ...
            system, beam, preparedPath, photonModel, blockSize);
    end
else
    for scenarioIdx = 1:N
        [detect(scenarioIdx), firstTime_s(scenarioIdx), ...
            firstEncounter_s(scenarioIdx), effectivePulses(scenarioIdx)] = ...
            simulateTarget(positions(scenarioIdx, :), ...
            velocities(scenarioIdx, :), initialIndex(scenarioIdx), ...
            system, beam, preparedPath, photonModel, blockSize);
    end
end

result = struct();
result.beam = string(beam.name);
result.beamLabel = string(beam.label);
result.wD_rad = beam.wD_rad;
result.wd_rad = beam.wd_rad;
result.scenarioId = bank.scenarioId;
result.detect = detect;
result.firstTime_s = firstTime_s;
result.firstEncounter_s = firstEncounter_s;
result.effectivePulses = effectivePulses;
result.initialPhase_u = bank.scanPhase_u;
result.pathPointCount = preparedPath.pathPointCount;
result.cycleLength = preparedPath.cycleLength;
result.windowPulses = photonModel.windowPulses;
result.elapsed_s = toc(timerValue);
end

function validateInputs(cfg, beam, bank)
hasConfig = isfield(cfg, 'system') && isfield(cfg, 'execution');
hasBank = isfield(bank, 'N') && isfield(bank, 'scenarioId') && ...
    isfield(bank, 'entryAzimuth_u') && ...
    isfield(bank, 'entryCosColatitude_u') && ...
    isfield(bank, 'velocityAzimuth_u') && ...
    isfield(bank, 'velocityCosPolar_u') && ...
    isfield(bank, 'scanPhase_u');
hasBeam = isfield(beam, 'name') && isfield(beam, 'label') && ...
    isfield(beam, 'wD_rad') && isfield(beam, 'wd_rad');
if ~hasConfig
    error("ea_sens_simulate_case:InvalidConfig", ...
        "cfg must contain system and execution fields.");
end
if ~hasBank
    error("ea_sens_simulate_case:InvalidScenarioBank", ...
        "Scenario bank is missing required normalized draws.");
end
if ~hasBeam
    error("ea_sens_simulate_case:InvalidBeam", ...
        "beam must contain name, label, wD_rad, and wd_rad.");
end
if ~any(strcmpi(string(beam.name), ["point", "line", "ring"]))
    error("ea_sens_simulate_case:UnknownBeam", ...
        "Beam name must be point, line, or ring.");
end
if bank.N ~= numel(bank.scenarioId)
    error("ea_sens_simulate_case:ScenarioCountMismatch", ...
        "bank.N must match the number of scenario IDs.");
end
if beam.wD_rad <= 0 || (~strcmpi(beam.name, "point") && beam.wd_rad <= 0)
    error("ea_sens_simulate_case:InvalidDivergence", ...
        "Beam divergence angles must be positive.");
end
end

function [P, V] = mapScenarioBank(bank, system)
azimuth = 2 * pi .* bank.entryAzimuth_u;
colatitude = acos(bank.entryCosColatitude_u);
range = system.maxRange_m;
P = range .* [ ...
    sin(colatitude) .* cos(azimuth), ...
    sin(colatitude) .* sin(azimuth), ...
    cos(colatitude)];

velocityAzimuth = 2 * pi .* bank.velocityAzimuth_u;
velocityCosPolar = 2 .* bank.velocityCosPolar_u - 1;
velocityPolar = acos(velocityCosPolar);
direction = [ ...
    sin(velocityPolar) .* cos(velocityAzimuth), ...
    sin(velocityPolar) .* sin(velocityAzimuth), ...
    cos(velocityPolar)];

outward = sum(P .* direction, 2) > 0;
direction(outward, :) = -direction(outward, :);
V = system.targetSpeed_m_s .* direction;
end

function prepared = prepareRoePath(system, wD_rad)
% Generate one continuous ROE cycle at constant angular scan speed.

persistent cachedKey cachedPath
key = sprintf("%.17g|%.17g|%.17g", ...
    wD_rad, system.prf_Hz, system.scanRate_rad_s);
if ~isempty(cachedKey) && strcmp(cachedKey, key)
    prepared = cachedPath;
    return;
end

angularStep = system.scanRate_rad_s / system.prf_Hz;
halfWidth = wD_rad / 2;
phiMin = min(max(halfWidth, eps), pi / 2);
phiMax = max(phiMin, pi / 2 - halfWidth);
pitch = wD_rad;

spiral = generateSphericalSpiral(phiMin, phiMax, pitch, angularStep);
startAzimuth = atan2(spiral(1, 2), spiral(1, 1));
endAzimuth = atan2(spiral(end, 2), spiral(end, 1));
lowerCircle = generateLatitudeCircle(phiMin, angularStep, startAzimuth);
upperCircle = generateLatitudeCircle(phiMax, angularStep, endAzimuth);

reverseSpiral = flipud(spiral);
if size(reverseSpiral, 1) > 2
    reverseInterior = reverseSpiral(2:end-1, :);
else
    reverseInterior = zeros(0, 3);
end
cycle = [ ...
    lowerCircle; ...
    spiral(2:end, :); ...
    upperCircle(2:end, :); ...
    reverseInterior];
cycle = cycle ./ vecnorm(cycle, 2, 2);

prepared = struct();
prepared.B = cycle;
prepared.pathPointCount = size(lowerCircle, 1) + size(spiral, 1) + ...
    size(upperCircle, 1);
prepared.cycleLength = size(cycle, 1);

cachedKey = key;
cachedPath = prepared;
end

function B = generateLatitudeCircle(phi, angularStep, azimuth0)
circumference = 2 * pi * max(sin(phi), eps);
pointCount = max(4, ceil(circumference / angularStep));
azimuth = azimuth0 + (0:pointCount-1).' .* (2 * pi / pointCount);
B = [ ...
    sin(phi) .* cos(azimuth), ...
    sin(phi) .* sin(azimuth), ...
    cos(phi) .* ones(pointCount, 1)];
end

function B = generateSphericalSpiral(phiMin, phiMax, pitch, angularStep)
if phiMax <= phiMin + eps
    B = [sin(phiMin), 0, cos(phiMin)];
    return;
end

b = pitch / (2 * pi);
thetaMax = (phiMax - phiMin) / b;
speed = @(theta) sqrt(b.^2 + sin(phiMin + b .* theta).^2);
arcLength = integral(speed, 0, thetaMax, ...
    "RelTol", 1e-8, "AbsTol", 1e-10);
capacity = max(2, ceil(arcLength / angularStep) + 3);
theta = zeros(capacity, 1);
count = 1;

while theta(count) < thetaMax
    phi = phiMin + b * theta(count);
    dTheta = angularStep / sqrt(b^2 + sin(phi)^2);
    midTheta = min(theta(count) + dTheta / 2, thetaMax);
    midPhi = phiMin + b * midTheta;
    dTheta = angularStep / sqrt(b^2 + sin(midPhi)^2);
    nextTheta = min(theta(count) + dTheta, thetaMax);
    count = count + 1;
    if count > numel(theta)
        theta = [theta; zeros(max(100, ceil(0.1 * capacity)), 1)]; %#ok<AGROW>
        capacity = numel(theta);
    end
    theta(count) = nextTheta;
end

theta = theta(1:count);
phi = phiMin + b .* theta;
B = [ ...
    sin(phi) .* cos(theta), ...
    sin(phi) .* sin(theta), ...
    cos(phi)];
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
    case "point"
        gaussianWidth = beam.wD_rad / 2;
        fovHalfAngle = min(3 * gaussianWidth, pi);
        omegaBeam = pointGaussianSolidAngle(gaussianWidth);
        omegaFov = 2 * pi * (1 - cos(fovHalfAngle));
        windowPulses = max(1, ceil(system.prf_Hz * ...
            (3 * beam.wD_rad) / system.scanRate_rad_s));
    case "line"
        gaussianWidth = beam.wd_rad / 2;
        omegaBeam = beam.wD_rad * gaussianWidth * sqrt(pi / 2);
        omegaFov = beam.wD_rad * (6 * beam.wd_rad);
        windowPulses = max(1, ceil(system.prf_Hz * beam.wD_rad / ...
            system.scanRate_rad_s));
    case "ring"
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

function [detected, firstTime_s, firstEncounter_s, effectivePulses] = ...
        simulateTarget(UP, UV, initialIndex, system, beam, path, ...
        photonModel, blockSize)
detected = false;
firstTime_s = nan;
firstEncounter_s = nan;
effectivePulses = nan;
carrySignal = zeros(0, 1);
carryNoise = zeros(0, 1);
carryStep = zeros(0, 1);

if ~isInsideHemisphere(UP, system.maxRange_m)
    return;
end

[B0, EL0, ES0] = pathFrame(path, initialIndex, 1);
range0 = norm(UP);
illumination0 = illuminationFactor(beam, system, UP ./ range0, ...
    range0, B0, EL0, ES0);
if illumination0 > 0
    firstEncounter_s = 0;
end
noise0 = lookupNoise(photonModel, range0);
signal0 = calcSignal(photonModel, range0, illumination0);
[detected, firstStep, effectivePulses, carrySignal, carryNoise, carryStep] = ...
    checkWindow(carrySignal, carryNoise, carryStep, signal0, noise0, 0, ...
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
    inside = isInsideHemisphere(positions, system.maxRange_m);
    firstOutside = find(~inside, 1, "first");
    if isempty(firstOutside)
        validCount = blockSize;
    else
        validCount = firstOutside - 1;
    end

    if validCount > 0
        validPositions = positions(1:validCount, :);
        ranges = vecnorm(validPositions, 2, 2);
        unitDirections = validPositions ./ ranges;
        illumination = illuminationFactor(beam, system, unitDirections, ...
            ranges, B(1:validCount, :), EL(1:validCount, :), ...
            ES(1:validCount, :));
        stepNumbers = stepOffset + (1:validCount).';
        if isnan(firstEncounter_s)
            firstIlluminated = find(illumination > 0, 1, "first");
            if ~isempty(firstIlluminated)
                firstEncounter_s = stepNumbers(firstIlluminated) / ...
                    system.prf_Hz;
            end
        end
        noise = lookupNoise(photonModel, ranges);
        signal = calcSignal(photonModel, ranges, illumination);
        [detected, firstStep, effectivePulses, carrySignal, ...
            carryNoise, carryStep] = checkWindow( ...
            carrySignal, carryNoise, carryStep, signal, noise, ...
            stepNumbers, photonModel.windowPulses, system.snrThreshold);
        if detected
            firstTime_s = firstStep / system.prf_Hz;
            return;
        end
    end

    if ~isempty(firstOutside)
        return;
    end
    elapsedBlockTime = blockSize / system.prf_Hz;
    basePosition = basePosition + elapsedBlockTime .* UV;
    beamStart = beamStart + blockSize;
    stepOffset = stepOffset + blockSize;
end
end

function inside = isInsideHemisphere(P, maxRange)
if isvector(P)
    P = P(:).';
end
range = vecnorm(P, 2, 2);
tolerance = max(1e-9, 1e-12 * maxRange);
inside = P(:, 3) >= -tolerance & range <= maxRange + tolerance;
end

function [B, EL, ES] = pathFrame(path, startIndex, count)
idx = mod((startIndex - 1) + (0:count-1), path.cycleLength) + 1;
B = path.B(idx, :);

horizontalNorm = hypot(B(:, 1), B(:, 2));
safeNorm = horizontalNorm;
pole = horizontalNorm < eps;
safeNorm(pole) = 1;
ES = [-B(:, 2) ./ safeNorm, B(:, 1) ./ safeNorm, ...
    zeros(count, 1)];
ES(pole, :) = repmat([1, 0, 0], nnz(pole), 1);
EL = cross(B, ES, 2);
EL = EL ./ vecnorm(EL, 2, 2);
end

function illumination = illuminationFactor(beam, system, U, ranges, B, EL, ES)
ranges = ranges(:);
if isvector(U)
    U = U(:).';
end
switch lower(char(beam.name))
    case "point"
        gaussianWidth = beam.wD_rad / 2;
        fovHalfAngle = 3 * gaussianWidth;
        cosTheta = sum(U .* B, 2);
        targetHalfDiagonal = atan((sqrt(2) * system.targetWidth_m / 2) ./ ...
            ranges);
        candidate = cosTheta >= cos(min(fovHalfAngle + ...
            targetHalfDiagonal, pi));
        theta = safeAcos(cosTheta(candidate));
        illumination = zeros(size(ranges));
        illumination(candidate) = pointTargetFactor(theta, ...
            ranges(candidate), system.targetWidth_m, gaussianWidth, ...
            fovHalfAngle, system.targetSampleGrid);
    case "ring"
        ringRadius = beam.wD_rad / 2;
        gaussianWidth = beam.wd_rad / 2;
        theta = safeAcos(sum(U .* B, 2));
        targetHalfDiagonal = atan((sqrt(2) * system.targetWidth_m / 2) ./ ...
            ranges);
        candidate = abs(theta - ringRadius) <= targetHalfDiagonal + ...
            system.gaussianCutoffWidths * gaussianWidth;
        illumination = zeros(size(ranges));
        illumination(candidate) = ringTargetFactor(theta(candidate), ...
            ranges(candidate), system.targetWidth_m, ringRadius, ...
            gaussianWidth, system.targetSampleGrid);
    case "line"
        c3 = sum(U .* B, 2);
        longOffset = atan2(sum(U .* EL, 2), c3);
        shortOffset = atan2(sum(U .* ES, 2), c3);
        targetHalf = atan((system.targetWidth_m / 2) ./ ranges);
        longHalfWidth = beam.wD_rad / 2;
        gaussianWidth = beam.wd_rad / 2;
        shortPrefilter = system.gaussianCutoffWidths * beam.wd_rad;
        candidate = c3 > 0 & ...
            abs(longOffset) <= longHalfWidth + targetHalf & ...
            abs(shortOffset) <= shortPrefilter + targetHalf;
        illumination = zeros(size(ranges));
        illumination(candidate) = lineTargetFactor( ...
            longOffset(candidate), shortOffset(candidate), ...
            ranges(candidate), system.targetWidth_m, longHalfWidth, ...
            gaussianWidth, system.targetSampleGrid);
end
end

function factor = pointTargetFactor(theta, z, targetWidth, ...
        gaussianWidth, fovHalf, gridN)
[sampleX, sampleY] = sampleGrid(gridN);
angularSide = 2 .* atan((targetWidth / 2) ./ z(:));
sampleTheta = sqrt((theta(:) + angularSide .* sampleX).^2 + ...
    (angularSide .* sampleY).^2);
intensity = (sampleTheta <= fovHalf) .* ...
    exp(-2 .* (sampleTheta ./ gaussianWidth).^2);
factor = mean(intensity, 2);
end

function factor = ringTargetFactor(theta, z, targetWidth, ...
        ringRadius, gaussianWidth, gridN)
[sampleX, sampleY] = sampleGrid(gridN);
angularSide = 2 .* atan((targetWidth / 2) ./ z(:));
sampleTheta = sqrt((theta(:) + angularSide .* sampleX).^2 + ...
    (angularSide .* sampleY).^2);
intensity = exp(-2 .* ((sampleTheta - ringRadius) ./ gaussianWidth).^2);
factor = mean(intensity, 2);
end

function factor = lineTargetFactor(longOffset, shortOffset, z, ...
        targetWidth, longHalfWidth, gaussianWidth, gridN)
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
noise = model.backscatterPhotons(idx) + ...
    model.backgroundPhotons + model.darkPhotons;
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
    model.targetReflectivity .* ...
    (model.receiverArea ./ (pi .* ranges(valid).^2)) .* ...
    twoWayTransmission .* model.systemEfficiency;
end

function [detected, firstStep, effectivePulses, carrySignal, ...
        carryNoise, carryStep] = checkWindow(carrySignal, carryNoise, ...
        carryStep, signal, noise, steps, windowPulses, threshold)
detected = false;
firstStep = nan;
effectivePulses = nan;
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
    tail = reshape(signalPositions, [], 1);
    lengths = reshape(1:min(windowPulses, numel(allSignal)), 1, []);
    starts = tail + 1 - lengths;
    validLength = starts >= 1;
    starts(~validLength) = 1;
    variance = allSignal + allNoise;
    score = zeros(size(allSignal));
    validScore = allSignal > 0 & variance > 0 & isfinite(variance);
    score(validScore) = allSignal(validScore).^2 ./ variance(validScore);
    cumulativeScore = [0; cumsum(score)];
    tailScore = cumulativeScore(tail + 1);
    windowScore = repmat(tailScore, 1, numel(lengths)) - ...
        reshape(cumulativeScore(starts), size(starts));
    windowSnr = sqrt(windowScore);
    windowSnr(~validLength | windowScore <= 0 | ~isfinite(windowSnr)) = -inf;
    hitRow = find(any(windowSnr >= threshold, 2), 1, "first");
    if ~isempty(hitRow)
        successfulLength = find( ...
            windowSnr(hitRow, :) >= threshold, 1, "first");
        successfulStart = starts(hitRow, successfulLength);
        detected = true;
        firstStep = allSteps(tail(hitRow));
        effectivePulses = nnz( ...
            score(successfulStart:tail(hitRow)) > 0);
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
omega = 2 * pi * integral(integrand, 0, pi, ...
    "RelTol", 1e-10, "AbsTol", 1e-14);
end

function omega = annularGaussianSolidAngle(radius, width)
u0 = -sqrt(2) * radius / width;
radialIntegral = radius * width * sqrt(pi) / (2 * sqrt(2)) .* ...
    erfc(u0) + width^2 / 4 .* exp(-2 * (radius / width)^2);
omega = 2 * pi * radialIntegral;
end

function photons = calcBackscatter(Ntx, eta, receiverArea, beta, ...
        alpha, z, gate, c)
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

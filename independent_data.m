clc;
clearvars;

% Independent trajectory/beam-position study for a stationary beam.
% The only Monte Carlo variable is the target's initial along-track phase.
% It is uniform over one inter-pulse travel interval, so every realization
% traverses the same complete path with a different pulse-sampling lattice.

scriptDir = fileparts(mfilename("fullpath"));
if strlength(scriptDir) == 0
    scriptDir = pwd;
end

rngSeed = 20260904;
rng(rngSeed, "twister");

beamNames = ["point", "line", "ring"];
beamLabels = ["Spot", "Line", "Annular"];
numberOfTrials = 200;
distanceFractions = (0:1/5:1).';
trajectoryAngles_deg = (0:15:90).';

% MC_EA.m currently uses the fast scan (2*pi rad/s), f = 5 kHz,
% v = 30 m/s, w_D as a full angular diameter/length, and w_d = 1 mrad.
% The scan rate is retained only to define the same fixed-window lengths;
% the beam direction itself remains [0 0 1] throughout this experiment.
system = struct();
system.range_m = 2000;
system.prf_Hz = 5e3;
system.referenceScanRate_rad_s = 2 * pi;
system.targetSpeed_m_s = (system.referenceScanRate_rad_s + 30 / system.range_m) * system.range_m;
system.wD_rad = 100e-3;
system.wd_rad = 1e-3;
system.windowMode = "fixed";
system.snrThreshold = 2;
system.pulseEnergy_J = 150e-6;
system.wavelength_m = 1550e-9;
system.receiverDiameter_m = 0.0508;
system.opticalTransmission = 0.80;
system.quantumEfficiency = 0.80;
system.rangeGate_s = 50e-9;
system.extinction_m_inv = 1.5e-5;
system.backscatter_m_inv_sr_inv = 0.3e-6;
system.skyRadiance_W_m2_sr_nm = 5e-9;
system.filterBandwidth_nm = 10;
system.darkCountRate_Hz = 400;
system.targetWidth_m = 0.297;
system.targetReflectivity = 0.5;
system.gaussianCutoffWidths = 3;
system.targetSampleGrid = 7;

% Under the repository convention, w_D is a full outer diameter.  Thus the
% requested tangent trajectory is one geometric radius from the beam axis.
system.geometricRadius_rad = system.wD_rad / 2;
system.geometricRadius_m = system.range_m * tan(system.geometricRadius_rad);
system.interPulseTravel_m = system.targetSpeed_m_s / system.prf_Hz;
system.trackHalfLength_m = system.range_m * tan(1.75 * system.wD_rad);
system.initialAlongInterval_m = [ ...
    -system.trackHalfLength_m - system.interPulseTravel_m, ...
    -system.trackHalfLength_m];
system.windowPulses = [ ...
    ceil(system.prf_Hz * 3 * system.wD_rad / ...
        system.referenceScanRate_rad_s), ... % point
    ceil(system.prf_Hz * system.wD_rad / ...
        system.referenceScanRate_rad_s), ... % line
    ceil(system.prf_Hz * system.wD_rad / ...
        system.referenceScanRate_rad_s)];     % ring

photonModels = cell(1, numel(beamNames));
for beamIndex = 1:numel(beamNames)
    photonModels{beamIndex} = buildPhotonModel(system, beamNames(beamIndex));
end

useParallel = license("test", "Distrib_Computing_Toolbox") && ...
    ~isempty(ver("parallel"));
if useParallel && isempty(gcp("nocreate"))
    try
        parpool("threads");
    catch
        parpool("local");
    end
end

fprintf("Independent study: N=%d, range=%g m, speed=%g m/s.\n", ...
    numberOfTrials, system.range_m, system.targetSpeed_m_s);

% Study 1: trajectory displacement.  For circular beams, the value is the
% closest distance from the trajectory to the optical axis.  For the line
% beam, it is the distance from the beam center to the trajectory/long-axis
% intersection; the long side and trajectory are fixed at 30 degrees.  The
% intersection/closest point is only an anchor that defines the infinite
% trajectory line; it is not an initialized or pulse-sampled target point.
distanceInitialAlong_m = sampleInitialAlongPositions( ...
    system.initialAlongInterval_m, numberOfTrials, numel(distanceFractions));
distanceMetric = nan(numberOfTrials, numel(distanceFractions), ...
    numel(beamNames));
distanceDetected = false(size(distanceMetric));
distanceFirstDetection_s = nan(size(distanceMetric));
distanceInitialPosition_m = nan(numberOfTrials, 2, ...
    numel(distanceFractions), numel(beamNames));

for parameterIndex = 1:numel(distanceFractions)
    fraction = distanceFractions(parameterIndex);
    for beamIndex = 1:numel(beamNames)
        beamName = beamNames(beamIndex);
        photonModel = photonModels{beamIndex};
        initialAlong = distanceInitialAlong_m(:, parameterIndex);
        metrics = nan(numberOfTrials, 1);
        detections = false(numberOfTrials, 1);
        firstTimes = nan(numberOfTrials, 1);
        initialPositions = nan(numberOfTrials, 2);
        if useParallel
            parfor trialIndex = 1:numberOfTrials
                [metrics(trialIndex), detections(trialIndex), ...
                    firstTimes(trialIndex), initialPositions(trialIndex, :)] = ...
                    simulateTrial(system, beamName, photonModel, ...
                    initialAlong(trialIndex), ...
                    "distance", fraction);
            end
        else
            for trialIndex = 1:numberOfTrials
                [metrics(trialIndex), detections(trialIndex), ...
                    firstTimes(trialIndex), initialPositions(trialIndex, :)] = ...
                    simulateTrial(system, beamName, photonModel, ...
                    initialAlong(trialIndex), ...
                    "distance", fraction);
            end
        end
        distanceMetric(:, parameterIndex, beamIndex) = metrics;
        distanceDetected(:, parameterIndex, beamIndex) = detections;
        distanceFirstDetection_s(:, parameterIndex, beamIndex) = firstTimes;
        distanceInitialPosition_m(:, :, parameterIndex, beamIndex) = ...
            initialPositions;
    end
    fprintf("  Distance case %d/%d complete.\n", ...
        parameterIndex, numel(distanceFractions));
end

% Study 2: trajectory direction.  Circular-beam paths keep a closest-axis
% distance of 2/5 geometric radius; their angle is measured from horizontal.
% The line-beam path crosses its center, and the parameter is the angle from
% the line long side (the horizontal x axis).  Again, the center/closest
% point anchors the trajectory line and is not forced onto the pulse grid.
angleInitialAlong_m = sampleInitialAlongPositions( ...
    system.initialAlongInterval_m, numberOfTrials, numel(trajectoryAngles_deg));
angleMetric = nan(numberOfTrials, numel(trajectoryAngles_deg), ...
    numel(beamNames));
angleDetected = false(size(angleMetric));
angleFirstDetection_s = nan(size(angleMetric));
angleInitialPosition_m = nan(numberOfTrials, 2, ...
    numel(trajectoryAngles_deg), numel(beamNames));

for parameterIndex = 1:numel(trajectoryAngles_deg)
    angle_deg = trajectoryAngles_deg(parameterIndex);
    for beamIndex = 1:numel(beamNames)
        beamName = beamNames(beamIndex);
        photonModel = photonModels{beamIndex};
        initialAlong = angleInitialAlong_m(:, parameterIndex);
        metrics = nan(numberOfTrials, 1);
        detections = false(numberOfTrials, 1);
        firstTimes = nan(numberOfTrials, 1);
        initialPositions = nan(numberOfTrials, 2);
        if useParallel
            parfor trialIndex = 1:numberOfTrials
                [metrics(trialIndex), detections(trialIndex), ...
                    firstTimes(trialIndex), initialPositions(trialIndex, :)] = ...
                    simulateTrial(system, beamName, photonModel, ...
                    initialAlong(trialIndex), ...
                    "angle", angle_deg);
            end
        else
            for trialIndex = 1:numberOfTrials
                [metrics(trialIndex), detections(trialIndex), ...
                    firstTimes(trialIndex), initialPositions(trialIndex, :)] = ...
                    simulateTrial(system, beamName, photonModel, ...
                    initialAlong(trialIndex), ...
                    "angle", angle_deg);
            end
        end
        angleMetric(:, parameterIndex, beamIndex) = metrics;
        angleDetected(:, parameterIndex, beamIndex) = detections;
        angleFirstDetection_s(:, parameterIndex, beamIndex) = firstTimes;
        angleInitialPosition_m(:, :, parameterIndex, beamIndex) = ...
            initialPositions;
    end
    fprintf("  Angle case %d/%d complete.\n", ...
        parameterIndex, numel(trajectoryAngles_deg));
end

independent = struct();
independent.version = "1.1";
independent.generatedAt = string(datetime("now", "Format", ...
    "yyyy-MM-dd HH:mm:ss Z"));
independent.rngSeed = rngSeed;
independent.beamNames = beamNames;
independent.beamLabels = beamLabels;
independent.numberOfTrials = numberOfTrials;
independent.system = system;
independent.metricDefinition = ...
    "D_max=max_t sqrt(sum_window(Nsig^2/(Nsig+Nbackscatter+Nbackground+Ndark)))";
independent.randomizationDefinition = ...
    "The target starts outside the beam; its initial along-line coordinate is uniform over one inter-pulse travel interval.";

independent.distance = struct( ...
    "fractionsOfGeometricRadius", distanceFractions, ...
    "offset_m", distanceFractions .* system.geometricRadius_m, ...
    "lineTrajectoryAngle_deg", 30, ...
    "initialAlongPosition_m", distanceInitialAlong_m, ...
    "initialPosition_m", distanceInitialPosition_m, ...
    "maxDetectionMetric", distanceMetric, ...
    "detected", distanceDetected, ...
    "firstDetectionTime_s", distanceFirstDetection_s);

independent.angle = struct( ...
    "values_deg", trajectoryAngles_deg, ...
    "circularBeamOffsetFraction", 2/5, ...
    "circularBeamOffset_m", (2/5) * system.geometricRadius_m, ...
    "lineIntersectionOffset_m", 0, ...
    "initialAlongPosition_m", angleInitialAlong_m, ...
    "initialPosition_m", angleInitialPosition_m, ...
    "maxDetectionMetric", angleMetric, ...
    "detected", angleDetected, ...
    "firstDetectionTime_s", angleFirstDetection_s);

outputFile = fullfile(scriptDir, "independent_data.mat");
save(outputFile, "independent", "-v7.3");
fprintf("Saved %s\n", outputFile);

function initialAlong_m = sampleInitialAlongPositions( ...
        initializationInterval_m, numberOfTrials, numberOfCases)
intervalWidth_m = diff(initializationInterval_m);
initialAlong_m = initializationInterval_m(1) + intervalWidth_m .* ...
    rand(numberOfTrials, numberOfCases);
end

function [maximumMetric, detected, firstDetectionTime_s, ...
        initialPosition_m] = simulateTrial(system, beamName, model, ...
        initialAlong_m, studyName, parameterValue)
% Define an infinite trajectory line first, then initialize the target on an
% exterior segment of that line and propagate it at exactly 1/PRF intervals.

switch studyName
    case "distance"
        trajectoryAngle = deg2rad(30);
        distance_m = parameterValue * system.geometricRadius_m;
        if beamName == "line"
            % The trajectory intersects the line long axis at [distance, 0].
            trajectoryAnchor_m = [distance_m, 0];
        else
            % The trajectory's closest-axis point is distance_m from origin.
            trajectoryAnchor_m = distance_m .* [-sin(trajectoryAngle), ...
                cos(trajectoryAngle)];
        end
    case "angle"
        trajectoryAngle = deg2rad(parameterValue);
        if beamName == "line"
            trajectoryAnchor_m = [0, 0];
        else
            distance_m = (2/5) * system.geometricRadius_m;
            trajectoryAnchor_m = distance_m .* [-sin(trajectoryAngle), ...
                cos(trajectoryAngle)];
        end
    otherwise
        error("independent_data:UnknownStudy", ...
            "Unknown study name: %s", studyName);
end

direction = [cos(trajectoryAngle), sin(trajectoryAngle)];
initialPosition_m = trajectoryAnchor_m + initialAlong_m .* direction;
if illuminationFactor(system, beamName, initialPosition_m) > 0
    error("independent_data:InitialPositionInsideBeam", ...
        "The target initial position must be outside the active beam field.");
end

simulationDuration_s = (system.trackHalfLength_m - initialAlong_m) / ...
    system.targetSpeed_m_s;
numberOfPulses = floor(simulationDuration_s * system.prf_Hz) + 1;
sampleTime_s = (0:numberOfPulses-1).' / system.prf_Hz;
xy_m = initialPosition_m + ...
    (system.targetSpeed_m_s .* sampleTime_s) .* direction;

illumination = illuminationFactor(system, beamName, xy_m);
signalPhotons = model.signalAtPeak .* illumination;
perPulseScore = zeros(size(signalPhotons));
hasSignal = signalPhotons > 0;
perPulseScore(hasSignal) = signalPhotons(hasSignal).^2 ./ ...
    (signalPhotons(hasSignal) + model.totalNoisePhotons);

windowScore = movsum(perPulseScore, [model.windowPulses - 1, 0], ...
    "Endpoints", "shrink");
windowMetric = sqrt(max(windowScore, 0));
maximumMetric = max(windowMetric);
firstIndex = find(windowMetric >= system.snrThreshold, 1, "first");
detected = ~isempty(firstIndex);
if detected
    firstDetectionTime_s = (firstIndex - 1) / system.prf_Hz;
else
    firstDetectionTime_s = nan;
end
end

function illumination = illuminationFactor(system, beamName, xy_m)
% Average normalized illumination over the same 7-by-7 target grid used by
% MC_point_snr, MC_line_snr, and MC_ring_snr.

x = xy_m(:, 1);
y = xy_m(:, 2);
z = system.range_m;
angularSide = 2 * atan((system.targetWidth_m / 2) / z);
nodes = ((1:system.targetSampleGrid) - 0.5) ./ ...
    system.targetSampleGrid - 0.5;
illumination = zeros(size(x));

switch beamName
    case "point"
        gaussianWidth = system.wD_rad / 2;
        fovHalfAngle = 3 * gaussianWidth;
        theta = atan2(hypot(x, y), z);
        targetHalfDiagonal = atan((sqrt(2) * system.targetWidth_m / 2) / z);
        candidate = theta <= fovHalfAngle + targetHalfDiagonal;
        interior = candidate & theta + targetHalfDiagonal <= fovHalfAngle;

        % Away from the hard FOV edge, the circular Gaussian samples factor
        % into radial and tangential one-dimensional means exactly.
        if any(interior)
            thetaInterior = theta(interior);
            radial = thetaInterior + angularSide .* nodes;
            tangent = angularSide .* nodes;
            illumination(interior) = mean(exp(-2 .* ...
                (radial ./ gaussianWidth).^2), 2) .* ...
                mean(exp(-2 .* (tangent ./ gaussianWidth).^2));
        end

        boundaryIndex = find(candidate & ~interior);
        if ~isempty(boundaryIndex)
            illumination(boundaryIndex) = pointBoundaryFactor( ...
                theta(boundaryIndex), angularSide, nodes, ...
                gaussianWidth, fovHalfAngle);
        end

    case "line"
        longOffset = atan2(x, z);
        shortOffset = atan2(y, z);
        longHalfWidth = system.wD_rad / 2;
        gaussianWidth = system.wd_rad / 2;
        targetHalf = atan((system.targetWidth_m / 2) / z);
        candidate = abs(longOffset) <= longHalfWidth + targetHalf & ...
            abs(shortOffset) <= ...
            system.gaussianCutoffWidths * system.wd_rad + targetHalf;
        if any(candidate)
            sampledLong = longOffset(candidate) + angularSide .* nodes;
            sampledShort = shortOffset(candidate) + angularSide .* nodes;
            longFactor = mean(abs(sampledLong) <= longHalfWidth, 2);
            shortFactor = mean(exp(-2 .* ...
                (sampledShort ./ gaussianWidth).^2), 2);
            illumination(candidate) = longFactor .* shortFactor;
        end

    case "ring"
        ringRadius = system.wD_rad / 2;
        gaussianWidth = system.wd_rad / 2;
        theta = atan2(hypot(x, y), z);
        targetHalfDiagonal = atan((sqrt(2) * system.targetWidth_m / 2) / z);
        candidate = abs(theta - ringRadius) <= targetHalfDiagonal + ...
            system.gaussianCutoffWidths * gaussianWidth;
        candidateIndex = find(candidate);
        if ~isempty(candidateIndex)
            illumination(candidateIndex) = ringTargetFactor( ...
                theta(candidateIndex), angularSide, nodes, ...
                ringRadius, gaussianWidth);
        end
end
end

function factor = pointBoundaryFactor(theta, angularSide, nodes, ...
        gaussianWidth, fovHalfAngle)
[sampleX, sampleY] = meshgrid(nodes, nodes);
dx = angularSide .* sampleX(:).';
dy = angularSide .* sampleY(:).';
sampleTheta = sqrt((theta(:) + dx).^2 + dy.^2);
intensity = (sampleTheta <= fovHalfAngle) .* ...
    exp(-2 .* (sampleTheta ./ gaussianWidth).^2);
factor = mean(intensity, 2);
end

function factor = ringTargetFactor(theta, angularSide, nodes, ...
        ringRadius, gaussianWidth)
[sampleX, sampleY] = meshgrid(nodes, nodes);
dx = angularSide .* sampleX(:).';
dy = angularSide .* sampleY(:).';
sampleTheta = sqrt((theta(:) + dx).^2 + dy.^2);
factor = mean(exp(-2 .* ...
    ((sampleTheta - ringRadius) ./ gaussianWidth).^2), 2);
end

function model = buildPhotonModel(system, beamName)
h = 6.62607015e-34;
c = 299792458;
photonEnergy = h * c / system.wavelength_m;
receiverArea = pi * (system.receiverDiameter_m / 2)^2;
systemEfficiency = system.opticalTransmission * system.quantumEfficiency;
transmittedPhotons = system.pulseEnergy_J / photonEnergy;

switch beamName
    case "point"
        gaussianWidth = system.wD_rad / 2;
        omegaBeam = pointGaussianSolidAngle(gaussianWidth);
        fovHalfAngle = min(3 * gaussianWidth, pi);
        omegaFov = 2 * pi * (1 - cos(fovHalfAngle));
        windowPulses = system.windowPulses(1);
    case "line"
        gaussianWidth = system.wd_rad / 2;
        omegaBeam = system.wD_rad * gaussianWidth * sqrt(pi / 2);
        omegaFov = system.wD_rad * (6 * system.wd_rad);
        windowPulses = system.windowPulses(2);
    case "ring"
        ringRadius = system.wD_rad / 2;
        gaussianWidth = system.wd_rad / 2;
        omegaBeam = annularGaussianSolidAngle(ringRadius, gaussianWidth);
        fovHalfWidth = 3 * system.wd_rad;
        thetaInner = max(ringRadius - fovHalfWidth, 0);
        thetaOuter = min(ringRadius + fovHalfWidth, pi);
        omegaFov = 2 * pi * (cos(thetaInner) - cos(thetaOuter));
        windowPulses = system.windowPulses(3);
end

backgroundPhotons = system.skyRadiance_W_m2_sr_nm * ...
    system.filterBandwidth_nm * receiverArea * omegaFov * ...
    system.rangeGate_s * systemEfficiency / photonEnergy;
darkPhotons = system.darkCountRate_Hz * system.rangeGate_s;

gateRange_m = c * system.rangeGate_s / 2;
offsets_m = linspace(-gateRange_m/2, gateRange_m/2, 101);
rangeSamples_m = system.range_m + offsets_m;
integrand = system.backscatter_m_inv_sr_inv .* ...
    exp(-2 .* system.extinction_m_inv .* rangeSamples_m) ./ ...
    rangeSamples_m.^2;
backscatterPhotons = transmittedPhotons * systemEfficiency * ...
    receiverArea * trapz(offsets_m, integrand);

targetArea = system.targetWidth_m^2;
twoWayTransmission = exp(-2 * system.extinction_m_inv * system.range_m);
signalAtPeak = transmittedPhotons * (targetArea / ...
    (omegaBeam * system.range_m^2)) * system.targetReflectivity * ...
    (receiverArea / (pi * system.range_m^2)) * twoWayTransmission * ...
    systemEfficiency;

model = struct( ...
    "beamName", beamName, ...
    "omegaBeam_sr", omegaBeam, ...
    "omegaFov_sr", omegaFov, ...
    "signalAtPeak", signalAtPeak, ...
    "backscatterPhotons", backscatterPhotons, ...
    "backgroundPhotons", backgroundPhotons, ...
    "darkPhotons", darkPhotons, ...
    "totalNoisePhotons", backscatterPhotons + backgroundPhotons + darkPhotons, ...
    "windowPulses", windowPulses);
end

function omega = pointGaussianSolidAngle(gaussianWidth)
integrand = @(theta) exp(-2 .* (theta ./ gaussianWidth).^2) .* sin(theta);
omega = 2 * pi * integral(integrand, 0, pi, ...
    "RelTol", 1e-10, "AbsTol", 1e-14);
end

function omega = annularGaussianSolidAngle(ringRadius, gaussianWidth)
u0 = -sqrt(2) * ringRadius / gaussianWidth;
radialIntegral = ringRadius * gaussianWidth * sqrt(pi) / (2 * sqrt(2)) * ...
    erfc(u0) + gaussianWidth^2 / 4 * ...
    exp(-2 * (ringRadius / gaussianWidth)^2);
omega = 2 * pi * radialIntegral;
end

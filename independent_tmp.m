clc;
clearvars;
close all;

% Common-trajectory comparison of spot, line, and annular beams.
% This script does not modify independent_data.m or independent_draw.m.
% It reuses their saved physical parameters, but performs a new simulation
% in which all three beams receive exactly the same target trajectories.

scriptDir = fileparts(mfilename("fullpath"));
if strlength(scriptDir) == 0
    scriptDir = pwd;
end

sourceDataFile = fullfile(scriptDir, "independent_data.mat");
if ~isfile(sourceDataFile)
    error("independent_tmp:MissingData", ...
        "Run independent_data.m first: %s", sourceDataFile);
end

loaded = load(sourceDataFile, "independent");
data = simulateCommonTrajectoryStudy(loaded.independent);
dataFile = fullfile(scriptDir, "independent_tmp_data.mat");
save(dataFile, "data", "-v7.3");
validateIndependentData(data, dataFile);

colors = [ ...
    0.4660, 0.6740, 0.1880; ... % Spot
    0.0000, 0.4470, 0.7410; ... % Line
    0.8500, 0.3250, 0.0980];    % Annular
markers = ["o", "s", "^"];
beamOffsets = [-0.20, 0, 0.20];
cloudWidth = 0.13;
fontSize = 13;
axesLineWidth = 2;
curveLineWidth = 2;

distanceX = 1:numel(data.distance.fractionsOfGeometricRadius);
angleX = 1:numel(data.angle.values_deg);
distanceTickLabels = {"0", "1/5", "2/5", "3/5", "4/5", "1"};
angleTickLabels = compose("%g^\\circ", data.angle.values_deg);

%% Figure 1: detection probability
probabilityFigure = figure("Color", "w", ...
    "Position", [60, 80, 1200, 520]);
probabilityLayout = tiledlayout(probabilityFigure, 1, 2, ...
    "TileSpacing", "compact", "Padding", "compact");

probabilityDistanceAxes = nexttile(probabilityLayout, 1);
probabilityLegendHandles = drawProbabilityPanel( ...
    probabilityDistanceAxes, distanceX, data.distance.detected, ...
    data.beamLabels, colors, markers, beamOffsets, data.numberOfTrials, ...
    fontSize, axesLineWidth, curveLineWidth);
xticks(probabilityDistanceAxes, distanceX);
xticklabels(probabilityDistanceAxes, distanceTickLabels);
xlabel(probabilityDistanceAxes, ...
    "Normalized center offset, {\itd}_{\rmc}/{\itr}_{\rmD}", ...
    "Interpreter", "tex");
ylabel(probabilityDistanceAxes, "Detection probability (%)");
title(probabilityDistanceAxes, "(a) Trajectory offset", ...
    "FontWeight", "normal");

probabilityAngleAxes = nexttile(probabilityLayout, 2);
drawProbabilityPanel(probabilityAngleAxes, angleX, data.angle.detected, ...
    data.beamLabels, colors, markers, beamOffsets, data.numberOfTrials, ...
    fontSize, axesLineWidth, curveLineWidth);
xticks(probabilityAngleAxes, angleX);
xticklabels(probabilityAngleAxes, angleTickLabels);
xlabel(probabilityAngleAxes, "Trajectory-beam angle", ...
    "Interpreter", "tex");
title(probabilityAngleAxes, "(b) Trajectory angle", ...
    "FontWeight", "normal");

probabilityLegend = legend(probabilityDistanceAxes, ...
    probabilityLegendHandles, cellstr(data.beamLabels), ...
    "Orientation", "horizontal", "Box", "off", ...
    "FontName", "Helvetica", "FontSize", fontSize);
probabilityLegend.Layout.Tile = "north";
probabilityLegend.NumColumns = numel(data.beamLabels);
probabilityLegend.ItemTokenSize = [45, 18];

probabilityPng = fullfile(scriptDir, "independent_tmp_probability.png");
probabilityFig = fullfile(scriptDir, "independent_tmp_probability.fig");
exportgraphics(probabilityFigure, probabilityPng, "Resolution", 300);
savefig(probabilityFigure, probabilityFig);

%% Figure 2: maximum sliding-window statistic D_max
dmaxFigure = figure("Color", "w", "Position", [80, 100, 1200, 520]);
dmaxLayout = tiledlayout(dmaxFigure, 1, 2, ...
    "TileSpacing", "compact", "Padding", "compact");

dmaxDistanceAxes = nexttile(dmaxLayout, 1);
dmaxLegendHandles = drawDmaxPanel(dmaxDistanceAxes, distanceX, ...
    data.distance.maxDetectionMetric, data.beamLabels, colors, ...
    beamOffsets, cloudWidth, fontSize, axesLineWidth, curveLineWidth);
xticks(dmaxDistanceAxes, distanceX);
xticklabels(dmaxDistanceAxes, distanceTickLabels);
xlabel(dmaxDistanceAxes, ...
    "Normalized center offset, {\itd}_{\rmc}/{\itr}_{\rmD}", ...
    "Interpreter", "tex");
ylabel(dmaxDistanceAxes, ...
    "Maximum sliding-window statistic, {\itD}_{\rmmax}", ...
    "Interpreter", "tex");
title(dmaxDistanceAxes, "(a) Trajectory offset", ...
    "FontWeight", "normal");

dmaxAngleAxes = nexttile(dmaxLayout, 2);
drawDmaxPanel(dmaxAngleAxes, angleX, data.angle.maxDetectionMetric, ...
    data.beamLabels, colors, beamOffsets, cloudWidth, fontSize, ...
    axesLineWidth, curveLineWidth);
xticks(dmaxAngleAxes, angleX);
xticklabels(dmaxAngleAxes, angleTickLabels);
xlabel(dmaxAngleAxes, "Trajectory-beam angle", "Interpreter", "tex");
title(dmaxAngleAxes, "(b) Trajectory angle", ...
    "FontWeight", "normal");

maximumY = max([dmaxDistanceAxes.YLim(2), dmaxAngleAxes.YLim(2)]);
ylim(dmaxDistanceAxes, [0, maximumY]);
ylim(dmaxAngleAxes, [0, maximumY]);
yline(dmaxDistanceAxes, data.system.snrThreshold, "--", ...
    "Threshold", "Color", [0.25, 0.25, 0.25], ...
    "LineWidth", 1.5, "LabelHorizontalAlignment", "left", ...
    "HandleVisibility", "off");
yline(dmaxAngleAxes, data.system.snrThreshold, "--", ...
    "Color", [0.25, 0.25, 0.25], "LineWidth", 1.5, ...
    "HandleVisibility", "off");

dmaxLegend = legend(dmaxDistanceAxes, dmaxLegendHandles, ...
    cellstr(data.beamLabels), "Orientation", "horizontal", ...
    "Box", "off", "FontName", "Helvetica", "FontSize", fontSize);
dmaxLegend.Layout.Tile = "north";
dmaxLegend.NumColumns = numel(data.beamLabels);
dmaxLegend.ItemTokenSize = [45, 18];

dmaxPng = fullfile(scriptDir, "independent_tmp_dmax.png");
dmaxFig = fullfile(scriptDir, "independent_tmp_dmax.fig");
exportgraphics(dmaxFigure, dmaxPng, "Resolution", 300);
savefig(dmaxFigure, dmaxFig);

fprintf("Saved %s\nSaved %s\n", probabilityPng, probabilityFig);
fprintf("Saved %s\nSaved %s\n", dmaxPng, dmaxFig);

distanceProbability = squeeze(mean(data.distance.detected, 1));
angleProbability = squeeze(mean(data.angle.detected, 1));
fprintf("Mean detection probability over offset cases [Spot Line Annular]:\n");
disp(mean(distanceProbability, 1));
fprintf("Mean detection probability over angle cases [Spot Line Annular]:\n");
disp(mean(angleProbability, 1));

function legendHandles = drawProbabilityPanel(ax, categoryX, detected, ...
        beamLabels, colors, markers, beamOffsets, numberOfTrials, ...
        fontSize, axesLineWidth, curveLineWidth)
hold(ax, "on");
numberOfBeams = numel(beamLabels);
legendHandles = gobjects(1, numberOfBeams);

for beamIndex = 1:numberOfBeams
    successes = squeeze(sum(detected(:, :, beamIndex), 1));
    probability = successes ./ numberOfTrials;
    [lowerBound, upperBound] = wilsonInterval(successes, numberOfTrials);
    shiftedX = categoryX + beamOffsets(beamIndex);
    legendHandles(beamIndex) = errorbar(ax, shiftedX, ...
        100 .* probability, ...
        100 .* (probability - lowerBound), ...
        100 .* (upperBound - probability), ...
        "-" + markers(beamIndex), ...
        "Color", colors(beamIndex, :), ...
        "MarkerFaceColor", colors(beamIndex, :), ...
        "MarkerSize", 5, "LineWidth", curveLineWidth, ...
        "CapSize", 7, "DisplayName", beamLabels(beamIndex));
end

grid(ax, "on");
box(ax, "on");
xlim(ax, [0.55, numel(categoryX) + 0.45]);
ylim(ax, [0, 105]);
yticks(ax, 0:20:100);
set(ax, "FontName", "Helvetica", "FontSize", fontSize, ...
    "LineWidth", axesLineWidth, "Layer", "top");
set([ax.XLabel, ax.YLabel], ...
    "FontName", "Helvetica", "FontSize", fontSize);
end

function [lowerBound, upperBound] = wilsonInterval(successes, trials)
% Two-sided 95% Wilson score interval for a binomial proportion.
z = 1.95996398454005;
probability = successes ./ trials;
denominator = 1 + z^2 / trials;
center = (probability + z^2 / (2 * trials)) ./ denominator;
halfWidth = z .* sqrt(probability .* (1 - probability) ./ trials + ...
    z^2 / (4 * trials^2)) ./ denominator;
lowerBound = max(0, center - halfWidth);
upperBound = min(1, center + halfWidth);
end

function legendHandles = drawDmaxPanel(ax, categoryX, metric, ...
        beamLabels, colors, beamOffsets, cloudWidth, fontSize, ...
        axesLineWidth, curveLineWidth)
hold(ax, "on");
numberOfBeams = numel(beamLabels);
legendHandles = gobjects(1, numberOfBeams);
medianValues = nan(numel(categoryX), numberOfBeams);

for beamIndex = 1:numberOfBeams
    for categoryIndex = 1:numel(categoryX)
        center = categoryX(categoryIndex) + beamOffsets(beamIndex);
        values = metric(:, categoryIndex, beamIndex);
        drawRaincloud(ax, center, values, colors(beamIndex, :), ...
            cloudWidth, categoryIndex + 17 * beamIndex);
        medianValues(categoryIndex, beamIndex) = median(values, "omitnan");
    end
    shiftedX = categoryX + beamOffsets(beamIndex);
    legendHandles(beamIndex) = plot(ax, shiftedX, ...
        medianValues(:, beamIndex), "-o", ...
        "Color", colors(beamIndex, :), ...
        "LineWidth", curveLineWidth, "MarkerSize", 4.5, ...
        "MarkerFaceColor", colors(beamIndex, :), ...
        "DisplayName", beamLabels(beamIndex));
end

grid(ax, "on");
box(ax, "on");
xlim(ax, [0.55, numel(categoryX) + 0.45]);
set(ax, "FontName", "Helvetica", "FontSize", fontSize, ...
    "LineWidth", axesLineWidth, "Layer", "top");
set([ax.XLabel, ax.YLabel], ...
    "FontName", "Helvetica", "FontSize", fontSize);
end

function drawRaincloud(ax, center, values, color, cloudWidth, phase)
values = values(isfinite(values));
if isempty(values)
    return;
end

[density, yGrid] = simpleKDE(values);
if max(density) > 0
    density = cloudWidth .* density ./ max(density);
end
patch(ax, [center - density, repmat(center, size(density))], ...
    [yGrid, fliplr(yGrid)], color, ...
    "FaceAlpha", 0.20, "EdgeColor", "none", ...
    "HandleVisibility", "off");

index = (1:numel(values)).';
jitter01 = mod(index * 0.618033988749895 + phase * 0.137, 1);
rainX = center + 0.015 + cloudWidth * jitter01;
scatter(ax, rainX, values, 8, color, "filled", ...
    "MarkerFaceAlpha", 0.22, "MarkerEdgeAlpha", 0.12, ...
    "HandleVisibility", "off");

quartiles = localQuantile(values, [0.25, 0.50, 0.75]);
plot(ax, [center, center], quartiles([1, 3]), "-", ...
    "Color", color, "LineWidth", 2.6, "HandleVisibility", "off");
plot(ax, center, quartiles(2), "o", "MarkerSize", 4.5, ...
    "MarkerFaceColor", "w", "MarkerEdgeColor", color, ...
    "LineWidth", 1.4, "HandleVisibility", "off");
end

function [density, yGrid] = simpleKDE(values)
values = values(:);
numberOfValues = numel(values);
sigma = std(values);
quartiles = localQuantile(values, [0.25, 0.75]);
robustSigma = min(sigma, (quartiles(2) - quartiles(1)) / 1.34);
if ~isfinite(robustSigma) || robustSigma <= 0
    robustSigma = max(abs(mean(values)), 1) * 1e-3;
end
bandwidth = max(0.9 * robustSigma * numberOfValues^(-1/5), ...
    max(abs(mean(values)), 1) * 1e-6);
yPadding = max(3 * bandwidth, ...
    0.03 * max(max(values) - min(values), bandwidth));
yGrid = linspace(min(values) - yPadding, ...
    max(values) + yPadding, 120);
standardized = (yGrid(:) - values.') ./ bandwidth;
density = mean(exp(-0.5 .* standardized.^2), 2).' ./ ...
    (bandwidth * sqrt(2 * pi));
end

function q = localQuantile(values, probabilities)
values = sort(values(:));
numberOfValues = numel(values);
positions = 1 + (numberOfValues - 1) .* probabilities(:).';
lowerIndex = floor(positions);
upperIndex = ceil(positions);
weights = positions - lowerIndex;
q = values(lowerIndex).' .* (1 - weights) + ...
    values(upperIndex).' .* weights;
end

function validateIndependentData(data, dataFile)
required = ["beamLabels", "numberOfTrials", "system", "distance", "angle"];
missing = required(~isfield(data, required));
if ~isempty(missing)
    error("independent_tmp:InvalidData", ...
        "Missing fields in %s: %s", dataFile, strjoin(missing, ", "));
end

requiredResultFields = ["detected", "maxDetectionMetric"];
for studyName = ["distance", "angle"]
    study = data.(studyName);
    missingResultFields = requiredResultFields( ...
        ~isfield(study, requiredResultFields));
    if ~isempty(missingResultFields)
        error("independent_tmp:InvalidData", ...
            "Missing %s fields in %s: %s", studyName, dataFile, ...
            strjoin(missingResultFields, ", "));
    end
    if size(study.detected, 1) ~= data.numberOfTrials || ...
            ~isequal(size(study.detected), size(study.maxDetectionMetric))
        error("independent_tmp:InvalidDimensions", ...
            "%s result dimensions in %s are inconsistent.", ...
            studyName, dataFile);
    end
end
end

function data = simulateCommonTrajectoryStudy(source)
% Recompute both studies with common physical trajectories for all beams.

system = source.system;
beamNames = source.beamNames;
beamLabels = source.beamLabels;
numberOfTrials = source.numberOfTrials;
distanceFractions = source.distance.fractionsOfGeometricRadius(:);
trajectoryAngles_deg = source.angle.values_deg(:);

rng(source.rngSeed + 101, "twister");
photonModels = cell(1, numel(beamNames));
for beamIndex = 1:numel(beamNames)
    photonModels{beamIndex} = buildPhotonModel(system, beamNames(beamIndex));
end

distanceInitialPhase_m = system.interPulseTravel_m .* ...
    rand(numberOfTrials, numel(distanceFractions));
distanceMetric = nan(numberOfTrials, numel(distanceFractions), ...
    numel(beamNames));
distanceDetected = false(size(distanceMetric));

fprintf("Running common-trajectory offset study...\n");
for parameterIndex = 1:numel(distanceFractions)
    fraction = distanceFractions(parameterIndex);
    phases = distanceInitialPhase_m(:, parameterIndex);
    for beamIndex = 1:numel(beamNames)
        for trialIndex = 1:numberOfTrials
            [distanceMetric(trialIndex, parameterIndex, beamIndex), ...
                distanceDetected(trialIndex, parameterIndex, beamIndex)] = ...
                simulateTrial(system, beamNames(beamIndex), ...
                photonModels{beamIndex}, phases(trialIndex), ...
                30, fraction);
        end
    end
end

angleInitialPhase_m = system.interPulseTravel_m .* ...
    rand(numberOfTrials, numel(trajectoryAngles_deg));
angleSignedOffsetFraction = 2 .* ...
    rand(numberOfTrials, numel(trajectoryAngles_deg)) - 1;
angleMetric = nan(numberOfTrials, numel(trajectoryAngles_deg), ...
    numel(beamNames));
angleDetected = false(size(angleMetric));

fprintf("Running common-trajectory angle study...\n");
for parameterIndex = 1:numel(trajectoryAngles_deg)
    angle_deg = trajectoryAngles_deg(parameterIndex);
    phases = angleInitialPhase_m(:, parameterIndex);
    signedOffsets = angleSignedOffsetFraction(:, parameterIndex);
    for beamIndex = 1:numel(beamNames)
        for trialIndex = 1:numberOfTrials
            [angleMetric(trialIndex, parameterIndex, beamIndex), ...
                angleDetected(trialIndex, parameterIndex, beamIndex)] = ...
                simulateTrial(system, beamNames(beamIndex), ...
                photonModels{beamIndex}, phases(trialIndex), ...
                angle_deg, signedOffsets(trialIndex));
        end
    end
end

data = struct();
data.version = "common-trajectory-1.0";
data.generatedAt = string(datetime("now", "Format", ...
    "yyyy-MM-dd HH:mm:ss Z"));
data.rngSeed = source.rngSeed + 101;
data.beamNames = beamNames;
data.beamLabels = beamLabels;
data.numberOfTrials = numberOfTrials;
data.system = system;
data.trajectoryDefinition = ...
    "All beam types use identical trajectory direction, closest approach, and pulse phase.";
data.distance = struct( ...
    "fractionsOfGeometricRadius", distanceFractions, ...
    "offset_m", distanceFractions .* system.geometricRadius_m, ...
    "trajectoryAngle_deg", 30, ...
    "initialPhase_m", distanceInitialPhase_m, ...
    "maxDetectionMetric", distanceMetric, ...
    "detected", distanceDetected);
data.angle = struct( ...
    "values_deg", trajectoryAngles_deg, ...
    "signedOffsetFraction", angleSignedOffsetFraction, ...
    "offsetDistribution", "Uniform(-r_D, r_D)", ...
    "initialPhase_m", angleInitialPhase_m, ...
    "maxDetectionMetric", angleMetric, ...
    "detected", angleDetected);
end

function [maximumMetric, detected] = simulateTrial( ...
        system, beamName, model, initialPhase_m, ...
        trajectoryAngle_deg, signedDistanceFraction)
% Beam axis is +z. The same tangent-plane trajectory is used for all beams.

trajectoryAngle = deg2rad(trajectoryAngle_deg);
distance_m = signedDistanceFraction * system.geometricRadius_m;

direction = [cos(trajectoryAngle), sin(trajectoryAngle)];
normal = [-sin(trajectoryAngle), cos(trajectoryAngle)];
closestPoint = distance_m .* normal;
startAlong_m = -system.trackHalfLength_m + initialPhase_m;
numberOfPulses = floor((system.trackHalfLength_m - startAlong_m) ./ ...
    system.interPulseTravel_m) + 1;
along_m = startAlong_m + (0:numberOfPulses-1).' .* ...
    system.interPulseTravel_m;
xy_m = closestPoint + along_m .* direction;

illumination = illuminationFactor(system, beamName, xy_m);
signalPhotons = model.signalAtPeak .* illumination;
perPulseScore = zeros(size(signalPhotons));
hasSignal = signalPhotons > 0;
perPulseScore(hasSignal) = signalPhotons(hasSignal).^2 ./ ...
    (signalPhotons(hasSignal) + model.totalNoisePhotons);
windowScore = movsum(perPulseScore, [model.windowPulses - 1, 0], ...
    "Endpoints", "shrink");
maximumMetric = sqrt(max(windowScore));
detected = maximumMetric >= system.snrThreshold;
end

function illumination = illuminationFactor(system, beamName, xy_m)
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
        thetaInterior = theta(interior);
        radial = thetaInterior + angularSide .* nodes;
        tangent = angularSide .* nodes;
        illumination(interior) = mean(exp(-2 .* ...
            (radial ./ gaussianWidth).^2), 2) .* ...
            mean(exp(-2 .* (tangent ./ gaussianWidth).^2));
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
    "signalAtPeak", signalAtPeak, ...
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

clc;
clearvars;
close all;

% Raincloud summary of the independent stationary-beam trajectory study.

scriptDir = fileparts(mfilename("fullpath"));
if strlength(scriptDir) == 0
    scriptDir = pwd;
end
dataFile = fullfile(scriptDir, "independent_data.mat");
if ~isfile(dataFile)
    error("independent_draw:MissingData", ...
        "Run independent_data.m first: %s", dataFile);
end

loaded = load(dataFile, "independent");
validateIndependentData(loaded.independent, dataFile);
data = loaded.independent;

colors = [ ...
    0.4660, 0.6740, 0.1880; ... % Spot
    0.0000, 0.4470, 0.7410; ... % Line
    0.8500, 0.3250, 0.0980];    % Annular
beamOffsets = [-0.22, 0, 0.22];
cloudWidth = 0.14;

figureHandle = figure("Color", "w", "Position", [60, 80, 1200, 520]);
layout = tiledlayout(figureHandle, 1, 2, ...
    "TileSpacing", "compact", "Padding", "compact");

distanceAxes = nexttile(layout, 1);
distanceX = 1:numel(data.distance.fractionsOfGeometricRadius);
drawPanel(distanceAxes, distanceX, data.distance.maxDetectionMetric, ...
    data.beamLabels, colors, beamOffsets, cloudWidth);
xticks(distanceAxes, distanceX);
xticklabels(distanceAxes, {"0", "1/5", "2/5", "3/5", "4/5", "1"});
xlabel(distanceAxes, ...
    "Normalized center offset, {\itd}_{\rmc}/{\itr}_{\rmD}", ...
    "Interpreter", "tex");
ylabel(distanceAxes, ...
    "Maximum sliding-window statistic, {\itD}_{\rmmax}", ...
    "Interpreter", "tex");
title(distanceAxes, "(a) Trajectory offset", ...
    "FontWeight", "normal");

angleAxes = nexttile(layout, 2);
angleX = 1:numel(data.angle.values_deg);
legendHandles = drawPanel(angleAxes, angleX, ...
    data.angle.maxDetectionMetric, data.beamLabels, colors, ...
    beamOffsets, cloudWidth);
xticks(angleAxes, angleX);
xticklabels(angleAxes, compose("%g^\\circ", data.angle.values_deg));
xlabel(angleAxes, "Trajectory-beam angle", "Interpreter", "tex");
title(angleAxes, "(b) Trajectory angle", ...
    "FontWeight", "normal");

maximumY = max([distanceAxes.YLim(2), angleAxes.YLim(2)]);
ylim(distanceAxes, [0, maximumY]);
ylim(angleAxes, [0, maximumY]);
yline(distanceAxes, data.system.snrThreshold, "--", ...
    "Threshold", "Color", [0.25, 0.25, 0.25], ...
    "LineWidth", 1.5, "LabelHorizontalAlignment", "left", ...
    "HandleVisibility", "off");
yline(angleAxes, data.system.snrThreshold, "--", ...
    "Color", [0.25, 0.25, 0.25], "LineWidth", 1.5, ...
    "HandleVisibility", "off");

legendHandle = legend(angleAxes, legendHandles, cellstr(data.beamLabels), ...
    "Orientation", "horizontal", "Box", "off", ...
    "FontName", "Helvetica", "FontSize", 13);
legendHandle.Layout.Tile = "north";
legendHandle.NumColumns = numel(data.beamLabels);
legendHandle.ItemTokenSize = [45, 18];

outputPng = fullfile(scriptDir, "independent.png");
outputFig = fullfile(scriptDir, "independent.fig");
exportgraphics(figureHandle, outputPng, "Resolution", 300);
savefig(figureHandle, outputFig);
fprintf("Saved %s\nSaved %s\n", outputPng, outputFig);

function legendHandles = drawPanel(ax, categoryX, metric, beamLabels, ...
        colors, beamOffsets, cloudWidth)
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
        "Color", colors(beamIndex, :), "LineWidth", 2, ...
        "MarkerSize", 4.5, "MarkerFaceColor", colors(beamIndex, :), ...
        "DisplayName", beamLabels(beamIndex));
end

grid(ax, "on");
box(ax, "on");
xlim(ax, [0.55, numel(categoryX) + 0.45]);
set(ax, "FontName", "Helvetica", "FontSize", 13, "LineWidth", 2, ...
    "Layer", "top");
set([ax.XLabel, ax.YLabel], ...
    "FontName", "Helvetica", "FontSize", 13);
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

% Deterministic low-discrepancy jitter keeps redraws identical without
% consuming or changing the Monte Carlo RNG stream.
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
n = numel(values);
sigma = std(values);
q = localQuantile(values, [0.25, 0.75]);
robustSigma = min(sigma, (q(2) - q(1)) / 1.34);
if ~isfinite(robustSigma) || robustSigma <= 0
    robustSigma = max(abs(mean(values)), 1) * 1e-3;
end
bandwidth = max(0.9 * robustSigma * n^(-1/5), ...
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
n = numel(values);
positions = 1 + (n - 1) .* probabilities(:).';
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
    error("independent_draw:InvalidData", ...
        "Missing fields in %s: %s", dataFile, strjoin(missing, ", "));
end
if size(data.distance.maxDetectionMetric, 1) ~= data.numberOfTrials || ...
        size(data.angle.maxDetectionMetric, 1) ~= data.numberOfTrials
    error("independent_draw:InvalidTrialCount", ...
        "Metric arrays in %s do not match numberOfTrials.", dataFile);
end
end

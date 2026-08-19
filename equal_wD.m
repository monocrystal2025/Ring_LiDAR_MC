clc;
clearvars;
close all;

% Equal-performance minimum-wD comparison for spot, line, and annular beams.
% Each heat-map cell represents a common performance requirement (t, P_D(t)).
% The displayed quantity is
%   kappa = minimum annular wD / minimum line wD.
% Only cells attainable by both line and annular beams are displayed. A value
% below one means that the annular beam reaches the common performance
% requirement with a smaller wD than the line beam; a value above one means
% that the line beam reaches it with a smaller wD.

scriptDir = string(fileparts(mfilename("fullpath")));
if strlength(scriptDir) == 0
    scriptDir = string(pwd);
end

eaResultDir = "D:\lzx\MC_RESULTS";
roiResultDir = fullfile(scriptDir, "MC_ROI_RESULTS_SNR_NEW");
assert(isfolder(eaResultDir), "EA result folder not found: %s", eaResultDir);
assert(isfolder(roiResultDir), "ROI result folder not found: %s", roiResultDir);

% Analysis settings. The small width wd is fixed at 1 mrad, as in the
% existing MC_overall_plot.m and P_D_T.m result summaries.
smallWdMrad = 1;
eaColumnIndex = 1;
timeThresholds_s = 0:1:60;
targetProbabilities = 0:0.02:1;

scenarioInfo = struct( ...
    "name", {"EA", "RA", "SP"}, ...
    "title", {"ROE", "ROI-RA", "ROI-SP"}, ...
    "resultDir", {eaResultDir, roiResultDir, roiResultDir}, ...
    "path", {"", "RA", "SP"}, ...
    "columnIndex", {eaColumnIndex, 1, 1});

numberOfScenarios = numel(scenarioInfo);
results = repmat(struct( ...
    "Kappa", [], ...
    "PointMinWD_mrad", [], ...
    "LineMinWD_mrad", [], ...
    "RingMinWD_mrad", [], ...
    "Series", []), 1, numberOfScenarios);
allRows = table();

for scenarioIndex = 1:numberOfScenarios
    info = scenarioInfo(scenarioIndex);
    fprintf("Loading %s results from %s ...\n", info.title, info.resultDir);

    series = loadScenarioSeries(info, smallWdMrad, timeThresholds_s);
    [kappa, pointMinimum, lineMinimum, ringMinimum] = ...
        equalPerformanceMinima(series, targetProbabilities);

    results(scenarioIndex).Kappa = kappa;
    results(scenarioIndex).PointMinWD_mrad = pointMinimum;
    results(scenarioIndex).LineMinWD_mrad = lineMinimum;
    results(scenarioIndex).RingMinWD_mrad = ringMinimum;
    results(scenarioIndex).Series = series;

    [timeGrid, probabilityGrid] = meshgrid( ...
        timeThresholds_s, targetProbabilities);
    rows = table( ...
        repmat(string(info.name), numel(timeGrid), 1), ...
        timeGrid(:), probabilityGrid(:), ...
        pointMinimum(:), lineMinimum(:), ringMinimum(:), kappa(:), ...
        'VariableNames', ["Scenario", "FirstWarningThreshold_s", ...
        "TargetPD", "PointMinWD_mrad", "LineMinWD_mrad", ...
        "RingMinWD_mrad", "KappaRingLine"]);
    allRows = [allRows; rows]; %#ok<AGROW>
end

% Report the finite kappa range represented by each heat-map subplot.
fprintf("\nKappa ranges (finite heat-map cells only):\n");
for scenarioIndex = 1:numberOfScenarios
    finiteKappa = results(scenarioIndex).Kappa( ...
        isfinite(results(scenarioIndex).Kappa));
    if isempty(finiteKappa)
        fprintf("  %s: no finite kappa values\n", ...
            scenarioInfo(scenarioIndex).title);
    else
        fprintf("  %s: min = %.6f, max = %.6f\n", ...
            scenarioInfo(scenarioIndex).title, ...
            min(finiteKappa), max(finiteKappa));
    end
end
fprintf("\n");

figureHandle = figure("Color", "w", "Position", [80, 80, 710, 660]);
layout = tiledlayout(figureHandle, 2, 2, ...
    "TileSpacing", "compact", "Padding", "compact");

% Tile (1,1) contains two coordinated views explaining one heat-map cell.
methodLayout = tiledlayout(layout, 1, 2, ...
    "TileSpacing", "compact", "Padding", "compact");
methodLayout.Layout.Tile = 1;
drawMethodSchematic(methodLayout, results(1).Series, timeThresholds_s);

heatmapAxes = gobjects(1, numberOfScenarios);
tileNumbers = [2, 3, 4];
for scenarioIndex = 1:numberOfScenarios
    heatmapAxes(scenarioIndex) = nexttile(layout, tileNumbers(scenarioIndex));
    plotKappaHeatmap(heatmapAxes(scenarioIndex), ...
        timeThresholds_s, targetProbabilities, ...
        results(scenarioIndex).Kappa, scenarioInfo(scenarioIndex).title);
end

colormap(figureHandle, ratioColormap(513));
set(heatmapAxes, "CLim", [0.5, 1.5]);

% One shared colorbar below the entire 2-by-2 layout.
colorbarHandle = colorbar(heatmapAxes(1));
colorbarHandle.Layout.Tile = "south";
colorbarHandle.Limits = [0.5, 1.5];
colorbarHandle.Ticks = 0.5:0.1:1.5;
colorbarHandle.Label.String = ...
    "\kappa = w_{D, min}^{Annular} / w_{D, min}^{Line}";
colorbarHandle.Label.Interpreter = "tex";
colorbarHandle.Label.FontSize = 12;
colorbarHandle.FontName = "Helvetica";
colorbarHandle.FontSize = 12;


outputBase = fullfile(scriptDir, "equal_wD");
outputPng = outputBase + ".png";
outputFig = outputBase + ".fig";
outputCsv = outputBase + ".csv";
outputMat = outputBase + ".mat";

writetable(allRows, outputCsv);
save(outputMat, "results", "scenarioInfo", "timeThresholds_s", ...
    "targetProbabilities", "smallWdMrad", "eaColumnIndex");
exportgraphics(figureHandle, outputPng, "Resolution", 300);
savefig(figureHandle, outputFig);

fprintf("Saved PNG: %s\n", outputPng);
fprintf("Saved FIG: %s\n", outputFig);
fprintf("Saved CSV: %s\n", outputCsv);
fprintf("Saved MAT: %s\n", outputMat);

function series = loadScenarioSeries(info, smallWdMrad, timeThresholds_s)
beamTypes = ["POINT", "LINE", "RING"];
series = repmat(struct("Beam", "", "WD_mrad", [], ...
    "PD", [], "SourceFiles", strings(0, 1)), 1, numel(beamTypes));

for beamIndex = 1:numel(beamTypes)
    beamType = beamTypes(beamIndex);
    [fileExpression, detectVariable, firstTimeVariable] = ...
        resultDefinition(info, beamType, smallWdMrad);
    files = dir(fullfile(info.resultDir, "*.mat"));
    fileNames = string({files.name});
    matching = ~cellfun("isempty", regexp(cellstr(fileNames), ...
        fileExpression, "once"));
    files = files(matching);
    if isempty(files)
        error("equal_wD:NoFiles", ...
            "No %s result files found for %s in %s.", ...
            beamType, info.title, info.resultDir);
    end

    wD = nan(numel(files), 1);
    sourceFiles = strings(numel(files), 1);
    for fileIndex = 1:numel(files)
        token = regexp(files(fileIndex).name, fileExpression, ...
            "tokens", "once");
        wD(fileIndex) = str2double(token{1});
        sourceFiles(fileIndex) = string(files(fileIndex).name);
    end
    [wD, order] = sort(wD);
    files = files(order);
    sourceFiles = sourceFiles(order);

    % Guard against accidentally mixing duplicate results for the same wD.
    if numel(unique(wD)) ~= numel(wD)
        error("equal_wD:DuplicateWD", ...
            "Duplicate %s wD values found for %s.", beamType, info.title);
    end

    probability = nan(numel(timeThresholds_s), numel(files));
    for fileIndex = 1:numel(files)
        resultFile = fullfile(files(fileIndex).folder, files(fileIndex).name);
        data = load(resultFile, detectVariable, firstTimeVariable);
        if ~isfield(data, detectVariable) || ...
                ~isfield(data, firstTimeVariable)
            error("equal_wD:MissingVariable", ...
                "Expected variables are missing in %s.", resultFile);
        end

        detect = data.(detectVariable);
        firstTime = data.(firstTimeVariable);
        if size(detect, 2) < info.columnIndex || ...
                size(firstTime, 2) < info.columnIndex
            error("equal_wD:MissingColumn", ...
                "Column %d is unavailable in %s.", ...
                info.columnIndex, resultFile);
        end
        detect = detect(:, info.columnIndex);
        firstTime = firstTime(:, info.columnIndex);
        if numel(detect) ~= numel(firstTime)
            error("equal_wD:SizeMismatch", ...
                "Detection and first-warning arrays differ in %s.", resultFile);
        end

        valid = detect ~= 0 & isfinite(firstTime) & firstTime >= 0;
        detectionTimes = firstTime(valid);
        upperEdges = timeThresholds_s + eps(timeThresholds_s);
        binCounts = histcounts(detectionTimes, ...
            [-inf, upperEdges, inf]);
        probability(:, fileIndex) = cumsum( ...
            binCounts(1:numel(timeThresholds_s))).' ./ numel(detect);
    end

    series(beamIndex).Beam = beamType;
    series(beamIndex).WD_mrad = wD;
    series(beamIndex).PD = probability;
    series(beamIndex).SourceFiles = sourceFiles;
end
end

function [fileExpression, detectVariable, firstTimeVariable] = ...
        resultDefinition(info, beamType, smallWdMrad)
beamCode = extractBetween(beamType, 1, 1);
numberPattern = "([0-9]+(?:\.[0-9]+)?)";
smallWdText = regexptranslate("escape", num2str(smallWdMrad, "%g"));

if info.name == "EA"
    prefix = "MC_1par_EA_" + beamType + "_D";
    pathSuffix = "";
else
    prefix = "MC_1par_ROI_NEW_" + beamType + "_" + info.path + "_D";
    pathSuffix = "_" + info.path;
end

if beamType == "POINT"
    fileExpression = "^" + prefix + numberPattern + "mrad\.mat$";
else
    fileExpression = "^" + prefix + numberPattern + "d" + ...
        smallWdText + "mrad\.mat$";
end
detectVariable = "detect_" + beamCode + pathSuffix;
firstTimeVariable = "first_time_" + beamCode + pathSuffix;
end

function [kappa, pointMinimum, lineMinimum, ringMinimum] = ...
        equalPerformanceMinima(series, targetProbabilities)
numberOfProbabilities = numel(targetProbabilities);
numberOfTimes = size(series(1).PD, 1);
minimumWD = nan(3, numberOfProbabilities, numberOfTimes);

for beamIndex = 1:3
    for timeIndex = 1:numberOfTimes
        for probabilityIndex = 1:numberOfProbabilities
            minimumWD(beamIndex, probabilityIndex, timeIndex) = ...
                firstMonotoneCrossing(series(beamIndex).WD_mrad, ...
                series(beamIndex).PD(timeIndex, :), ...
                targetProbabilities(probabilityIndex));
        end
    end
end

pointMinimum = reshape(minimumWD(1, :, :), ...
    numberOfProbabilities, numberOfTimes);
lineMinimum = reshape(minimumWD(2, :, :), ...
    numberOfProbabilities, numberOfTimes);
ringMinimum = reshape(minimumWD(3, :, :), ...
    numberOfProbabilities, numberOfTimes);

kappa = ringMinimum ./ lineMinimum;
commonPerformanceRegion = isfinite(lineMinimum) & isfinite(ringMinimum);
kappa(~commonPerformanceRegion) = nan;
end

function drawMethodSchematic(methodLayout, eaSeries, timeThresholds_s)
%DRAWMETHODSCHEMATIC Show time selection and minimum-wD inversion.
exampleTime_s = 30;
exampleProbability = 0.88;
timeIndex = find(timeThresholds_s == exampleTime_s, 1);
assert(~isempty(timeIndex), "The schematic example time is unavailable.");

lineSeries = eaSeries(2);
ringSeries = eaSeries(3);
lineMinimum = firstMonotoneCrossing(lineSeries.WD_mrad, ...
    lineSeries.PD(timeIndex, :), exampleProbability);
ringMinimum = firstMonotoneCrossing(ringSeries.WD_mrad, ...
    ringSeries.PD(timeIndex, :), exampleProbability);
kappaExample = ringMinimum / lineMinimum;

lineMinimumCurve = probabilityCurveAtWidth(lineSeries, lineMinimum);
ringMinimumCurve = probabilityCurveAtWidth(ringSeries, ringMinimum);
lineProbabilityVsWD = cummax(lineSeries.PD(timeIndex, :).');
ringProbabilityVsWD = cummax(ringSeries.PD(timeIndex, :).');

lineColor = [0.00, 0.45, 0.74];
ringColor = [0.85, 0.33, 0.10];
referenceColor = [0.28, 0.30, 0.32];

% Left: a compact P_D(t) view selects one heat-map coordinate.
timeAxes = nexttile(methodLayout, 1);
hold(timeAxes, "on");
grid(timeAxes, "on");
box(timeAxes, "on");
xline(timeAxes, exampleTime_s, ":", "Color", referenceColor, ...
    "LineWidth", 1.5);
yline(timeAxes, exampleProbability, ":", "Color", referenceColor, ...
    "LineWidth", 1.5);
lineMinimumHandle = plot(timeAxes, timeThresholds_s, ...
    lineMinimumCurve, "-", ...
    "Color", lineColor, "LineWidth", 2.0);
ringMinimumHandle = plot(timeAxes, timeThresholds_s, ...
    ringMinimumCurve, "-", ...
    "Color", ringColor, "LineWidth", 2.0);
plot(timeAxes, exampleTime_s, exampleProbability, "o", ...
    "Color", lineColor, "MarkerSize", 10, "LineWidth", 1.7);
plot(timeAxes, exampleTime_s, exampleProbability, "s", ...
    "Color", ringColor, "MarkerSize", 7, "LineWidth", 1.7);
text(timeAxes, exampleTime_s - 22.0, exampleProbability + 0.025, ...
    "selected cell", "Interpreter", "tex", ...
    "FontName", "Helvetica", "FontSize", 12.0, ...
    "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
legend(timeAxes, [lineMinimumHandle, ringMinimumHandle], ...
    {"Line", "Annular"}, "Location", "southwest", ...
    "Box", "off", "FontName", "Helvetica", "FontSize", 11);
xlim(timeAxes, [0, 60]);
ylim(timeAxes, [0, 1]);
xticks(timeAxes, [0, 30, 60]);
yticks(timeAxes, [0, 0.5, 1]);
xlabel(timeAxes, "First-warning time (s)");
ylabel(timeAxes, "P_D(t)", "Interpreter", "tex");
set(timeAxes, "FontName", "Helvetica", "FontSize", 12, ...
    "LineWidth", 1.0, "Layer", "top");

% Right: at fixed t*, invert P_D(t*;wD) to obtain the two minimum widths.
widthAxes = nexttile(methodLayout, 2);
hold(widthAxes, "on");
grid(widthAxes, "on");
box(widthAxes, "on");
plot(widthAxes, lineSeries.WD_mrad, lineProbabilityVsWD, "-o", ...
    "Color", lineColor, "LineWidth", 2, "MarkerSize", 2.8, ...
    "MarkerFaceColor", lineColor);
plot(widthAxes, ringSeries.WD_mrad, ringProbabilityVsWD, "-s", ...
    "Color", ringColor, "LineWidth", 2, "MarkerSize", 2.8, ...
    "MarkerFaceColor", ringColor);
yline(widthAxes, exampleProbability, ":", "Color", referenceColor, ...
    "LineWidth", 2);
plot(widthAxes, [lineMinimum, lineMinimum], [0, exampleProbability], ":", ...
    "Color", lineColor, "LineWidth", 1.3);
plot(widthAxes, [ringMinimum, ringMinimum], [0, exampleProbability], ":", ...
    "Color", ringColor, "LineWidth", 1.3);
plot(widthAxes, lineMinimum, exampleProbability, "o", ...
    "Color", lineColor, "MarkerFaceColor", "w", ...
    "MarkerSize", 7.5, "LineWidth", 1.5);
plot(widthAxes, ringMinimum, exampleProbability, "s", ...
    "Color", ringColor, "MarkerFaceColor", "w", ...
    "MarkerSize", 6.5, "LineWidth", 1.5);

% % text(widthAxes, lineMinimum, exampleProbability + 0.035, ...
%     compose("Line %.1f", lineMinimum), "Color", lineColor, ...
%     "FontName", "Helvetica", "FontSize",12, ...
%     "HorizontalAlignment", "center", "VerticalAlignment", "bottom");
% text(widthAxes, ringMinimum, exampleProbability - 0.045, ...
%     compose("Annular %.1f", ringMinimum), "Color", ringColor, ...
%     "FontName", "Helvetica", "FontSize", 12, ...
%     "HorizontalAlignment", "center", "VerticalAlignment", "top");
ratioText = "\kappa=70/100=0.7"; 
text(widthAxes, 0.95, 0.09, ratioText, ...
    "Units", "normalized", "Interpreter", "tex", ...
    "HorizontalAlignment", "right", "VerticalAlignment", "bottom", ...
    "FontName", "Helvetica", "FontSize", 11, ...
    "BackgroundColor", [1, 1, 1], ...
    "EdgeColor", [0, 0, 0], "Margin", 2);

maximumWidth = 1.8 * max(lineMinimum, ringMinimum);
xlim(widthAxes, [0, maximumWidth]);
ylim(widthAxes, [0, 1]);
xticks(widthAxes, [0, round(ringMinimum), round(lineMinimum), ...
    180]);
xtickangle(widthAxes, 0);
yticks(widthAxes, [0, 0.5, 1]);
yticklabels(widthAxes, []);
xlabel(widthAxes, "w_D (mrad)", "Interpreter", "tex");
% ylabel(widthAxes, "P_D(t^*)", "Interpreter", "tex");

set(widthAxes, "FontName", "Helvetica", "FontSize", 12, ...
    "LineWidth", 1.0, "Layer", "top");
end

function probabilityCurve = probabilityCurveAtWidth(series, wD)
%PROBABILITYCURVEATWIDTH Interpolate the monotone P_D values at one wD.
numberOfTimes = size(series.PD, 1);
probabilityCurve = nan(numberOfTimes, 1);
for timeIndex = 1:numberOfTimes
    monotoneProbability = cummax(series.PD(timeIndex, :).');
    probabilityCurve(timeIndex) = interp1(series.WD_mrad(:), ...
        monotoneProbability, wD, "linear", nan);
end
end

function minimumWD = firstMonotoneCrossing(wD, probability, target)
wD = wD(:);
probability = cummax(probability(:));
crossingIndex = find(probability >= target, 1, "first");

if isempty(crossingIndex)
    minimumWD = nan;
elseif crossingIndex == 1
    minimumWD = wD(1);
else
    x = probability(crossingIndex-1:crossingIndex);
    y = wD(crossingIndex-1:crossingIndex);
    if abs(diff(x)) <= eps(max(abs(x)))
        minimumWD = y(2);
    else
        minimumWD = interp1(x, y, target, "linear");
    end
end
end

function plotKappaHeatmap(ax, timeThresholds_s, targetProbabilities, ...
        kappa, panelTitle)
xEdges = cellEdges(timeThresholds_s);
yEdges = 100 .* cellEdges(targetProbabilities);
[xGrid, yGrid] = meshgrid(xEdges, yEdges);

% Use true-color cell data so that the yellow tolerance band is assigned
% exactly to 0.995 <= kappa <= 1.005, independently of colormap binning.
kappaRGB = reshape(ratioColors(kappa(:)), [size(kappa), 3]);
colorData = nan(size(kappa, 1) + 1, size(kappa, 2) + 1, 3);
colorData(1:end-1, 1:end-1, :) = kappaRGB;
colorData(end, 1:end-1, :) = kappaRGB(end, :, :);
colorData(1:end-1, end, :) = kappaRGB(:, end, :);
colorData(end, end, :) = kappaRGB(end, end, :);

alphaData = zeros(size(kappa, 1) + 1, size(kappa, 2) + 1);
alphaData(1:end-1, 1:end-1) = isfinite(kappa);
alphaData(end, 1:end-1) = alphaData(end-1, 1:end-1);
alphaData(1:end-1, end) = alphaData(1:end-1, end-1);
alphaData(end, end) = alphaData(end-1, end-1);

surface(ax, xGrid, yGrid, zeros(size(xGrid)), colorData, ...
    "FaceColor", "flat", "EdgeColor", [0.87, 0.87, 0.87], ...
    "AlphaData", alphaData, "FaceAlpha", "flat", "LineWidth", 0.12);
view(ax, 2);
set(ax, "YDir", "normal", "Color", [0.93, 0.93, 0.93]);
xlim(ax, [0, 60]);
ylim(ax, [0, 100]);
xticks(ax, 0:10:60);
yticks(ax, 0:20:100);
xtickformat(ax, "%g");
ytickformat(ax, "%g");
xlabel(ax, "First-warning time (s)");
ylabel(ax, "P_D(t) (%)", "Interpreter", "tex");
title(ax, panelTitle, "Interpreter", "none", "FontWeight", "normal");
box(ax, "on");
set(ax, "FontName", "Helvetica", "FontSize", 12, "LineWidth", 1.2, ...
    "Layer", "top");
end

function edges = cellEdges(centers)
centers = centers(:).';
if isscalar(centers)
    halfWidth = max(abs(centers(1)) * 0.05, 0.5);
    edges = centers(1) + [-halfWidth, halfWidth];
    return;
end
midpoints = (centers(1:end-1) + centers(2:end)) / 2;
edges = [centers(1) - (midpoints(1) - centers(1)), ...
    midpoints, centers(end) + (centers(end) - midpoints(end))];
end

function colorMap = ratioColormap(numberOfColors)
%RATIOCOLORMAP Jet-like blue/yellow/green map with a narrow yellow band.
numberOfColors = max(9, round(numberOfColors));
if mod(numberOfColors, 2) == 0
    numberOfColors = numberOfColors + 1;
end
ratioValues = linspace(0.5, 1.5, numberOfColors).';
colorMap = ratioColors(ratioValues);
end

function colors = ratioColors(ratioValues)
%RATIOCOLORS Three-way color assignment with a yellow tolerance band.
ratioValues = ratioValues(:);
colors = nan(numel(ratioValues), 3);

yellowLower = 0.995;
yellowUpper = 1.005;

blueRatioAnchors = [0.50, 0.65, 0.78, 0.90, yellowLower];
blueColorAnchors = [ ...
      0,  30, 180; ... % saturated deep royal blue
      0,  85, 255; ... % intense electric blue
      0, 165, 255; ... % vivid azure
      0, 220, 255; ... % brilliant cyan-blue
     70, 245, 255] / 255; % bright saturated cyan at the band boundary

greenRatioAnchors = [yellowUpper, 1.08, 1.18, 1.32, 1.50];
greenColorAnchors = [ ...
     65, 255,  75; ... % brilliant pure green at the band boundary
      0, 235,  70; ... % intense vivid green
      0, 205,  55; ... % saturated grass green
      0, 160,  40; ... % strong medium green
      0, 105,  30] / 255; % deep saturated forest green

belowBand = ratioValues < yellowLower;
yellowBand = ratioValues >= yellowLower & ratioValues <= yellowUpper;
aboveBand = ratioValues > yellowUpper;

colors(belowBand, :) = interp1(blueRatioAnchors, blueColorAnchors, ...
    ratioValues(belowBand), "linear", "extrap");
colors(yellowBand, :) = repmat([255, 240, 0] / 255, ...
    nnz(yellowBand), 1);
colors(aboveBand, :) = interp1(greenRatioAnchors, greenColorAnchors, ...
    ratioValues(aboveBand), "linear", "extrap");
end

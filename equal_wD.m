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

% Fixed asymmetric limits. Kappa=1 is at 80% of the colorbar length;
% each side has its own linear scale.
% See equal_wD_colormap_notes.md for references and interpretation.
finiteKappa = [results.Kappa];
finiteKappa = finiteKappa(isfinite(finiteKappa));
assert(~isempty(finiteKappa), "No finite kappa values to plot.");
kappaLimits = [0.45, 1.05];

figureHandle = figure("Color", "w", "Position", [80, 80, 850, 400]);
layout = tiledlayout(figureHandle, 1, 3, ...
    "TileSpacing", "compact", "Padding", "compact");

heatmapAxes = gobjects(1, numberOfScenarios);

for scenarioIndex = 1:numberOfScenarios
    heatmapAxes(scenarioIndex) = nexttile(layout, scenarioIndex);
    plotKappaHeatmap(heatmapAxes(scenarioIndex), ...
        timeThresholds_s, targetProbabilities, ...
        results(scenarioIndex).Kappa, scenarioInfo(scenarioIndex).title, ...
        kappaLimits);
end
% The center axes places the shared label above the south colorbar.
xlabel(heatmapAxes(2), "First-warning time (s)", "FontSize", 14);
ylabel(heatmapAxes(1), "Detection probability (%)", ...
    "Interpreter", "tex", "FontSize", 14);
for scenarioIndex = 2:numberOfScenarios
    yticklabels(heatmapAxes(scenarioIndex), []);
end

colormap(figureHandle, ratioColormap(1001, kappaLimits));
set(heatmapAxes, "CLim", [0, 1]);

% One shared colorbar below the entire 1-by-3 layout.
colorbarHandle = colorbar(heatmapAxes(1));
colorbarHandle.Layout.Tile = "south";
colorbarHandle.Limits = [0, 1];
% Positions are transformed, labels are the original physical ratios.
lowerTicks = [0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1];
upperTicks = [1, 1.025, 1.05];
kappaTicks = [lowerTicks, upperTicks(2:end)];
colorbarHandle.Ticks = ratioPosition(kappaTicks, kappaLimits);
colorbarHandle.TickLabels = compose("%.4g", kappaTicks);
colorbarHandle.Label.String = ...
    "\kappa = {\itw}_{D, Min}^{Annular} / w_{D, Min}^{Line}";
colorbarHandle.Label.Interpreter = "tex";
colorbarHandle.Label.FontSize = 14;
colorbarHandle.FontName = "Helvetica";
colorbarHandle.FontSize = 14;


outputBase = fullfile(scriptDir, "equal_wD");
outputPng = outputBase + ".png";
outputFig = outputBase + ".fig";
outputCsv = outputBase + ".csv";
outputMat = outputBase + ".mat";

writetable(allRows, outputCsv);
save(outputMat, "results", "scenarioInfo", "timeThresholds_s", ...
    "targetProbabilities", "smallWdMrad", "eaColumnIndex", "kappaLimits");
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
        kappa, panelTitle, kappaLimits)
xEdges = cellEdges(timeThresholds_s);
yEdges = 100 .* cellEdges(targetProbabilities);
[xGrid, yGrid] = meshgrid(xEdges, yEdges);

% True color reserves the compressed interval's midpoint for EXACTLY 1.
kappaRGB = reshape(ratioColors(kappa(:), kappaLimits), [size(kappa), 3]);
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
    "FaceColor", "flat", "EdgeColor", "none", ...
    "AlphaData", alphaData, "FaceAlpha", "flat");
view(ax, 2);
set(ax, "YDir", "normal", "Color", [0.93, 0.93, 0.93]);
xlim(ax, [0, 60]);
ylim(ax, [0, 100]);
xticks(ax, 0:15:60);
yticks(ax, 0:25:100);
xtickformat(ax, "%g");
ytickformat(ax, "%g");
xtickangle(ax, 0);
title(ax, panelTitle, "Interpreter", "none", "FontWeight", "normal", ...
    "FontSize", 14);
box(ax, "on");
set(ax, "FontName", "Helvetica", "FontSize", 14, "LineWidth", 1.2, ...
    "TitleFontSizeMultiplier", 1, "LabelFontSizeMultiplier", 1, ...
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

function colorMap = ratioColormap(numberOfColors, kappaLimits)
%RATIOCOLORMAP Sample the same transfer function used by the heat maps.
numberOfColors = max(9, round(numberOfColors));
if mod(numberOfColors, 2) == 0
    numberOfColors = numberOfColors + 1;
end
centerPosition = ratioPosition(1, kappaLimits);
positions = linspace(0, 1, numberOfColors).';
below = positions <= centerPosition;
ratioValues = zeros(size(positions));
ratioValues(below) = kappaLimits(1) + positions(below) ./ ...
    centerPosition .* (1 - kappaLimits(1));
ratioValues(~below) = 1 + (positions(~below) - centerPosition) ./ ...
    (1 - centerPosition) .* (kappaLimits(2) - 1);
% Explicit equality avoids floating-point reconstruction at the divider.
ratioValues(abs(positions - centerPosition) < eps) = 1;
colorMap = ratioColors(ratioValues, kappaLimits);
end

function positions = ratioPosition(ratioValues, kappaLimits)
%RATIOPOSITION Independently linear normalization on either side of 1.
positions = nan(size(ratioValues));
centerPosition = 0.8;
below = isfinite(ratioValues) & ratioValues <= 1;
above = isfinite(ratioValues) & ratioValues > 1;
positions(below) = centerPosition .* (ratioValues(below) - kappaLimits(1)) ./ ...
    (1 - kappaLimits(1));
positions(above) = centerPosition + (1 - centerPosition) .* (ratioValues(above) - 1) ./ ...
    (kappaLimits(2) - 1);
end

function colors = ratioColors(ratioValues, kappaLimits)
%RATIOCOLORS Muted jet with a short asymmetric transition compressed to 1.
% Match WD_T_P_2D.m's saturation and brightness adjustments exactly.
colorHSV = rgb2hsv(jet(1024));
colorHSV(:, 2) = 0.55 * colorHSV(:, 2);
colorHSV(:, 3) = 0.35 + 0.55 * colorHSV(:, 3);
baseColorMap = hsv2rgb(colorHSV);
% Normalized source-colormap positions, not physical wD or kappa values.
% Retain [0, 0.60) below 1 (blue/cyan/green) and (0.70, 1] above 1
% (orange/red): the lower side receives twice as much source-colormap range.
% Compress only 10% of the source map; use its yellow midpoint at exact 1,
% without assigning a finite kappa tolerance band.
compressedInterval = [0.4, 0.65];
ratioValues = ratioValues(:);
colors = nan(numel(ratioValues), 3);
positions = ratioPosition(ratioValues, kappaLimits);
centerPosition = ratioPosition(1, kappaLimits);
below = isfinite(ratioValues) & ratioValues < 1;
above = isfinite(ratioValues) & ratioValues > 1;
equal = ratioValues == 1;
sourcePositions = nan(size(ratioValues));
sourcePositions(below) = compressedInterval(1) .* ...
    min(max(positions(below) ./ centerPosition, 0), 1);
sourcePositions(above) = compressedInterval(2) + ...
    (1 - compressedInterval(2)) .* ...
    min(max((positions(above) - centerPosition) ./ ...
    (1 - centerPosition), 0), 1);
sourcePositions(equal) = mean(compressedInterval);
valid = isfinite(sourcePositions);
knots = linspace(0, 1, size(baseColorMap, 1));
colors(valid, :) = interp1(knots, baseColorMap, ...
    sourcePositions(valid), "linear");
end

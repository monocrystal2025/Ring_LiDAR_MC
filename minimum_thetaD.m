clc;
clearvars;
close all;

% Combined 4-by-3 minimum-wD and equal-performance comparison.
% Columns: ROE, ROI-RA, ROI-SP. Rows 1-3: point, line, annular minimum wD.
% Each heat-map cell represents a common performance requirement (t, P_D(t)).
% The fourth-row displayed quantity is
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
minimumWDMaps = cell(3, numberOfScenarios);
fasanDListMrad = 5:5:400;

for scenarioIndex = 1:numberOfScenarios
    info = scenarioInfo(scenarioIndex);
    fprintf("Loading %s results from %s ...\n", info.title, info.resultDir);

    series = loadScenarioSeries(info, smallWdMrad, timeThresholds_s);
    [kappa, pointMinimum, lineMinimum, ringMinimum] = ...
        equalPerformanceMinima(series, targetProbabilities);

    % Preserve WD_T_P_2D's width grid and percentage-based interpolation.
    for beamIndex = 1:3
        [hasWidth, widthIndices] = ismember(fasanDListMrad, ...
            series(beamIndex).WD_mrad);
        assert(all(hasWidth), "Missing widths in the 5:5:400 mrad grid.");
        probabilityByWidth = 100 .* series(beamIndex).PD(:, widthIndices);
        minimumMap = nan(numel(targetProbabilities), numel(timeThresholds_s));
        for timeIndex = 1:numel(timeThresholds_s)
            for probabilityIndex = 1:numel(targetProbabilities)
                minimumMap(probabilityIndex, timeIndex) = ...
                    firstMonotoneCrossing(fasanDListMrad, ...
                    probabilityByWidth(timeIndex, :), ...
                    2 * (probabilityIndex - 1));
            end
        end
        minimumWDMaps{beamIndex, scenarioIndex} = minimumMap;
    end
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

% One 4-by-3 grid; reserve a right-hand strip for the two group colorbars.
fontSize = 13;
fontName = "Helvetica";
figureHandle = figure("Color", "w", "Position", [50, 50, 1200, 1000]);
layout = tiledlayout(figureHandle, 4, 3, ...
    "TileSpacing", "compact", "Padding", "compact");
layout.Units = "normalized";
layout.OuterPosition = [0.015, 0.01, 0.85, 0.98];
heatmapAxes = gobjects(4, 3);

colorHSV = rgb2hsv(jet(1024));
colorHSV(:, 2) = 0.55 * colorHSV(:, 2);
colorHSV(:, 3) = 0.35 + 0.55 * colorHSV(:, 3);
widthColorMap = hsv2rgb(colorHSV);
kappaColorMap = ratioColormap(1001, kappaLimits);

for rowIndex = 1:4
    for scenarioIndex = 1:numberOfScenarios
        ax = nexttile(layout, (rowIndex - 1) * 3 + scenarioIndex);
        heatmapAxes(rowIndex, scenarioIndex) = ax;
        if rowIndex <= 3
            minimumMap = minimumWDMaps{rowIndex, scenarioIndex};
            imageHandle = imagesc(ax, timeThresholds_s, ...
                100 .* targetProbabilities, minimumMap);
            imageHandle.AlphaData = isfinite(minimumMap);
            axis(ax, "xy");
            set(ax, "ColorScale", "log");
            clim(ax, [1, 400]);
            colormap(ax, widthColorMap);
        else
            plotKappaHeatmap(ax, timeThresholds_s, targetProbabilities, ...
                results(scenarioIndex).Kappa, "", kappaLimits);
            colormap(ax, kappaColorMap);
            clim(ax, [0, 1]);
        end
        set(ax, "FontName", fontName, "FontSize", fontSize, ...
            "TitleFontSizeMultiplier", 1, "LabelFontSizeMultiplier", 1, ...
            "LineWidth", 1.2, "Layer", "top");
        xlim(ax, [0, 60]);
        ylim(ax, [0, 100]);
        xticks(ax, 0:10:60);
        yticks(ax, 0:25:100);
        box(ax, "on");
        title(ax, "");
        if rowIndex == 1
            title(ax, scenarioInfo(scenarioIndex).title, ...
                "FontWeight", "normal");
        end
        xlabel(ax, "");
        if rowIndex < 4
            ax.XTickLabel = [];
        end
        ylabel(ax, "");
        if scenarioIndex > 1
            ax.YTickLabel = [];
        end
    end
end

xlabel(layout, "First-warning time (s)", ...
    "FontName", fontName, "FontSize", fontSize);
ylabel(layout, "Detection probability (%)", ...
    "FontName", fontName, "FontSize", fontSize);

% Explicit positions keep both bars outside the same three subplot columns.
widthColorbar = colorbar(heatmapAxes(1, 3));
widthColorbar.Units = "normalized";
widthColorbar.Position = [0.88, 0.35, 0.015, 0.60];
widthColorbar.Ticks = [5, 15, 50, 150, 400];
widthColorbar.TickLabels = {'5', '15', '50', '150', '400'};
widthColorbar.Label.String = 'Minimum required {\itw}_D (mrad)';

kappaColorbar = colorbar(heatmapAxes(4, 3));
kappaColorbar.Units = "normalized";
kappaColorbar.Position = [0.88, 0.08, 0.015, 0.18];
% Fewer labels keep the single-row vertical bar readable at the same font size.
kappaTicks = [0.45, 0.65, 0.85, 1, 1.05];
kappaColorbar.Ticks = ratioPosition(kappaTicks, kappaLimits);
kappaColorbar.TickLabels = compose("%.4g", kappaTicks);
kappaColorbar.Label.String = '\kappa (Annular / Line)';

set([widthColorbar, kappaColorbar], ...
    "FontName", fontName, "FontSize", fontSize);
set([widthColorbar.Label, kappaColorbar.Label], ...
    "Interpreter", "tex", "FontName", fontName, "FontSize", fontSize);
set(findall(figureHandle, "Type", "text"), ...
    "FontName", fontName, "FontSize", fontSize);

drawnow;
% Use pixels in the common layout parent to align colorbar ends exactly.
topPosition = getpixelposition(heatmapAxes(1, 3));
thirdPosition = getpixelposition(heatmapAxes(3, 3));
bottomPosition = getpixelposition(heatmapAxes(4, 3));
set([widthColorbar, kappaColorbar], "Units", "pixels");
barPosition = widthColorbar.Position;
widthColorbar.Position = [barPosition(1), thirdPosition(2), ...
    barPosition(3), ...
    topPosition(2) + topPosition(4) - thirdPosition(2)];
kappaColorbar.Position = [barPosition(1), bottomPosition(2), ...
    barPosition(3), bottomPosition(4)];
set([widthColorbar, kappaColorbar], "Units", "normalized");

outputBase = fullfile(scriptDir, "minimum_thetaD");
outputPng = outputBase + ".png";
outputFig = outputBase + ".fig";
outputCsv = outputBase + ".csv";
outputMat = outputBase + ".mat";
writetable(allRows, outputCsv);
save(outputMat, "results", "scenarioInfo", "timeThresholds_s", ...
    "targetProbabilities", "smallWdMrad", "eaColumnIndex", "kappaLimits", ...
    "minimumWDMaps", "fasanDListMrad");
exportgraphics(figureHandle, outputPng, "Resolution", 300);
savefig(figureHandle, outputFig);
fprintf("Saved PNG: %s\nSaved FIG: %s\nSaved CSV: %s\nSaved MAT: %s\n", ...
    outputPng, outputFig, outputCsv, outputMat);
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
        error("minimum_thetaD:NoFiles", ...
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
        error("minimum_thetaD:DuplicateWD", ...
            "Duplicate %s wD values found for %s.", beamType, info.title);
    end

    probability = nan(numel(timeThresholds_s), numel(files));
    for fileIndex = 1:numel(files)
        resultFile = fullfile(files(fileIndex).folder, files(fileIndex).name);
        data = load(resultFile, detectVariable, firstTimeVariable);
        if ~isfield(data, detectVariable) || ...
                ~isfield(data, firstTimeVariable)
            error("minimum_thetaD:MissingVariable", ...
                "Expected variables are missing in %s.", resultFile);
        end

        detect = data.(detectVariable);
        firstTime = data.(firstTimeVariable);
        if size(detect, 2) < info.columnIndex || ...
                size(firstTime, 2) < info.columnIndex
            error("minimum_thetaD:MissingColumn", ...
                "Column %d is unavailable in %s.", ...
                info.columnIndex, resultFile);
        end
        detect = detect(:, info.columnIndex);
        firstTime = firstTime(:, info.columnIndex);
        if numel(detect) ~= numel(firstTime)
            error("minimum_thetaD:SizeMismatch", ...
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

% True color reserves warm gray for EXACTLY 1, without a finite neutral band.
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
xticks(ax, 0:10:60);
yticks(ax, 0:20:100);
xtickformat(ax, "%g");
ytickformat(ax, "%g");
xlabel(ax, "First-warning time (s)", "FontSize", 14);
ylabel(ax, "Detection probability (%)", "Interpreter", "tex", ...
    "FontSize", 14);
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
%RATIOCOLORS Custom blue / warm gray / orange-red diverging transfer function.
% Each branch has a broad lightness range. This is a custom RGB palette,
% not a claim of perceptual uniformity or a reproduction of a paper's map.
ratioValues = ratioValues(:);
colors = nan(numel(ratioValues), 3);
positions = ratioPosition(ratioValues, kappaLimits);
centerPosition = ratioPosition(1, kappaLimits);
below = isfinite(ratioValues) & ratioValues < 1;
above = isfinite(ratioValues) & ratioValues > 1;
equal = ratioValues == 1;
% Warm gray is distinct from both branches and the gray unavailable region.
blue = [8 48 107; 8 81 156; 33 113 181; 66 146 198; ...
    107 174 214; 158 202 225; 198 219 239] ./ 255;
warm = [254 230 206; 253 190 133; 253 141 60; 241 105 19; ...
    217 72 1; 166 54 3; 127 39 4] ./ 255;
knots = linspace(0, 1, size(blue, 1));
colors(below, :) = interp1(knots, blue, ...
    min(max(positions(below) ./ centerPosition, 0), 1), "linear");
colors(above, :) = interp1(knots, warm, ...
    min(max((positions(above) - centerPosition) ./ ...
    (1 - centerPosition), 0), 1), "linear");
colors(equal, :) = repmat([184 177 164] ./ 255, nnz(equal), 1);
end

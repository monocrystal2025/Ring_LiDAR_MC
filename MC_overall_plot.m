clc;
clearvars;
close all;

% Plot MC, ROI_RA, and ROI_SP results in one 2-by-3 summary figure.
scriptDir = fileparts(mfilename("fullpath"));
if strlength(scriptDir) == 0
    scriptDir = pwd;
end

dataRoot = fileparts(fileparts(scriptDir));
mcResultDir = fullfile(dataRoot, "MC_RESULTS");
roiResultDir = fullfile(scriptDir, "MC_ROI_RESULTS_SNR_NEW");

if ~isfolder(mcResultDir)
    error("MC result folder not found: %s", mcResultDir);
end
if ~isfolder(roiResultDir)
    error("ROI result folder not found: %s", roiResultDir);
end

fasanDMrad = 5:5:400;
fasanSmallDMrad = 1;
mcColumnIndex = 1;
beamTypes = ["POINT", "LINE", "RING"];
beamLabels = ["spot", "line", "annular"];
trackTypes = ["RA", "SP"];
columnTitles = ["MC", "ROI_RA", "ROI_SP"];
markers = ["^", "s", "o"];
colors = [ ...
    0.4660, 0.6740, 0.1880; ...
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980];

[mcProbability, mcFirstWarningTime] = loadMCResults( ...
    mcResultDir, fasanDMrad, fasanSmallDMrad, mcColumnIndex, beamTypes);
[roiProbability, roiFirstWarningTime] = loadROIResults( ...
    roiResultDir, fasanDMrad, fasanSmallDMrad, beamTypes, trackTypes);

% Preserve the point-beam interpolation used by plot_MC_EA_fasanD.m.
mcProbability(1, :) = fillMissingByInterpolation( ...
    fasanDMrad, mcProbability(1, :), "point probability");
mcFirstWarningTime(1, :) = fillMissingByInterpolation( ...
    fasanDMrad, mcFirstWarningTime(1, :), "point first warning time");

probability = cat(2, reshape(mcProbability, 3, 1, []), roiProbability);
firstWarningTime = cat(2, ...
    reshape(mcFirstWarningTime, 3, 1, []), roiFirstWarningTime);

figureHandle = figure("Color", "w", "Position", [100, 100, 1200, 650]);
layout = tiledlayout(figureHandle, 2, 3, ...
    "TileSpacing", "compact", "Padding", "compact");

numberOfBeams = numel(beamTypes);
numberOfColumns = numel(columnTitles);
probabilityAxes = gobjects(1, numberOfColumns);
timeAxes = gobjects(1, numberOfColumns);
legendLines = gobjects(1, numberOfBeams);

for columnIndex = 1:numberOfColumns
    probabilityAxes(columnIndex) = nexttile(layout, columnIndex);
    hold(probabilityAxes(columnIndex), "on");
    for beamIndex = 1:numberOfBeams
        yData = 100*reshape( ...
            probability(beamIndex, columnIndex, :), 1, []);
        lineHandle = plot(probabilityAxes(columnIndex), fasanDMrad, yData, ...
            "-" + markers(beamIndex), "LineWidth", 1.8, ...
            "MarkerSize", 4, "Color", colors(beamIndex, :), ...
            "DisplayName", beamLabels(beamIndex));
        if columnIndex == 1
            legendLines(beamIndex) = lineHandle;
        end
    end
    formatProbabilityAxes(probabilityAxes(columnIndex), ...
        columnTitles(columnIndex), columnIndex == 1);

    timeAxes(columnIndex) = nexttile(layout, numberOfColumns + columnIndex);
    hold(timeAxes(columnIndex), "on");
    for beamIndex = 1:numberOfBeams
        yData = reshape( ...
            firstWarningTime(beamIndex, columnIndex, :), 1, []);
        plot(timeAxes(columnIndex), fasanDMrad, yData, ...
            "-" + markers(beamIndex), "LineWidth", 1.8, ...
            "MarkerSize", 4, "Color", colors(beamIndex, :), ...
            "DisplayName", beamLabels(beamIndex));
    end
    formatTimeAxes(timeAxes(columnIndex), columnIndex == 1);
end

linkaxes(probabilityAxes, "xy");
linkaxes(timeAxes, "xy");

legendHandle = legend(probabilityAxes(1), legendLines, beamLabels, ...
    "Location", "best", "FontName", "Helvetica", "FontSize", 13);
legendHandle.Layout.Tile = "north";
legendHandle.NumColumns = numberOfBeams;
legendHandle.Box = "off";
legendHandle.ItemTokenSize = [50, 18];

outputPng = fullfile(scriptDir, "MC_overall_plot.png");
outputFig = fullfile(scriptDir, "MC_overall_plot.fig");
exportgraphics(figureHandle, outputPng, "Resolution", 300);
savefig(figureHandle, outputFig);

fprintf("Saved PNG: %s\n", outputPng);
fprintf("Saved FIG: %s\n", outputFig);

function [probability, firstWarningTime] = loadMCResults( ...
        resultDir, fasanDMrad, smallDMrad, columnIndex, beamTypes)
    numberOfBeams = numel(beamTypes);
    numberOfDiameters = numel(fasanDMrad);
    probability = nan(numberOfBeams, numberOfDiameters);
    firstWarningTime = nan(numberOfBeams, numberOfDiameters);

    for beamIndex = 1:numberOfBeams
        beamType = beamTypes(beamIndex);
        [detectVariable, timeVariable] = mcVariableNames(beamType);
        for diameterIndex = 1:numberOfDiameters
            resultFile = mcResultFileName(resultDir, beamType, ...
                fasanDMrad(diameterIndex), smallDMrad);
            if ~isfile(resultFile)
                warning("Missing MC result: %s", resultFile);
                continue;
            end

            result = load(resultFile, detectVariable, timeVariable);
            if ~isfield(result, detectVariable) || ...
                    ~isfield(result, timeVariable)
                warning("Missing expected variables in: %s", resultFile);
                continue;
            end
            if size(result.(detectVariable), 2) < columnIndex || ...
                    size(result.(timeVariable), 2) < columnIndex
                warning("Column %d is unavailable in: %s", ...
                    columnIndex, resultFile);
                continue;
            end

            [probability(beamIndex, diameterIndex), ...
                firstWarningTime(beamIndex, diameterIndex)] = ...
                summarizeMCResult( ...
                result.(detectVariable)(:, columnIndex), ...
                result.(timeVariable)(:, columnIndex));
        end
    end
end

function [probability, firstWarningTime] = loadROIResults( ...
        resultDir, fasanDMrad, smallDMrad, beamTypes, trackTypes)
    numberOfBeams = numel(beamTypes);
    numberOfTracks = numel(trackTypes);
    numberOfDiameters = numel(fasanDMrad);
    probability = nan(numberOfBeams, numberOfTracks, numberOfDiameters);
    firstWarningTime = nan(numberOfBeams, numberOfTracks, numberOfDiameters);

    for trackIndex = 1:numberOfTracks
        trackType = trackTypes(trackIndex);
        for beamIndex = 1:numberOfBeams
            beamType = beamTypes(beamIndex);
            [detectVariable, timeVariable] = roiVariableNames( ...
                beamType, trackType);
            for diameterIndex = 1:numberOfDiameters
                resultFile = roiResultFileName(resultDir, beamType, ...
                    trackType, fasanDMrad(diameterIndex), smallDMrad);
                if ~isfile(resultFile)
                    warning("Missing ROI result: %s", resultFile);
                    continue;
                end

                result = load(resultFile, detectVariable, timeVariable);
                if ~isfield(result, detectVariable) || ...
                        ~isfield(result, timeVariable)
                    warning("Missing expected variables in: %s", resultFile);
                    continue;
                end

                [probability(beamIndex, trackIndex, diameterIndex), ...
                    firstWarningTime(beamIndex, trackIndex, diameterIndex)] = ...
                    summarizeMCResult(result.(detectVariable), ...
                    result.(timeVariable));
            end
        end
    end
end

function resultFile = mcResultFileName( ...
        resultDir, beamType, diameterMrad, smallDMrad)
    diameterText = num2str(diameterMrad, "%g");
    if beamType == "POINT"
        fileName = "MC_1par_EA_POINT_D" + diameterText + "mrad.mat";
    else
        smallDiameterText = num2str(smallDMrad, "%g");
        fileName = "MC_1par_EA_" + beamType + "_D" + diameterText + ...
            "d" + smallDiameterText + "mrad.mat";
    end
    resultFile = fullfile(resultDir, fileName);
end

function resultFile = roiResultFileName( ...
        resultDir, beamType, trackType, diameterMrad, smallDMrad)
    diameterText = num2str(diameterMrad, "%g");
    if beamType == "POINT"
        fileName = "MC_1par_ROI_NEW_POINT_" + trackType + ...
            "_D" + diameterText + "mrad.mat";
    else
        smallDiameterText = num2str(smallDMrad, "%g");
        fileName = "MC_1par_ROI_NEW_" + beamType + "_" + trackType + ...
            "_D" + diameterText + "d" + smallDiameterText + "mrad.mat";
    end
    resultFile = fullfile(resultDir, fileName);
end

function [detectVariable, timeVariable] = mcVariableNames(beamType)
    beamCode = beamCodeForType(beamType);
    detectVariable = "detect_" + beamCode;
    timeVariable = "first_time_" + beamCode;
end

function [detectVariable, timeVariable] = roiVariableNames( ...
        beamType, trackType)
    beamCode = beamCodeForType(beamType);
    detectVariable = "detect_" + beamCode + "_" + trackType;
    timeVariable = "first_time_" + beamCode + "_" + trackType;
end

function beamCode = beamCodeForType(beamType)
    switch beamType
        case "POINT"
            beamCode = "P";
        case "LINE"
            beamCode = "L";
        case "RING"
            beamCode = "R";
        otherwise
            error("Unsupported beam type: %s", beamType);
    end
end

function [probability, firstWarningTime] = summarizeMCResult( ...
        detect, firstTime)
    detect = detect(:);
    firstTime = firstTime(:);
    if numel(detect) ~= numel(firstTime)
        error("Detection and first-warning-time arrays have different sizes.");
    end

    isDetected = detect ~= 0;
    probability = mean(isDetected);
    hasValidTime = isDetected & isfinite(firstTime);
    if any(hasValidTime)
        firstWarningTime = mean(firstTime(hasValidTime));
    else
        firstWarningTime = nan;
    end
end

function yFilled = fillMissingByInterpolation(x, yRaw, dataName)
    originalSize = size(yRaw);
    x = x(:);
    yRaw = yRaw(:);
    yFilled = yRaw;
    isValid = isfinite(x) & isfinite(yRaw);

    if ~any(isValid)
        warning("No valid %s data found. The curve will remain NaN.", dataName);
        yFilled = reshape(yFilled, originalSize);
        return;
    end

    if nnz(isValid) == 1
        yFilled(:) = yRaw(isValid);
        warning("Only one valid %s data point was found; missing values were filled with it.", ...
            dataName);
        yFilled = reshape(yFilled, originalSize);
        return;
    end

    isMissing = ~isfinite(yFilled);
    if any(isMissing)
        yFilled(isMissing) = interp1(x(isValid), yRaw(isValid), ...
            x(isMissing), "pchip", "extrap");
        fprintf("Filled %d missing %s points by pchip interpolation.\n", ...
            nnz(isMissing), dataName);
    end
    yFilled = reshape(yFilled, originalSize);
end

function formatProbabilityAxes(axesHandle, panelTitle, showYLabels)
    grid(axesHandle, "on");
    box(axesHandle, "on");
    title(axesHandle, panelTitle, "Interpreter", "none", ...
        "FontName", "Helvetica", "FontSize", 13, "FontWeight", "normal");
    xlabel(axesHandle, "");
    xticks(axesHandle, 0:100:500);
    yticks(axesHandle, 0:25:100);
    ylim(axesHandle, [0, 100]);
    if showYLabels
        ylabel(axesHandle, "{\it P}_{\rmdetect} (%)", "Interpreter", "tex");
    else
        ylabel(axesHandle, "");
    end
    set(axesHandle, "FontName", "Helvetica", "FontSize", 13, ...
        "LineWidth", 2, "XScale", "linear");
    set([axesHandle.XLabel, axesHandle.YLabel], ...
        "FontName", "Helvetica", "FontSize", 13);
end

function formatTimeAxes(axesHandle, showYLabels)
    grid(axesHandle, "on");
    box(axesHandle, "on");
    xlabel(axesHandle, "{\it w}_{\rmD} (mrad)", "Interpreter", "tex");
    xticks(axesHandle, 0:100:500);
    if showYLabels
        ylabel(axesHandle, "{\it t}_{\rm warn} (s)", ...
            "Interpreter", "tex");
    else
        ylabel(axesHandle, "");
    end
    set(axesHandle, "FontName", "Helvetica", "FontSize", 13, ...
        "LineWidth", 2, "XScale", "linear");
    set([axesHandle.XLabel, axesHandle.YLabel], ...
        "FontName", "Helvetica", "FontSize", 13);
end

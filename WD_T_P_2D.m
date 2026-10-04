clc;
clear;
close all;

% Two-dimensional maps of minimum wD for equal performance requirements.
% Columns: ROE, ROI-RA, and ROI-SP.
% Rows: point, line, and annular illumination fields.
% The horizontal coordinate is the first-warning time limit, the vertical
% coordinate is the target detection probability, and color is the minimum
% wD that reaches the corresponding performance requirement.
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

EAResultDir = 'D:\lzx\MC_RESULTS';
ROIResultDir = fullfile(scriptDir, 'MC_ROI_RESULTS_SNR_NEW');

fasanDListMrad = 5:5:400;
fasanDMrad = 1;
columnIndex = 1;

timeStep = 1;
timeLimit = 60;
timeGrid = 0:timeStep:timeLimit;
targetProbabilityGrid = 0:2:100;

beamTypes = {'POINT', 'LINE', 'RING'};
columnLabels = {'ROE', 'ROI-RA', 'ROI-SP'};
pathTypes = {'', 'RA', 'SP'};

numberOfRows = numel(beamTypes);
numberOfColumns = numel(columnLabels);
minimumWDMaps = cell(numberOfRows, numberOfColumns);

for rowIndex = 1:numberOfRows
    beamType = beamTypes{rowIndex};

    for columnIndexPlot = 1:numberOfColumns
        pathType = pathTypes{columnIndexPlot};
        probabilityByWidth = zeros(numel(timeGrid), ...
            numel(fasanDListMrad));

        for widthIndex = 1:numel(fasanDListMrad)
            fasanDText = num2str(fasanDListMrad(widthIndex), '%g');
            fasanDTextSmall = num2str(fasanDMrad, '%g');

            if columnIndexPlot == 1
                resultFile = EAResultFile(EAResultDir, beamType, ...
                    fasanDText, fasanDTextSmall);
                detectVariable = ['detect_', beamType(1)];
                firstTimeVariable = ['first_time_', beamType(1)];
            else
                resultFile = ROIResultFile(ROIResultDir, beamType, ...
                    pathType, fasanDText, fasanDTextSmall);
                detectVariable = ['detect_', beamType(1), '_', pathType];
                firstTimeVariable = ...
                    ['first_time_', beamType(1), '_', pathType];
            end

            probabilityByWidth(:, widthIndex) = 100 * ...
                loadTimeLimitedProbability(resultFile, detectVariable, ...
                firstTimeVariable, columnIndex, timeGrid).';
        end

        minimumWDMap = nan(numel(targetProbabilityGrid), numel(timeGrid));
        for timeIndex = 1:numel(timeGrid)
            for probabilityIndex = 1:numel(targetProbabilityGrid)
                minimumWDMap(probabilityIndex, timeIndex) = ...
                    firstMonotoneCrossing(fasanDListMrad, ...
                    probabilityByWidth(timeIndex, :), ...
                    targetProbabilityGrid(probabilityIndex));
            end
        end

        minimumWDMaps{rowIndex, columnIndexPlot} = minimumWDMap;
    end
end

figureHandle = figure('Color', 'w', 'Position', [50, 50, 650, 450]);
layout = tiledlayout(figureHandle, numberOfRows, numberOfColumns, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for rowIndex = 1:numberOfRows
    for columnIndexPlot = 1:numberOfColumns
        tileIndex = (rowIndex - 1) * numberOfColumns + columnIndexPlot;
        axesHandle = nexttile(layout, tileIndex);
        minimumWDMap = minimumWDMaps{rowIndex, columnIndexPlot};
        imageHandle = imagesc(axesHandle, timeGrid, ...
            targetProbabilityGrid, minimumWDMap);
        imageHandle.AlphaData = isfinite(minimumWDMap);
        % Unattainable requirements (NaN) reveal the black axes background.
        axesHandle.Color = [1, 1, 1];
        axis(axesHandle, 'xy');
        box(axesHandle, 'on');
        set(axesHandle, 'ColorScale', 'log');
        clim(axesHandle, [1, 400]);
        xlim(axesHandle, [timeGrid(1), timeGrid(end)]);
        ylim(axesHandle, ...
            [targetProbabilityGrid(1), targetProbabilityGrid(end)]);
        xticks(axesHandle, 0:15:60);
        xtickangle(axesHandle, 0);
        yticks(axesHandle, 0:25:100);

        if rowIndex == 1
            title(axesHandle, columnLabels{columnIndexPlot}, ...
                'FontName', 'Helvetica', 'FontSize', 13, ...
                'FontWeight', 'normal');
        end

        if rowIndex < numberOfRows
            axesHandle.XTickLabel = [];
        end

        if columnIndexPlot > 1
            axesHandle.YTickLabel = [];
        end

        set(axesHandle, 'FontName', 'Helvetica', 'FontSize', 13, ...
            'TitleFontSizeMultiplier', 1, ...
            'LabelFontSizeMultiplier', 1, ...
            'LineWidth', 1.2, 'Layer', 'top');
    end
end

xlabel(layout, 'First-warning time (s)', ...
    'FontName', 'Helvetica', 'FontSize', 13, 'FontWeight', 'normal');
ylabel(layout, 'Detection probability (%)', ...
    'FontName', 'Helvetica', 'FontSize', 13, 'FontWeight', 'normal');

% Keep jet's continuous hue progression with softer, muted colors
% inspired by equal_wD.png. Lift dark ends and soften peak brightness.
colorHSV = rgb2hsv(jet(1024));
colorHSV(:, 2) = 0.55 * colorHSV(:, 2);
colorHSV(:, 3) = 0.35 + 0.55 * colorHSV(:, 3);
smoothColorMap = hsv2rgb(colorHSV);
colormap(figureHandle, smoothColorMap);
colorbarHandle = colorbar;
colorbarHandle.Layout.Tile = 'east';
colorbarHandle.Label.String = 'Minimum required {\it\theta}_D (mrad)';
colorbarHandle.Label.Interpreter = 'tex';
colorbarHandle.Label.FontName = 'Helvetica';
colorbarHandle.Label.FontSize = 13;
colorbarHandle.FontName = 'Helvetica';
colorbarHandle.FontSize = 13;
% Rounded integer ticks spaced approximately evenly on the logarithmic scale.
colorbarHandle.Ticks = [1, 5, 20, 100, 400];
colorbarHandle.TickLabels = compose('%d', colorbarHandle.Ticks);

exportgraphics(figureHandle, fullfile(scriptDir, 'WD_T-P_2D.png'), ...
    'Resolution', 300);
savefig(figureHandle, fullfile(scriptDir, 'WD_T-P_2D.fig'));

function resultFile = EAResultFile(resultDir, beamType, DText, dText)
if strcmp(beamType, 'POINT')
    fileName = ['MC_1par_EA_POINT_D', DText, 'mrad.mat'];
else
    fileName = ['MC_1par_EA_', beamType, '_D', DText, ...
        'd', dText, 'mrad.mat'];
end
resultFile = fullfile(resultDir, fileName);
end

function resultFile = ROIResultFile(resultDir, beamType, pathType, ...
        DText, dText)
if strcmp(beamType, 'POINT')
    fileName = ['MC_1par_ROI_NEW_POINT_', pathType, ...
        '_D', DText, 'mrad.mat'];
else
    fileName = ['MC_1par_ROI_NEW_', beamType, '_', pathType, ...
        '_D', DText, 'd', dText, 'mrad.mat'];
end
resultFile = fullfile(resultDir, fileName);
end

function probability = loadTimeLimitedProbability(resultFile, ...
        detectVariable, firstTimeVariable, columnIndex, timeGrid)
assert(isfile(resultFile), 'Missing result file: %s', resultFile);
result = load(resultFile, detectVariable, firstTimeVariable);
assert(isfield(result, detectVariable), ...
    'Missing variable %s in %s.', detectVariable, resultFile);
assert(isfield(result, firstTimeVariable), ...
    'Missing variable %s in %s.', firstTimeVariable, resultFile);

detect = result.(detectVariable);
firstTime = result.(firstTimeVariable);
assert(size(detect, 2) >= columnIndex && ...
    size(firstTime, 2) >= columnIndex, ...
    'Column %d is unavailable in %s.', columnIndex, resultFile);

probability = timeLimitedProbability(detect(:, columnIndex), ...
    firstTime(:, columnIndex), timeGrid);
end

function probability = timeLimitedProbability(detect, firstTime, timeGrid)
detect = detect(:);
firstTime = firstTime(:);
assert(numel(detect) == numel(firstTime), ...
    'Detection flags and first-warning times must have equal lengths.');

validDetection = detect ~= 0 & isfinite(firstTime) & firstTime >= 0;
detectionTimes = sort(firstTime(validDetection));
numberOfTrials = numel(detect);
numberOfDetections = numel(detectionTimes);

probability = zeros(size(timeGrid));
detectionIndex = 1;
for timeIndex = 1:numel(timeGrid)
    while detectionIndex <= numberOfDetections && ...
            detectionTimes(detectionIndex) <= timeGrid(timeIndex)
        detectionIndex = detectionIndex + 1;
    end
    probability(timeIndex) = (detectionIndex - 1) / numberOfTrials;
end
end


function minimumWD = firstMonotoneCrossing(wD, probability, target)
wD = wD(:);
probability = cummax(probability(:));
crossingIndex = find(probability >= target, 1, 'first');

if isempty(crossingIndex)
    minimumWD = nan;
elseif crossingIndex == 1
    minimumWD = wD(1);
else
    probabilityPair = probability(crossingIndex-1:crossingIndex);
    wDPair = wD(crossingIndex-1:crossingIndex);
    if abs(diff(probabilityPair)) <= eps(max(abs(probabilityPair)))
        minimumWD = wDPair(2);
    else
        minimumWD = interp1(probabilityPair, wDPair, target, 'linear');
    end
end
end

clc;
clear;
close all;

% Time-limited detection probability for EA and two ROI scan paths.
% Rows correspond to beam widths; columns correspond to EA, ROI-RA, and
% ROI-SP. P_D(t) is the fraction of all Monte Carlo trials that have
% generated their first warning no later than time t.
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

EA_resultDir = 'D:\lzx\MC_RESULTS';
ROI_resultDir = fullfile(scriptDir, 'MC_ROI_RESULTS_SNR_NEW');
fasan_DList = [50, 125, 195] * 1e-3;
fasan_d = 1e-3;
columnIndex = 1;

timeStep = 0.1;
timeLimit = 60;
timeGrid = 0:timeStep:timeLimit;

beamTypes = {'POINT', 'LINE', 'RING'};
beamLabels = {'Spot', 'Line', 'Annular'};
lineStyles = {'--', ':', '-'};
% Match plot_MC_EA_fasanD.m: spot=green, line=blue, annular=orange.
beamColors = [ ...
    0.4660, 0.6740, 0.1880; ...
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980];
columnLabels = {'ROE', 'ROI-RA', 'ROI-SP'};
pathTypes = {'', 'RA', 'SP'};
dText = num2str(fasan_d * 1e3, '%g');

figureHandle = figure('Color', 'w', 'Position', [50, 50, 720, 600]);
layout = tiledlayout(figureHandle, 3, 3, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
sharedLegendHandles = gobjects(numel(beamTypes), 1);

for rowIndex = 1:numel(fasan_DList)
    DText = num2str(fasan_DList(rowIndex) * 1e3, '%g');

    for columnIndexPlot = 1:numel(columnLabels)
        axesHandle = nexttile(layout, ...
            (rowIndex - 1) * numel(columnLabels) + columnIndexPlot);
        hold(axesHandle, 'on');

        for beamIndex = 1:numel(beamTypes)
            beamType = beamTypes{beamIndex};
            pathType = pathTypes{columnIndexPlot};

            if columnIndexPlot == 1
                resultFile = EAResultFile(EA_resultDir, beamType, ...
                    DText, dText);
                detectVariable = ['detect_', beamType(1)];
                firstTimeVariable = ['first_time_', beamType(1)];
            else
                resultFile = ROIResultFile(ROI_resultDir, beamType, ...
                    pathType, DText, dText);
                detectVariable = ['detect_', beamType(1), '_', pathType];
                firstTimeVariable = ...
                    ['first_time_', beamType(1), '_', pathType];
            end

            probability = loadTimeLimitedProbability(resultFile, ...
                detectVariable, firstTimeVariable, columnIndex, timeGrid);
            lineHandle = plot(axesHandle, timeGrid, 100 * probability, ...
                'LineStyle', lineStyles{beamIndex}, ...
                'Color', beamColors(beamIndex, :), ...
                'LineWidth', 3, 'DisplayName', beamLabels{beamIndex});
            if rowIndex == 1 && columnIndexPlot == 1
                sharedLegendHandles(beamIndex) = lineHandle;
            end
        end

        grid(axesHandle, 'on');
        box(axesHandle, 'on');
        xlim(axesHandle, [0, timeLimit]);
        ylim(axesHandle, [0, 100]);
        xticks(axesHandle, [0, 20, 40, 60]);
        yticks(axesHandle, [0, 25, 50, 75, 100]);

        if rowIndex == 1
            title(axesHandle, columnLabels{columnIndexPlot}, ...
                'FontName', 'Helvetica', 'FontSize', 13);
        end
        if rowIndex == numel(fasan_DList)
            xLabelHandle = xlabel(axesHandle, 'First-warning time (s)');
        else
            axesHandle.XTickLabel = [];
            xLabelHandle = gobjects(0);
        end
        if columnIndexPlot == 1
            yLabelHandle = ylabel(axesHandle, ...
                '{\it P}_{D}({\itt}) (%)', 'Interpreter', 'tex');
        else
            axesHandle.YTickLabel = [];
            yLabelHandle = gobjects(0);
        end

        text(axesHandle, 0.99, 0.20, ...
            sprintf('{\\itw}_{\\rmD}=%g mrad', ...
            fasan_DList(rowIndex) * 1e3), ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'bottom', ...
            'Interpreter', 'tex', ...
            'BackgroundColor', 'none', ...
            'EdgeColor', 'none', ...
            'Margin', 3, ...
            'FontName', 'Helvetica', ...
            'FontSize', 12, ...
            'FontWeight', 'bold');
        set(axesHandle, 'FontName', 'Helvetica', ...
            'FontSize', 14, 'LineWidth', 2);
        set([xLabelHandle; yLabelHandle], ...
            'FontName', 'Helvetica', 'FontSize', 14);
    end
end

sharedLegend = legend(sharedLegendHandles, beamLabels, ...
    'Orientation', 'horizontal', ...
    'Box', 'off', ...
    'FontName', 'Helvetica', ...
    'FontSize', 14);
sharedLegend.Layout.Tile = 'north';
sharedLegend.NumColumns = numel(beamLabels);
sharedLegend.ItemTokenSize = [45, 18];

exportgraphics(figureHandle, fullfile(scriptDir, 'P_D-T.png'), ...
    'Resolution', 300);
savefig(figureHandle, fullfile(scriptDir, 'P_D-T.fig'));

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
detectionTimes = firstTime(validDetection);
numberOfTrials = numel(detect);

probability = arrayfun( ...
    @(time) nnz(detectionTimes <= time) / numberOfTrials, timeGrid);
end

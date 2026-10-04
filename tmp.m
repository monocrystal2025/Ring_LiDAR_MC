% ROE, ROI-RA, and ROI-SP: one row showing 1-Pdetect only.
% Read MAT files, open a figure, and export blind.png at 300 dpi.

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end
resultDir = 'D:\matlab\STATIC_MC_RESULTS';
% blind_SNR_ROI.m saves to fullfile(pwd, 'STATIC_MC_ROI_RESULTS_NEW').
% Resolve that workspace folder relative to this plotting script.
roiResultDir = fullfile(scriptDir, 'STATIC_MC_ROI_RESULTS_NEW');
% Match the current blind_snr_EA/blind_SNR_ROI runs, excluding older files.
fasanDMrad = 5:5:400;
fasanSmallDMrad = 1;
mcColumnIndex = 1; % The two stored columns are duplicates, not two trials.
beamTypes = ["POINT", "LINE", "RING"];
beamCodes = ["P", "L", "R"];
beamLabels = ["Spot", "Line", "Annular"];
columnTitles = ["ROE", "ROI-RA", "ROI-SP"];
trackTypes = ["", "RA", "SP"];
resultDirs = {resultDir, roiResultDir, roiResultDir};
colors = [ ...
    0.4660, 0.6740, 0.1880; ...
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980];

numberOfBeams = numel(beamTypes);
numberOfColumns = numel(columnTitles);
probability = nan(numberOfBeams, numberOfColumns, numel(fasanDMrad));
firstWarningTime = nan(size(probability));
loadedCount = zeros(1, numberOfColumns);
for columnIndex = 1:numberOfColumns
    if ~isfolder(resultDirs{columnIndex})
        error('tmp:MissingResultFolder', ...
            '%s result folder not found: %s', ...
            columnTitles(columnIndex), resultDirs{columnIndex});
    end
    [columnProbability, columnTime, loadedCount(columnIndex)] = loadMCResults( ...
        resultDirs{columnIndex}, fasanDMrad, fasanSmallDMrad, mcColumnIndex, ...
        beamTypes, beamCodes, trackTypes(columnIndex));
    if loadedCount(columnIndex) == 0
        error('tmp:NoResults', 'No usable %s results found in: %s', ...
            columnTitles(columnIndex), resultDirs{columnIndex});
    end
    probability(:, columnIndex, :) = reshape( ...
        columnProbability, numberOfBeams, 1, []);
    firstWarningTime(:, columnIndex, :) = reshape( ...
        columnTime, numberOfBeams, 1, []);
end
% Missing results stay NaN: no interpolation/extrapolation of measured data.

figureHandle = figure('Color', 'w', 'Position', [100, 100, 850, 350], ...
    'Name', 'Static blind-search SNR - ROE / ROI', 'NumberTitle', 'off');
layout = tiledlayout(figureHandle, 1, 3, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
probabilityAxes = gobjects(1, numberOfColumns);
legendLines = gobjects(1, numberOfBeams);

for columnIndex = 1:numberOfColumns
    probabilityAxes(columnIndex) = nexttile(layout, columnIndex);
    hold(probabilityAxes(columnIndex), 'on');
    for beamIndex = 1:numberOfBeams
        % Complement only the plotted values; keep probability/data unchanged.
        yProbability = 100 * (1 - reshape( ...
            probability(beamIndex, columnIndex, :), 1, []));
        lineHandle = plot(probabilityAxes(columnIndex), fasanDMrad, ...
            yProbability, '-', 'LineWidth', 3, ...
            'Color', colors(beamIndex, :), ...
            'DisplayName', beamLabels(beamIndex));
        if columnIndex == 1
            legendLines(beamIndex) = lineHandle;
        end
    end

    formatMCAxes(probabilityAxes(columnIndex), fasanDMrad, ...
        columnTitles(columnIndex));
    ylim(probabilityAxes(columnIndex), [0, 100]);
    yticks(probabilityAxes(columnIndex), 0:25:100);
    if columnIndex > 1
        yticklabels(probabilityAxes(columnIndex), {});
    end
end
linkaxes(probabilityAxes, 'x');
xlabel(layout, '{\it\theta}_{\rmD} (mrad)', 'Interpreter', 'tex', ...
    'FontName', 'Helvetica', 'FontSize', 14);
ylabel(layout, 'Blind-zone ratio (%)', 'Interpreter', 'tex', ...
    'FontName', 'Helvetica', 'FontSize', 14);

legendHandle = legend(probabilityAxes(1), legendLines, beamLabels, ...
    'Location', 'best', 'FontName', 'Helvetica', 'FontSize', 14);
legendHandle.Layout.Tile = 'north';
legendHandle.NumColumns = numel(beamTypes);
legendHandle.Box = 'off';
legendHandle.ItemTokenSize = [50, 18];
drawnow;

for columnIndex = 1:numberOfColumns
    fprintf('Loaded %d/%d %s result files from %s.\n', ...
        loadedCount(columnIndex), numberOfBeams * numel(fasanDMrad), ...
        columnTitles(columnIndex), resultDirs{columnIndex});
end
fprintf('One row: 1 - detection probability for ROE, ROI-RA, and ROI-SP.\n');

outputPng = fullfile(scriptDir, 'blind.png');
exportgraphics(figureHandle, outputPng, 'Resolution', 300);
fprintf('Saved 300 dpi PNG: %s\n', outputPng);

function [probability, firstWarningTime, loadedCount] = loadMCResults( ...
        resultDir, fasanDMrad, smallDMrad, columnIndex, ...
        beamTypes, beamCodes, trackType)
    probability = nan(numel(beamTypes), numel(fasanDMrad));
    firstWarningTime = nan(size(probability));
    loadedCount = 0;

    for beamIndex = 1:numel(beamTypes)
        detectVariable = "detect_" + beamCodes(beamIndex);
        timeVariable = "first_time_" + beamCodes(beamIndex);
        if strlength(trackType) == 0
            fileStem = "STATIC_MC_1par_EA_" + beamTypes(beamIndex);
        else
            fileStem = "STATIC_MC_1par_ROI_" + beamTypes(beamIndex) + ...
                "_" + trackType;
            detectVariable = detectVariable + "_" + trackType;
            timeVariable = timeVariable + "_" + trackType;
        end
        for diameterIndex = 1:numel(fasanDMrad)
            if beamTypes(beamIndex) == "POINT"
                fileName = sprintf('%s_D%gmrad.mat', ...
                    fileStem, fasanDMrad(diameterIndex));
            else
                fileName = sprintf('%s_D%gd%gmrad.mat', ...
                    fileStem, fasanDMrad(diameterIndex), smallDMrad);
            end
            resultFile = fullfile(resultDir, fileName);
            if ~isfile(resultFile)
                warning('tmp:MissingResult', 'Missing MC result: %s', resultFile);
                continue;
            end

            result = load(resultFile, detectVariable, timeVariable);
            if ~isfield(result, detectVariable) || ~isfield(result, timeVariable)
                warning('tmp:MissingVariables', ...
                    'Missing detection/time variables in: %s', resultFile);
                continue;
            end
            detectArray = result.(detectVariable);
            timeArray = result.(timeVariable);
            if ~(isnumeric(detectArray) || islogical(detectArray)) || ...
                    ~isnumeric(timeArray) || ~ismatrix(detectArray) || ...
                    ~ismatrix(timeArray) || isempty(detectArray) || ...
                    size(detectArray, 1) ~= size(timeArray, 1) || ...
                    size(detectArray, 2) < columnIndex || ...
                    size(timeArray, 2) < columnIndex
                warning('tmp:InvalidArrays', ...
                    'Invalid detection/time arrays in: %s', resultFile);
                continue;
            end
            detect = detectArray(:, columnIndex);
            firstTime = timeArray(:, columnIndex);
            if any(~isfinite(detect) | (detect ~= 0 & detect ~= 1))
                warning('tmp:InvalidDetection', ...
                    'Detection values must be finite 0 or 1 in: %s', resultFile);
                continue;
            end

            % Same statistics as summarizeMCResult in MC_overall_plot.m.
            % first_time_* already contains seconds: do not divide by f again.
            isDetected = detect ~= 0;
            probability(beamIndex, diameterIndex) = mean(isDetected);
            hasValidTime = isDetected & isfinite(firstTime);
            if any(hasValidTime)
                firstWarningTime(beamIndex, diameterIndex) = ...
                    mean(firstTime(hasValidTime));
            end
            loadedCount = loadedCount + 1;
        end
    end
end

function formatMCAxes(axesHandle, fasanDMrad, panelTitle)
    grid(axesHandle, 'on');
    box(axesHandle, 'on');
    title(axesHandle, panelTitle, 'Interpreter', 'none', ...
        'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
    xlim(axesHandle, [0, max(fasanDMrad)]);
    xticks(axesHandle, 0:100:max(fasanDMrad));
    xtickangle(axesHandle, 0);
    set(axesHandle, 'FontName', 'Helvetica', ...
        'FontUnits', 'points', 'FontSize', 14, ...
        'TitleFontSizeMultiplier', 1, 'LabelFontSizeMultiplier', 1, ...
        'LineWidth', 2, 'XScale', 'linear');
    % Set title size after the axes to avoid automatic title enlargement.
    set([axesHandle.Title, axesHandle.XLabel, axesHandle.YLabel], ...
        'FontName', 'Helvetica', 'FontUnits', 'points', 'FontSize', 14);
end

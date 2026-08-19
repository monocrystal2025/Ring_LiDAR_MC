clc;
clear;
close all;

% Mechanism comparison for EA, ROI raster (RA), and ROI spiral (SP) scans.
% Row 1 shows the spatial coverage rate over one complete scan cycle.
% Row 2 shows the mean first effective encounter time.
% Row 3 shows the mean effective pulse count for detected trials.
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

eaResultDir = 'D:\lzx\MC_RESULTS';
roiResultDir = fullfile(scriptDir, 'MC_ROI_RESULTS_SNR_NEW');
staticEaResultDir = 'D:\matlab\STATIC_MC_RESULTS';
staticRoiResultDir = fullfile(scriptDir, 'STATIC_MC_ROI_RESULTS_NEW');
assert(isfolder(eaResultDir), 'EA result folder not found: %s', eaResultDir);
assert(isfolder(roiResultDir), 'ROI result folder not found: %s', roiResultDir);
assert(isfolder(staticEaResultDir), ...
    'Static EA result folder not found: %s', staticEaResultDir);
assert(isfolder(staticRoiResultDir), ...
    'Static ROI result folder not found: %s', staticRoiResultDir);

fasanDList = (5:5:400) * 1e-3;
fasanDMrad = fasanDList * 1e3;
fasanDInner = 1e-3;
columnIndex = 1;

columnCases = struct( ...
    'family', {'EA', 'ROI', 'ROI'}, ...
    'path', {'', 'RA', 'SP'}, ...
    'title', {'ROE', 'ROI-RA', 'ROI-SP'}, ...
    'resultDir', {eaResultDir, roiResultDir, roiResultDir}, ...
    'staticResultDir', ...
    {staticEaResultDir, staticRoiResultDir, staticRoiResultDir});

nColumns = numel(columnCases);
coverageRate = cell(1, nColumns);
meanEncounterTime = cell(1, nColumns);
meanEffectivePulses = cell(1, nColumns);

for columnIdx = 1:nColumns
    coverageRate{columnIdx} = loadCoverageColumn( ...
        columnCases(columnIdx), fasanDMrad, fasanDInner);
    [meanEffectivePulses{columnIdx}, meanEncounterTime{columnIdx}] = ...
        loadMechanismColumn(columnCases(columnIdx), fasanDMrad, ...
        fasanDInner, columnIndex);
end

figureHandle = figure('Color', 'w', 'Position', [50, 50, 620, 580]);
layout = tiledlayout(figureHandle, 3, 3, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

beamLabels = {'Spot', 'Line', 'Annular'};
beamStyles = {'--', ':', '-'};
beamColors = [ ...
    0.4660, 0.6740, 0.1880; ...
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980];

coverageAxes = gobjects(1, nColumns);
encounterAxes = gobjects(1, nColumns);
pulseAxes = gobjects(1, nColumns);
pulsePlots = gobjects(1, numel(beamLabels));

for columnIdx = 1:nColumns
    coverageAxes(columnIdx) = nexttile(layout, columnIdx);
    plotBeamSeries(coverageAxes(columnIdx), fasanDMrad, ...
        coverageRate{columnIdx}, beamLabels, beamStyles, beamColors);
    configureMechanismAxes(coverageAxes(columnIdx), 100);
    ylim(coverageAxes(columnIdx), [70, 100]);
    yticks(coverageAxes(columnIdx), 60:10:100);
    xticklabels(coverageAxes(columnIdx), []);
    if columnIdx ~= 1
        yticklabels(coverageAxes(columnIdx), []);
    end
    title(coverageAxes(columnIdx), columnCases(columnIdx).title, ...
        'FontName', 'Helvetica', 'FontSize', 13, 'FontWeight', 'normal');
    if columnIdx == 1
        ylabel(coverageAxes(columnIdx), ...
            '{\it P}_{\rmcov} (%)', ...
            'FontName', 'Helvetica');
    end

    encounterAxes(columnIdx) = nexttile(layout, 3 + columnIdx);
    plotBeamSeries(encounterAxes(columnIdx), fasanDMrad, ...
        meanEncounterTime{columnIdx}, beamLabels, beamStyles, beamColors);
    configureMechanismAxes(encounterAxes(columnIdx), 40);
    yticks(encounterAxes(columnIdx), 0:10:40);
    xticklabels(encounterAxes(columnIdx), []);
    if columnIdx ~= 1
        yticklabels(encounterAxes(columnIdx), []);
    end
    if columnIdx == 1
        ylabel(encounterAxes(columnIdx), ...
            '{\itt}_{\rmenc} (s)', ...
            'FontName', 'Helvetica');
    end

    pulseAxes(columnIdx) = nexttile(layout, 6 + columnIdx);
    pulsePlotsCurrent = plotBeamSeries(pulseAxes(columnIdx), fasanDMrad, ...
        meanEffectivePulses{columnIdx}, beamLabels, beamStyles, beamColors);
    configureMechanismAxes(pulseAxes(columnIdx), 1000);
    set(pulseAxes(columnIdx), 'YScale', 'log');
    ylim(pulseAxes(columnIdx), [1, 1000]);
    yticks(pulseAxes(columnIdx), [1, 10, 100, 1000]);
    if columnIdx ~= 1
        yticklabels(pulseAxes(columnIdx), []);
    end
    set(pulseAxes(columnIdx), ...
        'XGrid', 'on', 'YGrid', 'on', ...
        'XMinorGrid', 'off', 'YMinorGrid', 'off');
    xlabel(pulseAxes(columnIdx), ...
        '{\itw}_{\rmD} (mrad)', 'Interpreter', 'tex');
    if columnIdx == 1
        ylabel(pulseAxes(columnIdx), ...
            '{\itN}_{\rmeff}', ...
            'FontName', 'Helvetica');
        pulsePlots = pulsePlotsCurrent;
    end
end

totalLegend = legend(pulseAxes(1), pulsePlots, beamLabels, ...
    'Orientation', 'horizontal', 'Box', 'off', ...
    'FontName', 'Helvetica', 'FontSize', 13);
totalLegend.Layout.Tile = 'north';

exportgraphics(figureHandle, fullfile(scriptDir, 'mechanism.png'), ...
    'Resolution', 300);
savefig(figureHandle, fullfile(scriptDir, 'mechanism.fig'));

function coverageRate = loadCoverageColumn( ...
        caseInfo, fasanDMrad, fasanDInner)
beamNames = {'POINT', 'LINE', 'RING'};
beamCodes = {'P', 'L', 'R'};
nD = numel(fasanDMrad);
nBeams = numel(beamNames);
coverageRate = nan(nD, nBeams);
dText = num2str(fasanDInner * 1e3, '%g');

for beamIdx = 1:nBeams
    beamName = beamNames{beamIdx};
    beamCode = beamCodes{beamIdx};

    for k = 1:nD
        DText = num2str(fasanDMrad(k), '%g');
        if strcmp(caseInfo.family, 'EA')
            if strcmp(beamName, 'POINT')
                fileName = sprintf( ...
                    'STATIC_MC_1par_EA_POINT_D%smrad.mat', DText);
            else
                fileName = sprintf( ...
                    'STATIC_MC_1par_EA_%s_D%sd%smrad.mat', ...
                    beamName, DText, dText);
            end
            detectVariable = ['detect_', beamCode];
        else
            if strcmp(beamName, 'POINT')
                fileName = sprintf( ...
                    'STATIC_MC_1par_ROI_POINT_%s_D%smrad.mat', ...
                    caseInfo.path, DText);
            else
                fileName = sprintf( ...
                    'STATIC_MC_1par_ROI_%s_%s_D%sd%smrad.mat', ...
                    beamName, caseInfo.path, DText, dText);
            end
            detectVariable = sprintf('detect_%s_%s', ...
                beamCode, caseInfo.path);
        end

        resultFile = fullfile(caseInfo.staticResultDir, fileName);
        assert(isfile(resultFile), ...
            'Missing static result file: %s', resultFile);
        result = load(resultFile, detectVariable);
        assert(isfield(result, detectVariable), ...
            'Missing variable %s in %s.', detectVariable, resultFile);

        detect = result.(detectVariable);
        assert(isnumeric(detect) || islogical(detect), ...
            'Variable %s in %s must be numeric or logical.', ...
            detectVariable, resultFile);
        validDetect = isfinite(detect);
        assert(any(validDetect, 'all'), ...
            'Variable %s in %s contains no finite samples.', ...
            detectVariable, resultFile);
        coverageRate(k, beamIdx) = ...
            100 * mean(detect(validDetect) ~= 0);
    end
end
end

function [meanEffectivePulses, meanEncounterTime] = ...
        loadMechanismColumn(caseInfo, fasanDMrad, fasanDInner, columnIndex)
beamNames = {'POINT', 'LINE', 'RING'};
beamCodes = {'P', 'L', 'R'};
nD = numel(fasanDMrad);
nBeams = numel(beamNames);
meanEffectivePulses = nan(nD, nBeams);
meanEncounterTime = nan(nD, nBeams);
dText = num2str(fasanDInner * 1e3, '%g');

for beamIdx = 1:nBeams
    beamName = beamNames{beamIdx};
    beamCode = beamCodes{beamIdx};

    for k = 1:nD
        DText = num2str(fasanDMrad(k), '%g');
        if strcmp(caseInfo.family, 'EA')
            if strcmp(beamName, 'POINT')
                fileName = sprintf('MC_1par_EA_POINT_D%smrad.mat', DText);
            else
                fileName = sprintf('MC_1par_EA_%s_D%sd%smrad.mat', ...
                    beamName, DText, dText);
            end
            variableSuffix = beamCode;
        else
            if strcmp(beamName, 'POINT')
                fileName = sprintf('MC_1par_ROI_NEW_POINT_%s_D%smrad.mat', ...
                    caseInfo.path, DText);
            else
                fileName = sprintf( ...
                    'MC_1par_ROI_NEW_%s_%s_D%sd%smrad.mat', ...
                    beamName, caseInfo.path, DText, dText);
            end
            variableSuffix = sprintf('%s_%s', beamCode, caseInfo.path);
        end

        resultFile = fullfile(caseInfo.resultDir, fileName);
        assert(isfile(resultFile), 'Missing result file: %s', resultFile);

        detectVariable = ['detect_', variableSuffix];
        pulseVariable = ['effective_pulses_', variableSuffix];
        encounterVariable = ['first_encounter_time_', variableSuffix];
        result = load(resultFile, detectVariable, pulseVariable, ...
            encounterVariable);

        assert(isfield(result, detectVariable) && ...
            isfield(result, pulseVariable) && ...
            isfield(result, encounterVariable), ...
            'Missing mechanism variables in %s.', resultFile);

        detect = selectResultColumn(result.(detectVariable), ...
            columnIndex, detectVariable, resultFile);
        effectivePulses = selectResultColumn(result.(pulseVariable), ...
            columnIndex, pulseVariable, resultFile);
        firstEncounterTime = selectResultColumn( ...
            result.(encounterVariable), columnIndex, ...
            encounterVariable, resultFile);

        [meanEffectivePulses(k, beamIdx), ...
            meanEncounterTime(k, beamIdx)] = summarizeMechanism( ...
            detect, effectivePulses, firstEncounterTime);
    end
end
end

function values = selectResultColumn(values, columnIndex, variableName, fileName)
assert(ismatrix(values) && size(values, 2) >= columnIndex, ...
    'Variable %s in %s does not contain column %d.', ...
    variableName, fileName, columnIndex);
values = values(:, columnIndex);
end

function plotHandles = plotBeamSeries(axesHandle, x, values, ...
        labels, lineStyles, colors)
plotHandles = gobjects(1, size(values, 2));
hold(axesHandle, 'on');
for beamIdx = 1:size(values, 2)
    plotHandles(beamIdx) = plot(axesHandle, x, values(:, beamIdx), ...
        'LineStyle', lineStyles{beamIdx}, ...
        'Color', colors(beamIdx, :), ...
        'LineWidth', 3, ...
        'DisplayName', labels{beamIdx});
end
end

function configureMechanismAxes(axesHandle, yUpper)
grid(axesHandle, 'on');
box(axesHandle, 'on');
xlim(axesHandle, [0, 400]);
ylim(axesHandle, [0, yUpper]);
xticks(axesHandle, 0:100:400);
set(axesHandle, 'FontName', 'Helvetica', ...
    'FontSize', 13, 'LineWidth', 2, ...
    'XTickLabelRotation', 0);
set([axesHandle.XLabel, axesHandle.YLabel], ...
    'FontName', 'Helvetica', 'FontSize', 13);
end

function [meanEffectivePulses, meanEncounterTime] = ...
        summarizeMechanism(detect, effectivePulses, firstEncounterTime)
detect = detect(:);
effectivePulses = effectivePulses(:);
firstEncounterTime = firstEncounterTime(:);

assert(numel(detect) == numel(effectivePulses) && ...
    numel(detect) == numel(firstEncounterTime), ...
    'Detection, pulse-count, and encounter-time arrays must be equal in length.');

validPulseCount = detect ~= 0 & isfinite(effectivePulses);
validEncounterTime = isfinite(firstEncounterTime) & firstEncounterTime >= 0;

if any(validPulseCount)
    meanEffectivePulses = mean(effectivePulses(validPulseCount));
else
    meanEffectivePulses = nan;
end

if any(validEncounterTime)
    meanEncounterTime = mean(firstEncounterTime(validEncounterTime));
else
    meanEncounterTime = nan;
end
end

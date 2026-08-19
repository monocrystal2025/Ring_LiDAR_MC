clc; clear; close all;

% Plot MC_ROI_NEW result curves using the style of plot_MC_EA_fasanD.m.
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

resultDir = fullfile(scriptDir, 'MC_ROI_RESULTS_SNR_NEW');
if ~isfolder(resultDir)
    error('Result folder not found: %s', resultDir);
end

fasan_D_mrad = 5:5:400;
fasan_d_mrad = 1;
trackTypes = {'RA', 'SP'};
beamTypes = {'POINT', 'LINE', 'RING'};
beamLabels = {'Spot', 'Line', 'Annular'};
markers = {'^', 's', 'o'};
colors = [ ...
    0.4660, 0.6740, 0.1880; ...
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980];

nTrack = numel(trackTypes);
nBeam = numel(beamTypes);
probability = nan(nBeam, nTrack, numel(fasan_D_mrad));
firstWarningTime = nan(nBeam, nTrack, numel(fasan_D_mrad));

for trackIdx = 1:nTrack
    trackType = trackTypes{trackIdx};
    for beamIdx = 1:nBeam
        beamType = beamTypes{beamIdx};
        [detectVar, firstTimeVar] = resultVariableNames( ...
            beamType, trackType);

        for k = 1:numel(fasan_D_mrad)
            resultFile = resultFileName(resultDir, beamType, ...
                trackType, fasan_D_mrad(k), fasan_d_mrad);

            if ~isfile(resultFile)
                warning('Missing ROI result: %s', resultFile);
                continue;
            end

            S = load(resultFile, detectVar, firstTimeVar);
            if ~isfield(S, detectVar) || ~isfield(S, firstTimeVar)
                warning('Missing expected variables in: %s', resultFile);
                continue;
            end

            [probability(beamIdx, trackIdx, k), ...
                firstWarningTime(beamIdx, trackIdx, k)] = ...
                summarizeMCResult(S.(detectVar), S.(firstTimeVar));
        end
    end
end

fig = figure('Color', 'w', 'Position', [100, 100, 830, 550]);
t = tiledlayout(fig, 2, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

probAxes = gobjects(1, nTrack);
timeAxes = gobjects(1, nTrack);
legendLines = gobjects(1, nBeam);

for trackIdx = 1:nTrack
    probAxes(trackIdx) = nexttile(t, trackIdx);
    hold(probAxes(trackIdx), 'on');
    for beamIdx = 1:nBeam
        y = 100 * reshape(probability(beamIdx, trackIdx, :), 1, []);
        h = plot(probAxes(trackIdx), fasan_D_mrad, y, ...
            ['-' markers{beamIdx}], 'LineWidth', 1.8, ...
            'MarkerSize', 4, 'Color', colors(beamIdx, :), ...
            'DisplayName', beamLabels{beamIdx});
        if trackIdx == 1
            legendLines(beamIdx) = h;
        end
    end
    formatProbabilityAxes(probAxes(trackIdx), trackTypes{trackIdx});
    xticklabels(probAxes(trackIdx), []);
    xlabel(probAxes(trackIdx), '');
    if trackIdx == nTrack
        yticklabels(probAxes(trackIdx), []);
        ylabel(probAxes(trackIdx), '');
    end

    timeAxes(trackIdx) = nexttile(t, nTrack + trackIdx);
    hold(timeAxes(trackIdx), 'on');
    for beamIdx = 1:nBeam
        y = reshape(firstWarningTime(beamIdx, trackIdx, :), 1, []);
        plot(timeAxes(trackIdx), fasan_D_mrad, y, ...
            ['-' markers{beamIdx}], 'LineWidth', 1.8, ...
            'MarkerSize', 4, 'Color', colors(beamIdx, :), ...
            'DisplayName', beamLabels{beamIdx});
    end
    formatTimeAxes(timeAxes(trackIdx));
    if trackIdx == nTrack
        yticklabels(timeAxes(trackIdx), []);
        ylabel(timeAxes(trackIdx), '');
    end
end

linkaxes(probAxes, 'xy');
linkaxes(timeAxes, 'xy');

lgd = legend(probAxes(1), legendLines, beamLabels, ...
    'Location', 'best', 'FontName', 'Helvetica', 'FontSize', 18);
lgd.Layout.Tile = 'north';
lgd.NumColumns = nBeam;
lgd.Box = 'off';
lgd.ItemTokenSize = [50, 18];

outPng = fullfile(scriptDir, 'MC_ROI_fasan_D_summary.png');
outFig = fullfile(scriptDir, 'MC_ROI_fasan_D_summary.fig');
exportgraphics(fig, outPng, 'Resolution', 300);
savefig(fig, outFig);

disp(['Saved PNG: ', outPng]);
disp(['Saved FIG: ', outFig]);

function resultFile = resultFileName(resultDir, beamType, ...
    trackType, D_mrad, d_mrad)
DText = num2str(D_mrad, '%g');

if strcmp(beamType, 'POINT')
    fileName = sprintf('MC_1par_ROI_NEW_POINT_%s_D%smrad.mat', ...
        trackType, DText);
else
    dText = num2str(d_mrad, '%g');
    fileName = sprintf('MC_1par_ROI_NEW_%s_%s_D%sd%smrad.mat', ...
        beamType, trackType, DText, dText);
end

resultFile = fullfile(resultDir, fileName);
end

function [detectVar, firstTimeVar] = resultVariableNames( ...
    beamType, trackType)
switch beamType
    case 'POINT'
        beamCode = 'P';
    case 'LINE'
        beamCode = 'L';
    case 'RING'
        beamCode = 'R';
    otherwise
        error('Unsupported beam type: %s', beamType);
end

detectVar = sprintf('detect_%s_%s', beamCode, trackType);
firstTimeVar = sprintf('first_time_%s_%s', beamCode, trackType);
end

function [probability, firstWarningTime] = summarizeMCResult( ...
    detect, firstTime)
detect = detect(:);
firstTime = firstTime(:);

if numel(detect) ~= numel(firstTime)
    error('Detection and first-warning-time arrays have different sizes.');
end

isDetected = detect ~= 0;
probability = mean(isDetected);

validTime = isDetected & isfinite(firstTime);
if any(validTime)
    firstWarningTime = mean(firstTime(validTime));
else
    firstWarningTime = nan;
end
end

function formatProbabilityAxes(ax, trackType)
grid(ax, 'on');
box(ax, 'on');
title(ax, trackType, 'FontName', 'Helvetica', 'FontSize', 18, ...
    'FontWeight', 'normal');
xlabel(ax, '{\it w}_{\rmD} (mrad)', 'Interpreter', 'tex');
ylabel(ax, '{\it P}_{\rmdetect} (%)', 'Interpreter', 'tex');
xticks(ax, 0:100:500);
yticks(ax, 0:25:100);
ylim(ax, [0, 100]);
set(ax, 'FontName', 'Helvetica', 'FontSize', 18, ...
    'LineWidth', 2, 'XScale', 'linear');
set([ax.XLabel, ax.YLabel], 'FontName', 'Helvetica', 'FontSize', 18);
end

function formatTimeAxes(ax)
grid(ax, 'on');
box(ax, 'on');
xlabel(ax, '{\it w}_{\rmD} (mrad)', 'Interpreter', 'tex');
ylabel(ax, ['{\it t' char(773) '}_{\rm warn} (s)'], ...
    'Interpreter', 'tex');
xticks(ax, 0:100:500);
set(ax, 'FontName', 'Helvetica', 'FontSize', 18, ...
    'LineWidth', 2, 'XScale', 'linear');
set([ax.XLabel, ax.YLabel], 'FontName', 'Helvetica', 'FontSize', 18);
end

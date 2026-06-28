clc; clear; close all;

% Plot MC_EA result curves from generated .mat files.
resultDir = 'D:\lzx\MC_RESULTS\';
fasan_D_LIST = 5e-3:10e-3:500e-3;
fasan_d = 0.8e-3;
colIdx = 1;

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

fasan_D_mrad = fasan_D_LIST * 1e3;
probLine = nan(size(fasan_D_LIST));
probRing = nan(size(fasan_D_LIST));
timeLine = nan(size(fasan_D_LIST));
timeRing = nan(size(fasan_D_LIST));

for k = 1:numel(fasan_D_LIST)
    DText = num2str(fasan_D_mrad(k), '%g');
    dText = num2str(fasan_d * 1e3, '%g');

    ringFile = fullfile(resultDir, ...
        ['MC_1par_EA_RING_D', DText, 'd', dText, 'mrad.mat']);
    lineFile = fullfile(resultDir, ...
        ['MC_1par_EA_LINE_D', DText, 'd', dText, 'mrad.mat']);

    if isfile(ringFile)
        S = load(ringFile, 'detect_R', 'first_time_R');
        [probRing(k), timeRing(k)] = summarizeMCResult( ...
            S.detect_R(:, colIdx), S.first_time_R(:, colIdx));
    else
        warning('Missing ring result: %s', ringFile);
    end

    if isfile(lineFile)
        S = load(lineFile, 'detect_L', 'first_time_L');
        [probLine(k), timeLine(k)] = summarizeMCResult( ...
            S.detect_L(:, colIdx), S.first_time_L(:, colIdx));
    else
        warning('Missing line result: %s', lineFile);
    end
end

fig = figure('Color', 'w', 'Position', [100, 100, 1100, 800]);

subplot(2, 1, 1);
pLineProb = plot(fasan_D_mrad, probLine, '-o', 'LineWidth', 1.8, ...
    'MarkerSize', 5, 'DisplayName', '线光');
hold on;
pRingProb = plot(fasan_D_mrad, probRing, '-s', 'LineWidth', 1.8, ...
    'MarkerSize', 5, 'DisplayName', '环光');
markCurveExtremum(fasan_D_mrad, probLine, 'max', pLineProb.Color, '线光最大值');
markCurveExtremum(fasan_D_mrad, probRing, 'max', pRingProb.Color, '环光最大值');
grid on;
box on;
xlabel('fasan_D (mrad)', 'Interpreter', 'none');
ylabel('概率');
title('概率 - fasan_D', 'Interpreter', 'none');
legend('Location', 'best');

subplot(2, 1, 2);
pLineTime = plot(fasan_D_mrad, timeLine, '-o', 'LineWidth', 1.8, ...
    'MarkerSize', 5, 'DisplayName', '线光');
hold on;
pRingTime = plot(fasan_D_mrad, timeRing, '-s', 'LineWidth', 1.8, ...
    'MarkerSize', 5, 'DisplayName', '环光');
markCurveExtremum(fasan_D_mrad, timeLine, 'min', pLineTime.Color, '线光最小值');
markCurveExtremum(fasan_D_mrad, timeRing, 'min', pRingTime.Color, '环光最小值');
grid on;
box on;
xlabel('fasan_D (mrad)', 'Interpreter', 'none');
ylabel('首次预警时长 (s)');
title('首次预警时长 - fasan_D', 'Interpreter', 'none');
legend('Location', 'best');

sgtitle('MC_EA 仿真结果', 'Interpreter', 'none');

outPng = fullfile(scriptDir, 'MC_EA_fasan_D_summary.png');
outFig = fullfile(scriptDir, 'MC_EA_fasan_D_summary.fig');
exportgraphics(fig, outPng, 'Resolution', 300);
savefig(fig, outFig);

disp(['Saved PNG: ', outPng]);
disp(['Saved FIG: ', outFig]);

function [probability, firstWarningTime] = summarizeMCResult(detect, firstTime)
detect = detect(:);
firstTime = firstTime(:);

isDetected = detect ~= 0;
probability = mean(isDetected);

validTime = isDetected & isfinite(firstTime);
if any(validTime)
    firstWarningTime = mean(firstTime(validTime));
else
    firstWarningTime = nan;
end
end

function markCurveExtremum(x, y, mode, color, labelText)
x = x(:);
y = y(:);
valid = isfinite(x) & isfinite(y);
if ~any(valid)
    return;
end

validIdx = find(valid);
switch lower(mode)
    case 'max'
        [value, localIdx] = max(y(valid));
        verticalAlignment = 'bottom';
        ySign = 1;
    case 'min'
        [value, localIdx] = min(y(valid));
        verticalAlignment = 'top';
        ySign = -1;
    otherwise
        error('Unknown mode: %s. Use max or min.', mode);
end

idx = validIdx(localIdx);
xValue = x(idx);
xValid = x(valid);
xRange = max(xValid) - min(xValid);
if xRange == 0
    xRange = 1;
end

yLimits = ylim;
yOffset = 0.035 * (max(yLimits) - min(yLimits));
if yOffset == 0
    yOffset = 0.05 * max(abs(value), 1);
end

plot(xValue, value, 'p', ...
    'MarkerSize', 15, ...
    'MarkerFaceColor', color, ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.1, ...
    'HandleVisibility', 'off');

text(xValue + 0.012 * xRange, value + ySign * yOffset, ...
    sprintf('%s\nD=%.0f mrad, %.3g', labelText, xValue, value), ...
    'Color', color, ...
    'FontSize', 10, ...
    'FontWeight', 'bold', ...
    'VerticalAlignment', verticalAlignment, ...
    'BackgroundColor', 'w', ...
    'Margin', 2, ...
    'Interpreter', 'none', ...
    'HandleVisibility', 'off');
end

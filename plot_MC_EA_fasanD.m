clc; clear; close all;

% Plot MC_EA result curves from generated .mat files.
resultDir = 'D:\lzx\MC_RESULTS\';
% resultDir = 'G:\MC_RESULTS\正常结果\MC_RESULTS\';
fasan_D_LIST = 5e-3:5e-3:400e-3;
fasan_d = 1e-3;
colIdx = 1;

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

fasan_D_mrad = fasan_D_LIST * 1e3;
probLine = nan(size(fasan_D_LIST));
probRing = nan(size(fasan_D_LIST));
probPointRaw = nan(size(fasan_D_LIST));
timeLine = nan(size(fasan_D_LIST));
timeRing = nan(size(fasan_D_LIST));
timePointRaw = nan(size(fasan_D_LIST));

for k = 1:numel(fasan_D_LIST)
    DText = num2str(fasan_D_mrad(k), '%g');
    dText = num2str(fasan_d * 1e3, '%g');

    ringFile = fullfile(resultDir, ...
        ['MC_1par_EA_RING_D', DText, 'd', dText, 'mrad.mat']);
    lineFile = fullfile(resultDir, ...
        ['MC_1par_EA_LINE_D', DText, 'd', dText, 'mrad.mat']);
    pointFile = fullfile(resultDir, ...
        ['MC_1par_EA_POINT_D', DText, 'mrad.mat']);

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

    if isfile(pointFile)
        S = load(pointFile, 'detect_P', 'first_time_P');
        [probPointRaw(k), timePointRaw(k)] = summarizeMCResult( ...
            S.detect_P(:, colIdx), S.first_time_P(:, colIdx));
    end
end

probPoint = fillMissingByInterpolation(fasan_D_mrad, probPointRaw, 'point probability');
timePoint = fillMissingByInterpolation(fasan_D_mrad, timePointRaw, 'point first warning time');

fig = figure('Position', [100, 100, 530, 600]);
t = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

axProb = nexttile(t);
pLineProb = plot(fasan_D_mrad, 100 * probLine, '-o', 'LineWidth', 1.8, ...
    'MarkerSize', 5, 'DisplayName', 'line');
hold on;
pRingProb = plot(fasan_D_mrad, 100 * probRing, '-s', 'LineWidth', 1.8, ...
    'MarkerSize', 5, 'DisplayName', 'annular');
pPointProb = plot(fasan_D_mrad, 100 * probPoint, '-^', 'LineWidth', 1.8, ...
    'MarkerSize', 5, 'Color', [0.4660, 0.6740, 0.1880], ...
    'DisplayName', 'spot');
grid on;
box on;
xlabel('{\it w}_{\rmD} (mrad)', 'Interpreter', 'tex');
ylabel('{\it P}_{\rmdetect} (%)', 'Interpreter', 'tex');
xticks(0:100:500);
yticks(0:25:100);
set(axProb, 'FontName', 'Helvetica', 'FontSize', 13, 'LineWidth', 2,'XScale','linear');
set([axProb.XLabel, axProb.YLabel], 'FontName', 'Helvetica', 'FontSize', 13);

axTime = nexttile(t);
pLineTime = plot(fasan_D_mrad, timeLine, '-o', 'LineWidth', 1.8, ...
    'MarkerSize', 5, 'DisplayName', 'line');
hold on;
pRingTime = plot(fasan_D_mrad, timeRing, '-s', 'LineWidth', 1.8, ...
    'MarkerSize', 5, 'DisplayName', 'annular');
pPointTime = plot(fasan_D_mrad, timePoint, '-^', 'LineWidth', 1.8, ...
    'MarkerSize', 5, 'Color', pPointProb.Color, ...
    'DisplayName', 'spot');
grid on;
box on;
xlabel('{\it w}_{\rmD} (mrad)', 'Interpreter', 'tex');
ylabel(['{\it t' char(773) '}_{\rm warn} (s)'], 'Interpreter', 'tex');
xticks(0:100:500);
set(axTime, 'FontName', 'Helvetica', 'FontSize', 13, 'LineWidth', 2);
set([axTime.XLabel, axTime.YLabel], 'FontName', 'Helvetica', 'FontSize', 13);
lgd = legend(axProb, [pPointProb, pLineProb, pRingProb], ...
    {'spot', 'line', 'annular'}, 'Location', 'best', ...
    'FontName', 'Helvetica', 'FontSize', 13);
lgd.Layout.Tile = 'north';
lgd.NumColumns = 3;
lgd.Box = 'off';
lgd.ItemTokenSize = [50, 18];

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

function yFilled = fillMissingByInterpolation(x, yRaw, dataName)
x = x(:);
yRaw = yRaw(:);
yFilled = yRaw;
valid = isfinite(x) & isfinite(yRaw);

if ~any(valid)
    warning('No valid %s data found. The curve will remain NaN.', dataName);
    return;
end

if nnz(valid) == 1
    yFilled(:) = yRaw(valid);
    warning('Only one valid %s data point found. Filled missing values with that value.', dataName);
    return;
end

missing = ~isfinite(yFilled);
if any(missing)
    yFilled(missing) = interp1(x(valid), yRaw(valid), x(missing), 'pchip', 'extrap');
    fprintf('Filled %d missing %s points by pchip interpolation.\n', nnz(missing), dataName);
end

yFilled = reshape(yFilled, size(yRaw));
end

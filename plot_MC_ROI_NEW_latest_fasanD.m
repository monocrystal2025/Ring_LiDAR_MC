clc; clear; close all;

% Plot the latest MC_ROI_NEW results saved in MC_ROI_RESULTS_SNR_NEW.
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

resultDir = fullfile(scriptDir, 'MC_ROI_RESULTS_SNR_NEW');
if ~isfolder(resultDir)
    error('Result folder not found: %s', resultDir);
end

cases = struct( ...
    'beam', {'ring', 'ring', 'line', 'line'}, ...
    'path', {'RA', 'SP', 'RA', 'SP'}, ...
    'label', {'环光-栅格轨迹', '环光-螺旋轨迹', '线光-栅格轨迹', '线光-螺旋轨迹'}, ...
    'detectVar', {'detect_R_RA', 'detect_R_SP', 'detect_L_RA', 'detect_L_SP'}, ...
    'firstVar', {'first_time_R_RA', 'first_time_R_SP', 'first_time_L_RA', 'first_time_L_SP'}, ...
    'lineStyle', {'-', '--', '-', '--'}, ...
    'marker', {'o', 's', '^', 'd'});

colors = [ ...
    0.0000, 0.4470, 0.7410; ...
    0.3010, 0.7450, 0.9330; ...
    0.8500, 0.3250, 0.0980; ...
    0.9290, 0.6940, 0.1250];

allRows = table();
series = repmat(struct('D_mrad', [], 'probability', [], ...
    'meanFirstWarningTime', []), size(cases));

for idx = 1:numel(cases)
    [D_mrad, probability, meanFirstWarningTime, nSamples, nDetected, ...
        nFiniteFirst, sourceFile] = loadCaseSummary(resultDir, cases(idx));

    series(idx).D_mrad = D_mrad;
    series(idx).probability = probability;
    series(idx).meanFirstWarningTime = meanFirstWarningTime;

    rows = table( ...
        repmat(string(cases(idx).beam), numel(D_mrad), 1), ...
        repmat(string(cases(idx).path), numel(D_mrad), 1), ...
        repmat(string(cases(idx).label), numel(D_mrad), 1), ...
        D_mrad(:), probability(:), meanFirstWarningTime(:), ...
        nSamples(:), nDetected(:), nFiniteFirst(:), string(sourceFile(:)), ...
        'VariableNames', {'Beam', 'Path', 'Label', 'D_mrad', ...
        'Probability', 'MeanFirstWarningTime_s', 'NSamples', ...
        'NDetected', 'NFiniteFirstTime', 'SourceFile'});

    allRows = [allRows; rows]; %#ok<AGROW>
end

fig = figure('Color', 'w', 'Position', [80, 80, 1120, 820]);
layout = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile(layout, 1);
hold on;
for idx = 1:numel(cases)
    plot(series(idx).D_mrad, series(idx).probability, ...
        'LineStyle', cases(idx).lineStyle, ...
        'Marker', cases(idx).marker, ...
        'Color', colors(idx, :), ...
        'MarkerFaceColor', 'w', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 5.5, ...
        'DisplayName', cases(idx).label);
end
grid on; box on;
xlabel('fasan\_D (mrad)', 'Interpreter', 'tex');
ylabel('检测概率');
title('检测概率 - fasan\_D', 'Interpreter', 'tex');
legend('Location', 'best');

nexttile(layout, 2);
hold on;
for idx = 1:numel(cases)
    plot(series(idx).D_mrad, series(idx).meanFirstWarningTime, ...
        'LineStyle', cases(idx).lineStyle, ...
        'Marker', cases(idx).marker, ...
        'Color', colors(idx, :), ...
        'MarkerFaceColor', 'w', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 5.5, ...
        'DisplayName', cases(idx).label);
end
grid on; box on;
xlabel('fasan\_D (mrad)', 'Interpreter', 'tex');
ylabel('平均首次预警时长 (s)');
title('平均首次预警时长 - fasan\_D', 'Interpreter', 'tex');
legend('Location', 'best');

title(layout, 'MC\_ROI\_NEW 最新仿真结果：线光与环光对比', 'Interpreter', 'tex');

outBase = fullfile(resultDir, 'roi_new_latest_line_ring_fasan_D_summary');
outCsv = [outBase, '.csv'];
outPng = [outBase, '.png'];
outFig = [outBase, '.fig'];

writetable(sortrows(allRows, {'Beam', 'Path', 'D_mrad'}), outCsv);
exportgraphics(fig, outPng, 'Resolution', 300);
savefig(fig, outFig);

fprintf('Saved CSV: %s\n', outCsv);
fprintf('Saved PNG: %s\n', outPng);
fprintf('Saved FIG: %s\n', outFig);
fprintf('Plotted %d result curves across %d source files.\n', ...
    numel(cases), height(allRows));

function [D_mrad, probability, meanFirstWarningTime, nSamples, nDetected, ...
    nFiniteFirst, sourceFile] = loadCaseSummary(resultDir, caseInfo)
filePattern = sprintf('MC_1par_ROI_NEW_%s_%s_D*d*mrad.mat', ...
    upper(caseInfo.beam), upper(caseInfo.path));
files = dir(fullfile(resultDir, filePattern));
if isempty(files)
    error('No result files found for %s-%s in %s.', ...
        caseInfo.beam, caseInfo.path, resultDir);
end

D_mrad = nan(numel(files), 1);
probability = nan(numel(files), 1);
meanFirstWarningTime = nan(numel(files), 1);
nSamples = zeros(numel(files), 1);
nDetected = zeros(numel(files), 1);
nFiniteFirst = zeros(numel(files), 1);
sourceFile = strings(numel(files), 1);

for k = 1:numel(files)
    sourceFile(k) = files(k).name;
    D_mrad(k) = parseDFromFileName(files(k).name);

    data = load(fullfile(files(k).folder, files(k).name), ...
        caseInfo.detectVar, caseInfo.firstVar);
    if ~isfield(data, caseInfo.detectVar) || ~isfield(data, caseInfo.firstVar)
        error('Missing expected variables in %s.', files(k).name);
    end

    detect = data.(caseInfo.detectVar)(:);
    firstTime = data.(caseInfo.firstVar)(:);
    isDetected = detect ~= 0;
    validFirstTime = isDetected & isfinite(firstTime);

    nSamples(k) = numel(detect);
    nDetected(k) = nnz(isDetected);
    nFiniteFirst(k) = nnz(validFirstTime);
    probability(k) = mean(isDetected);

    if any(validFirstTime)
        meanFirstWarningTime(k) = mean(firstTime(validFirstTime));
    end
end

[D_mrad, order] = sort(D_mrad);
probability = probability(order);
meanFirstWarningTime = meanFirstWarningTime(order);
nSamples = nSamples(order);
nDetected = nDetected(order);
nFiniteFirst = nFiniteFirst(order);
sourceFile = sourceFile(order);
end

function D_mrad = parseDFromFileName(fileName)
tokens = regexp(fileName, '_D([0-9]+(?:\.[0-9]+)?)d', 'tokens', 'once');
if isempty(tokens)
    error('Cannot parse fasan_D from file name: %s', fileName);
end
D_mrad = str2double(tokens{1});
end

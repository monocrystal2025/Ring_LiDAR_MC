function output = lubang(varargin)
%LUBANG Plot four normalized annular-beam robustness studies in a 2x2 grid.
%
%   Run the four data functions first:
%       prf_omiga_data
%       energy_budget_data
%       atmosphere_data
%       receiver_noise_data
%
%   LUBANG then loads their MAT outputs from LUBANG_data and plots kappa
%   against four normalized physical directions. Filled markers indicate
%   kappaHigh95 < 1; open markers indicate that robust advantage has not
%   been established at approximately 95% confidence. Confidence is
%   encoded by marker fill rather than long error bars so that the trends
%   remain readable.

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end
dataDir = fullfile(scriptDir, 'LUBANG_data');
defaultFiles = [ ...
    string(fullfile(dataDir, 'prf_omiga_data.mat')), ...
    string(fullfile(dataDir, 'energy_budget_data.mat')), ...
    string(fullfile(dataDir, 'atmosphere_data.mat')), ...
    string(fullfile(dataDir, 'receiver_noise_data.mat'))];

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'DataFiles', defaultFiles, ...
    @(x) (isstring(x) || iscellstr(x)) && numel(x) == 4);
parse(parser, varargin{:});
dataFiles = string(parser.Results.DataFiles);

results = cell(4, 1);
for studyIndex = 1:4
    assert(isfile(dataFiles(studyIndex)), 'lubang:MissingDataFile', ...
        ['Missing data file: %s\nRun the corresponding data function ', ...
        'before plotting.'], dataFiles(studyIndex));
    saved = load(dataFiles(studyIndex), 'result');
    assert(isfield(saved, 'result'), 'lubang:MissingResult', ...
        'Variable result is missing from %s.', dataFiles(studyIndex));
    validateResult(saved.result, dataFiles(studyIndex));
    results{studyIndex} = saved.result;
end

figureHandle = figure('Color', 'w', 'Position', [45, 45, 750, 680], ...
    'Name', 'Normalized annular-beam robustness analysis', ...
    'DefaultAxesFontSize', 14, 'DefaultTextFontSize', 14);
layout = tiledlayout(figureHandle, 2, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
axesHandles = gobjects(2, 2);

for studyIndex = 1:4
    ax = nexttile(layout, studyIndex);
    axesHandles(studyIndex) = ax;
    drawStudy(ax, results{studyIndex});
end

output = struct();
output.figure = figureHandle;
output.axes = axesHandles;
output.dataFiles = dataFiles;
output.results = results;
output.pngFile = fullfile(scriptDir, 'lubang.png');
exportgraphics(figureHandle, output.pngFile, 'Resolution', 300);
end

function drawStudy(ax, result)
plotData = result.plotData;
seriesNames = unique(plotData.series, 'stable');
colors = lines(numel(seriesNames));
markers = {'o', 's', '^', 'd', 'v', '>'};
lineStyles = {'-', '--', '-.', ':'};
hold(ax, 'on');
legendHandles = gobjects(numel(seriesNames), 1);

for seriesIndex = 1:numel(seriesNames)
    member = plotData.series == seriesNames(seriesIndex);
    x = plotData.x(member);
    y = plotData.kappa(member);
    robust = plotData.robustAdvantage(member);
    [x, order] = sort(x);
    y = y(order);
    robust = robust(order);
    valid = isfinite(x) & isfinite(y);
    marker = markers{1 + mod(seriesIndex - 1, numel(markers))};
    lineStyle = lineStyles{1 + mod(seriesIndex - 1, numel(lineStyles))};

    if ~any(valid)
        legendHandles(seriesIndex) = plot(ax, nan, nan, ...
            'LineStyle', lineStyle, 'Marker', marker, ...
            'Color', colors(seriesIndex, :), ...
            'DisplayName', seriesNames(seriesIndex));
        continue;
    end

    % Keep NaN values in the line. MATLAB then breaks the curve at a
    % parameter point where the target probability was not reached,
    % instead of drawing an unsupported bridge across the missing region.
    legendHandles(seriesIndex) = plot(ax, x, y, ...
        'LineStyle', lineStyle, 'Marker', marker, ...
        'Color', colors(seriesIndex, :), 'LineWidth', 3, ...
        'MarkerSize', 6.5, 'MarkerFaceColor', 'w', ...
        'DisplayName', seriesNames(seriesIndex));
    if any(robust & valid)
        filled = robust & valid;
        plot(ax, x(filled), y(filled), marker, ...
            'Color', colors(seriesIndex, :), ...
            'MarkerFaceColor', colors(seriesIndex, :), ...
            'MarkerSize', 6.5, 'LineWidth', 1.2, ...
            'HandleVisibility', 'off');
    end
end

yline(ax, 1, 'k--', '\kappa=1', 'Interpreter', 'tex', ...
    'LabelHorizontalAlignment', 'left', 'LineWidth', 2, ...
    'FontSize', 14);
xline(ax, 1, ':', "\eta=1", 'Interpreter', 'tex', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'right', ...
    'Color', [0.35, 0.35, 0.35], 'LineWidth', 2.0, ...
    'FontSize', 14);
hold(ax, 'off');
grid(ax, 'on');
box(ax, 'on');
xlabel(ax, studyXLabel(result.metadata.studyName), ...
    'Interpreter', 'latex', 'FontSize', 16);
ylabel(ax, '\kappa', 'Interpreter', 'tex', 'FontSize', 14);
legendHandle = legend(ax, legendHandles, ...
    'Location', studyLegendLocation(result.metadata.studyName), ...
    'Box', 'off', ...
    'FontSize', 12);
set(ax, 'FontName', 'Helvetica', 'FontSize', 14, ...
    'LineWidth', 1.5, 'Layer', 'top');
ax.XLabel.FontSize = 16;
applyUsefulXLimits(ax, plotData.x, result.metadata.studyName);
applyUsefulLimits(ax, plotData.kappa);
showUnavailableMarkers(ax, plotData);
moveLegendDown(legendHandle, 0.035);
end

function applyUsefulXLimits(ax, xValues, studyName)
% All four normalized studies use the same two-decade presentation range.
% xValues and studyName remain inputs to keep this helper's call explicit.
assert(any(isfinite(xValues) & xValues > 0), ...
    'lubang:MissingPositiveX', 'No positive normalized x values exist.');
assert(strlength(string(studyName)) > 0, ...
    'lubang:MissingStudyName', 'The study name is empty.');
set(ax, 'XScale', 'log');
xlim(ax, [0.1, 10]);
end

function location = studyLegendLocation(~)
location = 'west';
end

function moveLegendDown(legendHandle, normalizedShift)
legendHandle.Units = 'normalized';
position = legendHandle.Position;
position(2) = max(0.01, position(2) - normalizedShift);
legendHandle.Position = position;
end

function label = studyXLabel(studyName)
switch string(studyName)
    case "scan_sampling"
        label = ['$\eta_S$'];
    case "energy_budget"
        label = ['$\eta_E$'];
    case "atmospheric_quality"
        label = ['$\eta_A$'];
    case "receiver_noise_quality"
        label = ['$\eta_R$'];
    otherwise
        error('lubang:UnknownStudyName', ...
            'No x-axis formula is defined for study %s.', studyName);
end
end

function showUnavailableMarkers(ax, plotData)
% A missing kappa is a calculated operating point at which one or both
% beams never reached the target probability. Show it explicitly without
% assigning a fabricated kappa value or joining it into a trend line.
unavailable = isfinite(plotData.x) & ~isfinite(plotData.kappa);
if ~any(unavailable)
    return;
end
x = unique(plotData.x(unavailable));
limits = ylim(ax);
y = limits(1) + 0.025 * diff(limits);
hold(ax, 'on');
plot(ax, x, repmat(y, size(x)), 'x', ...
    'Color', [0.35, 0.35, 0.35], 'MarkerSize', 7, ...
    'LineWidth', 1.3, 'HandleVisibility', 'off');
text(ax, 0.99, 0.865, '\times: P_D target not reached', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom', 'Interpreter', 'tex', ...
    'Color', [0.10, 0.10, 0.10], 'FontSize', 12);
hold(ax, 'off');
end

function applyUsefulLimits(ax, kappa)
assert(isnumeric(kappa), 'lubang:InvalidKappa', ...
    'Kappa data must be numeric.');
ylim(ax, [0.4, 1.2]);
yticks(ax, 0.4:0.2:1.2);
end

function validateResult(result, dataFile)
required = {'metadata', 'config', 'plotData', 'summary'};
assert(all(isfield(result, required)), 'lubang:InvalidResult', ...
    'The result structure in %s is incomplete.', dataFile);
requiredPlotFields = {'series', 'x', 'kappa', 'kappaLow95', ...
    'kappaHigh95', 'robustAdvantage'};
assert(all(isfield(result.plotData, requiredPlotFields)), ...
    'lubang:InvalidPlotData', ...
    'The plotData structure in %s is incomplete.', dataFile);
assert(abs(result.config.smallWidth_mrad - 1.0) < 1e-12, ...
    'lubang:UnexpectedSmallWidth', ...
    'Expected wd=1.0 mrad in %s.', dataFile);
end

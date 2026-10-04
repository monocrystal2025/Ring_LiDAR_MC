function figureHandles = plot_sensitivity_results(outputRoot)
%PLOT_SENSITIVITY_RESULTS Create the four-panel manuscript sensitivity plot.

arguments
    outputRoot (1,1) string = string(fileparts(mfilename('fullpath')))
end

resultsDir = fullfile(outputRoot, 'results');
figuresDir = fullfile(outputRoot, 'figures');
effectsFile = fullfile(resultsDir, 'sensitivity_screen_effects.csv');
refinedFile = fullfile(resultsDir, 'sensitivity_refined_summary.csv');
selectedFile = fullfile(resultsDir, 'selected_factors.csv');
if ~isfile(effectsFile) || ~isfile(refinedFile) || ~isfile(selectedFile)
    error('plot_sensitivity_results:MissingResults', ...
        'Run run_sensitivity_analysis in full or smoke mode before plotting.');
end
if ~isfolder(figuresDir)
    mkdir(figuresDir);
end

effects = readtable(effectsFile, TextType='string');
refined = readtable(refinedFile, TextType='string');
selected = readtable(selectedFile, TextType='string');
probabilityFactor = selected.ProbabilityFactor(1);
warningFactor = selected.WarningTimeFactor(1);

fig = figure('Color', 'w', 'Position', [80, 60, 1450, 900]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(layout, 1);
plotEffectBars(ax1, effects, 'ProbabilityEffect_pp');
ylabel(ax1, 'Change in \DeltaP_{annular-line} (percentage points)', ...
    'Interpreter', 'tex');
title(ax1, '(a) Sensitivity of annular-line detection advantage');

ax2 = nexttile(layout, 2);
plotEffectBars(ax2, effects, 'WarningEffect_s');
ylabel(ax2, 'Change in \Deltat_{line-annular} (s)', 'Interpreter', 'tex');
title(ax2, '(b) Sensitivity of annular-line warning-time advantage');

ax3 = nexttile(layout, 3);
plotRefinedCurves(ax3, refined, probabilityFactor, 'probability');
title(ax3, "(c) Detection response: " + factorDisplayName(effects, probabilityFactor));

ax4 = nexttile(layout, 4);
plotRefinedCurves(ax4, refined, warningFactor, 'warning');
title(ax4, "(d) Warning-time response: " + factorDisplayName(effects, warningFactor));

title(layout, 'Paired Monte Carlo sensitivity analysis of spot, line, and annular beams', ...
    'FontWeight', 'bold');

pngFile = fullfile(figuresDir, 'sensitivity_analysis_main.png');
figFile = fullfile(figuresDir, 'sensitivity_analysis_main.fig');
exportgraphics(fig, pngFile, 'Resolution', 300);
savefig(fig, figFile);
figureHandles = fig;
end

function plotEffectBars(ax, effects, variableName)
factorNames = unique(effects.Factor, 'stable');
trajectories = ["RA", "SP"];
values = nan(numel(factorNames), numel(trajectories));
labels = strings(numel(factorNames), 1);
for factorIdx = 1:numel(factorNames)
    rowsForFactor = effects.Factor == factorNames(factorIdx);
    labels(factorIdx) = effects.FactorLabel(find(rowsForFactor, 1, 'first'));
    for trajectoryIdx = 1:numel(trajectories)
        row = rowsForFactor & effects.Trajectory == trajectories(trajectoryIdx);
        if any(row)
            values(factorIdx, trajectoryIdx) = effects.(variableName)(row);
        end
    end
end

bars = bar(ax, values, 'grouped');
bars(1).FaceColor = [0.20, 0.55, 0.80];
bars(2).FaceColor = [0.90, 0.45, 0.15];
hold(ax, 'on');
if strcmp(variableName, 'ProbabilityEffect_pp')
    lowName = 'ProbabilityEffect_CI_Low_pp';
    highName = 'ProbabilityEffect_CI_High_pp';
else
    lowName = 'WarningEffect_CI_Low_s';
    highName = 'WarningEffect_CI_High_s';
end
if all(ismember({lowName, highName}, effects.Properties.VariableNames))
    lowValues = nan(size(values));
    highValues = nan(size(values));
    for factorIdx = 1:numel(factorNames)
        rowsForFactor = effects.Factor == factorNames(factorIdx);
        for trajectoryIdx = 1:numel(trajectories)
            row = rowsForFactor & effects.Trajectory == trajectories(trajectoryIdx);
            if any(row)
                lowValues(factorIdx, trajectoryIdx) = effects.(lowName)(row);
                highValues(factorIdx, trajectoryIdx) = effects.(highName)(row);
            end
        end
    end
    drawnow;
    for trajectoryIdx = 1:numel(trajectories)
        errorbar(ax, bars(trajectoryIdx).XEndPoints, ...
            values(:, trajectoryIdx), ...
            values(:, trajectoryIdx) - lowValues(:, trajectoryIdx), ...
            highValues(:, trajectoryIdx) - values(:, trajectoryIdx), ...
            'k', 'LineStyle', 'none', 'LineWidth', 0.8, ...
            'CapSize', 4, 'HandleVisibility', 'off');
    end
end
yline(ax, 0, 'k-', 'LineWidth', 0.8);
grid(ax, 'on');
box(ax, 'on');
ax.XTick = 1:numel(labels);
ax.XTickLabel = labels;
ax.XTickLabelRotation = 35;
legend(ax, trajectories, 'Location', 'best', 'Box', 'off');
end

function plotRefinedCurves(ax, refined, factorName, metric)
rows = refined.Factor == factorName;
data = refined(rows, :);
if isempty(data)
    error('plot_sensitivity_results:MissingRefinedFactor', ...
        'No refined results found for factor %s.', factorName);
end

beamNames = ["point", "line", "ring"];
beamLabels = ["Spot", "Line", "Annular"];
beamColors = [0.35, 0.62, 0.18; 0.00, 0.45, 0.74; 0.85, 0.33, 0.10];
trajectories = ["RA", "SP"];
lineStyles = ["-", "--"];
markers = ["^", "o", "s"];
hold(ax, 'on');

for beamIdx = 1:numel(beamNames)
    for trajectoryIdx = 1:numel(trajectories)
        curve = data(data.Beam == beamNames(beamIdx) & ...
            data.Trajectory == trajectories(trajectoryIdx), :);
        curve = sortrows(curve, 'LevelValue');
        if metric == "probability"
            y = 100 .* curve.PDetect;
            lowError = 100 .* (curve.PDetect - curve.PDetect_CI_Low);
            highError = 100 .* (curve.PDetect_CI_High - curve.PDetect);
        else
            y = curve.MeanWarning_s;
            lowError = curve.MeanWarning_s - curve.MeanWarning_CI_Low_s;
            highError = curve.MeanWarning_CI_High_s - curve.MeanWarning_s;
        end
        displayName = beamLabels(beamIdx) + "-" + trajectories(trajectoryIdx);
        errorbar(ax, curve.LevelValue, y, lowError, highError, ...
            'Color', beamColors(beamIdx, :), ...
            'LineStyle', lineStyles(trajectoryIdx), ...
            'Marker', markers(beamIdx), 'MarkerFaceColor', 'w', ...
            'LineWidth', 1.5, 'MarkerSize', 5, ...
            'DisplayName', displayName);
    end
end

grid(ax, 'on');
box(ax, 'on');
xlabel(ax, factorDisplayNameFromRefined(data));
if metric == "probability"
    ylabel(ax, 'P_{detect} (%)', 'Interpreter', 'tex');
else
    ylabel(ax, 'Mean first-warning time (s)');
end
legend(ax, 'Location', 'best', 'NumColumns', 2, 'Box', 'off');
end

function label = factorDisplayName(effects, factorName)
row = find(effects.Factor == factorName, 1, 'first');
label = effects.FactorLabel(row);
end

function label = factorDisplayNameFromRefined(data)
label = data.FactorLabel(1);
unit = string(data.Unit(1));
if strlength(unit) > 0 && unit ~= "1"
    label = label + " (" + unit + ")";
end
end

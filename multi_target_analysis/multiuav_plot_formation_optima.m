function figureFiles = multiuav_plot_formation_optima( ...
    summary, optima, optimumComparisons, figuresDir, threshold)
%MULTIUAV_PLOT_FORMATION_OPTIMA Plot full manuscript-grid optimum search.

arguments
    summary table
    optima table
    optimumComparisons table
    figuresDir (1,1) string
    threshold (1,1) double {mustBePositive}
end

if ~isfolder(figuresDir)
    mkdir(figuresDir);
end
paths = unique(summary.Trajectory, 'stable');
targetCounts = unique(summary.TargetCount, 'stable');
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1160, 760]);
layout = tiledlayout(numel(targetCounts), numel(paths), ...
    'TileSpacing', 'compact', 'Padding', 'compact');
styles = {'-', '-'};
colors = [0.10, 0.35, 0.75; 0.85, 0.25, 0.15];
beamNames = ["line", "ring"];

for targetIndex = 1:numel(targetCounts)
    for pathIndex = 1:numel(paths)
        axesHandle = nexttile(layout);
        hold(axesHandle, 'on');
        for beamIndex = 1:numel(beamNames)
            select = summary.Trajectory == paths(pathIndex) & ...
                summary.TargetCount == targetCounts(targetIndex) & ...
                summary.Scenario == "sync_formation" & ...
                summary.Threshold_s == threshold & ...
                summary.Beam == beamNames(beamIndex);
            curve = sortrows(summary(select, :), 'wD_mrad');
            [curveX, curveY] = breakLargeGaps(curve.wD_mrad, curve.PAll);
            plot(axesHandle, curveX, curveY, ...
                styles{beamIndex}, 'Color', colors(beamIndex, :), ...
                'LineWidth', 1.6, 'DisplayName', beamNames(beamIndex));
            optimum = optima.Trajectory == paths(pathIndex) & ...
                optima.TargetCount == targetCounts(targetIndex) & ...
                optima.Beam == beamNames(beamIndex);
            scatter(axesHandle, optima.BestWD_mrad(optimum), ...
                optima.MaxPAll(optimum), 65, colors(beamIndex, :), ...
                'filled', 'HandleVisibility', 'off');
        end
        grid(axesHandle, 'on');
        xlim(axesHandle, [min(summary.wD_mrad), max(summary.wD_mrad)]);
        ylim(axesHandle, [0, 1]);
        xlabel(axesHandle, 'w_D (mrad)');
        ylabel(axesHandle, sprintf('P(all detected by %g s)', threshold));
        title(axesHandle, sprintf('%s / M=%d', ...
            paths(pathIndex), targetCounts(targetIndex)));
        legend(axesHandle, 'Location', 'best');
    end
end
title(layout, 'Synchronized-formation w_D optimization');
curveFile = fullfile(figuresDir, 'formationopt_full_wD_curves.png');
exportgraphics(figureHandle, curveFile, 'Resolution', 180);
close(figureHandle);

comparison = sortrows(optimumComparisons, {'TargetCount', 'Trajectory'});
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 900, 520]);
x = 1:height(comparison);
errorbar(x, comparison.DeltaCVPAll_pp, ...
    comparison.DeltaCVPAll_pp - comparison.DeltaCVPAllCILow_pp, ...
    comparison.DeltaCVPAllCIHigh_pp - comparison.DeltaCVPAll_pp, ...
    'o', 'LineWidth', 1.5, 'MarkerFaceColor', [0.85, 0.25, 0.15]);
yline(0, 'k--', 'HandleVisibility', 'off');
grid on;
xticks(x);
xticklabels(comparison.Trajectory + ", M=" + comparison.TargetCount);
xlabel('Trajectory and target count');
ylabel('Cross-validated \DeltaP_{all} (ring - line), percentage points');
title('Comparison after independently re-optimizing line and ring w_D');
deltaFile = fullfile(figuresDir, 'formationopt_reoptimized_delta.png');
exportgraphics(figureHandle, deltaFile, 'Resolution', 180);
close(figureHandle);

figureFiles = [string(curveFile); string(deltaFile)];
end

function [xPlot, yPlot] = breakLargeGaps(x, y)
x = x(:);
y = y(:);
if numel(x) < 3
    xPlot = x;
    yPlot = y;
    return;
end
typicalStep = median(diff(x));
gaps = find(diff(x) > 1.5 * typicalStep);
xPlot = x;
yPlot = y;
for gapIndex = numel(gaps):-1:1
    insertion = gaps(gapIndex) + 1;
    xPlot = [xPlot(1:insertion-1); NaN; xPlot(insertion:end)]; %#ok<AGROW>
    yPlot = [yPlot(1:insertion-1); NaN; yPlot(insertion:end)]; %#ok<AGROW>
end
end

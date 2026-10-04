function figureFiles = multiuav_plot_results(summary, comparisons, figuresDir)
%MULTIUAV_PLOT_RESULTS Create compact preview figures without UI display.

arguments
    summary table
    comparisons table
    figuresDir (1,1) string
end

if ~isfolder(figuresDir)
    mkdir(figuresDir);
end
figureFiles = strings(0, 1);
if isempty(summary)
    return;
end

threshold = min(5, max(summary.Threshold_s));
paths = unique(summary.Trajectory, 'stable');
scenarios = unique(summary.Scenario, 'stable');
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1120, 760]);
layout = tiledlayout(numel(scenarios), numel(paths), ...
    'TileSpacing', 'compact', 'Padding', 'compact');
for scenarioIndex = 1:numel(scenarios)
    for pathIndex = 1:numel(paths)
        axesHandle = nexttile(layout);
        hold(axesHandle, 'on');
        select = summary.Scenario == scenarios(scenarioIndex) & ...
            summary.Trajectory == paths(pathIndex) & ...
            summary.Threshold_s == threshold & ...
            ((summary.Beam == "line" & summary.wD_mrad == 195) | ...
            (summary.Beam == "ring" & summary.wD_mrad == 125));
        subset = summary(select, :);
        beamNames = ["line", "ring"];
        styles = {'-s', '-o'};
        colors = [0.10, 0.35, 0.75; 0.85, 0.25, 0.15];
        for beamIndex = 1:numel(beamNames)
            rows = subset.Beam == beamNames(beamIndex);
            beamTable = sortrows(subset(rows, :), 'TargetCount');
            if isempty(beamTable)
                continue;
            end
            plot(axesHandle, beamTable.TargetCount, beamTable.PAll, ...
                styles{beamIndex}, 'Color', colors(beamIndex, :), ...
                'LineWidth', 1.6, 'MarkerFaceColor', colors(beamIndex, :), ...
                'DisplayName', sprintf('%s, w_D=%g mrad', ...
                beamNames(beamIndex), beamTable.wD_mrad(1)));
        end
        grid(axesHandle, 'on');
        ylim(axesHandle, [0, 1]);
        xlabel(axesHandle, 'Target count M');
        ylabel(axesHandle, sprintf('P(all detected by %g s)', threshold));
        title(axesHandle, sprintf('%s / %s', ...
            paths(pathIndex), replace(scenarios(scenarioIndex), '_', ' ')));
        legend(axesHandle, 'Location', 'best');
    end
end
title(layout, 'Multi-UAV joint capture preview');
fileName = fullfile(figuresDir, 'preview_joint_capture.png');
exportgraphics(figureHandle, fileName, 'Resolution', 180);
close(figureHandle);
figureFiles(end + 1, 1) = string(fileName);

if ~isempty(comparisons) && ismember('DeltaPAll_pp', comparisons.Properties.VariableNames)
    select = comparisons.Comparison == "SingleTargetReference" & ...
        comparisons.Threshold_s == threshold;
    delta = comparisons(select, :);
    if ~isempty(delta)
        figureHandle = figure('Visible', 'off', 'Color', 'w', ...
            'Position', [100, 100, 900, 520]);
        group = delta.Trajectory + "/" + ...
            replace(delta.Scenario, '_', ' ');
        groupNames = unique(group, 'stable');
        hold on;
        for groupIndex = 1:numel(groupNames)
            groupTable = sortrows(delta(group == groupNames(groupIndex), :), ...
                'TargetCount');
            plot(groupTable.TargetCount, groupTable.DeltaPAll_pp, '-o', ...
                'LineWidth', 1.4, 'MarkerFaceColor', 'auto', ...
                'DisplayName', groupNames(groupIndex));
        end
        yline(0, 'k--', 'HandleVisibility', 'off');
        grid on;
        xlabel('Target count M');
        ylabel('\DeltaP_{all} (ring - line), percentage points');
        title(sprintf('Annular advantage by %g s', threshold));
        legend('Location', 'best');
        fileName = fullfile(figuresDir, 'preview_ring_minus_line.png');
        exportgraphics(figureHandle, fileName, 'Resolution', 180);
        close(figureHandle);
        figureFiles(end + 1, 1) = string(fileName);
    end
end
end

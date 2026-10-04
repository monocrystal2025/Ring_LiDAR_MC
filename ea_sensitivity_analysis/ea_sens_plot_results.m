function handles = ea_sens_plot_results(outputRoot)
%EA_SENS_PLOT_RESULTS Create the paper-oriented ROE sensitivity figure.
%
% Panel (a) is a tornado chart of endpoint effects on kappa. Panel (b)
% displays the two-factor interaction surface when available; before a full
% run it falls back to baseline time-limited detection curves.

arguments
    outputRoot (1,1) string = ""
end

if strlength(outputRoot) == 0
    outputRoot = string(fileparts(mfilename("fullpath")));
end
resultsDir = fullfile(outputRoot, "results");
figuresDir = fullfile(outputRoot, "figures");
if ~isfolder(figuresDir)
    mkdir(figuresDir);
end

configFile = fullfile(resultsDir, "ea_sensitivity_screen_results.mat");
effectsFile = fullfile(resultsDir, "ea_sensitivity_screen_effects.csv");
timelyFile = fullfile(resultsDir, "ea_sensitivity_screen_timely.csv");
interactionFile = fullfile(resultsDir, ...
    "ea_sensitivity_interaction_iso_performance.csv");
if ~isfile(configFile) || ~isfile(effectsFile)
    error("ea_sens_plot_results:MissingResults", ...
        "Run screen or full sensitivity analysis before plotting.");
end

saved = load(configFile, "cfg");
cfg = saved.cfg;
effects = readtable(effectsFile, TextType="string");

fig = figure("Color", "w", "Position", [100, 100, 1180, 470]);
layout = tiledlayout(fig, 1, 2, "TileSpacing", "compact", ...
    "Padding", "compact");

ax1 = nexttile(layout);
plotTornado(ax1, effects);

ax2 = nexttile(layout);
if isfile(interactionFile)
    interaction = readtable(interactionFile, TextType="string");
    plotInteraction(ax2, interaction, cfg);
else
    timely = readtable(timelyFile, TextType="string");
    plotBaselineCurves(ax2, timely, cfg);
end

title(layout, sprintf( ...
    "ROE sensitivity at P_D(%.0f s) = %.2f", ...
    cfg.primaryTimelyThreshold_s, cfg.primaryProbability), ...
    "FontName", "Helvetica", "FontSize", 14);

outPng = fullfile(figuresDir, "ea_sensitivity_main.png");
outFig = fullfile(figuresDir, "ea_sensitivity_main.fig");
exportgraphics(fig, outPng, "Resolution", 300);
savefig(fig, outFig);

handles = struct("mainFigure", fig, "png", string(outPng), ...
    "fig", string(outFig));
end

function plotTornado(ax, effects)
valid = effects(isfinite(effects.Effect), :);
if isempty(valid)
    axis(ax, "off");
    text(ax, 0.5, 0.5, "No finite endpoint effects", ...
        "HorizontalAlignment", "center");
    return;
end
[~, order] = sort(valid.AbsEffect, "ascend");
valid = valid(order, :);
y = 1:height(valid);
colors = repmat([0.20, 0.48, 0.72], height(valid), 1);
colors(valid.Effect > 0, :) = repmat([0.85, 0.33, 0.10], ...
    nnz(valid.Effect > 0), 1);

bars = barh(ax, y, valid.Effect, 0.72, "FaceColor", "flat");
bars.CData = colors;
xline(ax, 0, "k-", "LineWidth", 1.2);
yticks(ax, y);
yticklabels(ax, valid.FactorLabel);
xlabel(ax, "\Delta\kappa (high - low)", "Interpreter", "tex");
title(ax, "(a) Endpoint sensitivity");
grid(ax, "on");
box(ax, "on");
set(ax, "FontName", "Helvetica", "FontSize", 11, "LineWidth", 1.2);
end

function plotInteraction(ax, interaction, cfg)
rows = interaction(abs(interaction.TimeThreshold_s - ...
    cfg.primaryTimelyThreshold_s) < 1e-12 & ...
    abs(interaction.TargetProbability - ...
    cfg.primaryProbability) < 1e-12, :);
if isempty(rows)
    axis(ax, "off");
    text(ax, 0.5, 0.5, "Primary interaction result unavailable", ...
        "HorizontalAlignment", "center");
    return;
end

x = unique(rows.Factor1Value, "sorted");
y = unique(rows.Factor2Value, "sorted");
z = nan(numel(y), numel(x));
for rowIdx = 1:height(rows)
    xIdx = find(x == rows.Factor1Value(rowIdx), 1);
    yIdx = find(y == rows.Factor2Value(rowIdx), 1);
    z(yIdx, xIdx) = rows.KappaRingLine(rowIdx);
end

[xDisplay, xUnit] = displayFactorValues( ...
    x, rows.Factor1(1), rows.Factor1Unit(1));
[yDisplay, yUnit] = displayFactorValues( ...
    y, rows.Factor2(1), rows.Factor2Unit(1));
imageHandle = imagesc(ax, 1:numel(x), 1:numel(y), z);
imageHandle.AlphaData = isfinite(z);
ax.Color = [0.86, 0.86, 0.86];
set(ax, "YDir", "normal");
xticks(ax, 1:numel(x));
xticklabels(ax, compose("%.3g", xDisplay));
yticks(ax, 1:numel(y));
yticklabels(ax, compose("%.3g", yDisplay));
hold(ax, "on");
if all(isfinite(z), "all") && min(z, [], "all") <= 1 && ...
        max(z, [], "all") >= 1
    contour(ax, 1:numel(x), 1:numel(y), z, [1, 1], ...
        "w-", "LineWidth", 2);
end
for yIdx = 1:numel(y)
    for xIdx = 1:numel(x)
        if isfinite(z(yIdx, xIdx))
            label = sprintf("%.2f", z(yIdx, xIdx));
        else
            label = "N/A";
        end
        text(ax, xIdx, yIdx, label, ...
            "HorizontalAlignment", "center", ...
            "FontName", "Helvetica", "FontSize", 10, ...
            "Color", "k");
    end
end
hold(ax, "off");
colormap(ax, turbo);
finiteValues = z(isfinite(z));
if ~isempty(finiteValues)
    valueMin = min(finiteValues);
    valueMax = max(finiteValues);
    if abs(valueMax - valueMin) < 1e-12
        padding = max(0.05, 0.1 * abs(valueMin));
        clim(ax, [max(0, valueMin - padding), valueMax + padding]);
    else
        clim(ax, [valueMin, valueMax]);
    end
end
cb = colorbar(ax);
cb.Label.String = "\kappa = w_{D,annular}/w_{D,line}";
xlabel(ax, axisLabel(rows.Factor1Label(1), xUnit));
ylabel(ax, axisLabel(rows.Factor2Label(1), yUnit));
title(ax, "(b) Dominant-factor interaction");
box(ax, "on");
set(ax, "FontName", "Helvetica", "FontSize", 11, "LineWidth", 1.2);
end

function plotBaselineCurves(ax, timely, cfg)
rows = timely(timely.CaseId == "baseline" & ...
    abs(timely.TimeThreshold_s - ...
    cfg.primaryTimelyThreshold_s) < 1e-12, :);
hold(ax, "on");
beamOrder = ["line", "ring"];
styles = ["-o", "-s"];
labels = ["Line", "Annular"];
for beamIdx = 1:numel(beamOrder)
    beamRows = rows(strcmpi(rows.Beam, beamOrder(beamIdx)), :);
    beamRows = sortrows(beamRows, "WD_mrad");
    plot(ax, beamRows.WD_mrad, beamRows.PTimely, styles(beamIdx), ...
        "LineWidth", 1.8, "MarkerSize", 5, ...
        "DisplayName", labels(beamIdx));
end
yline(ax, cfg.primaryProbability, "k--", "Target");
hold(ax, "off");
xlabel(ax, "w_D (mrad)");
ylabel(ax, sprintf("P_D(%.0f s)", cfg.primaryTimelyThreshold_s));
title(ax, "(b) Baseline time-limited detection");
legend(ax, "Location", "best", "Box", "off");
grid(ax, "on");
box(ax, "on");
ylim(ax, [0, 1]);
set(ax, "FontName", "Helvetica", "FontSize", 11, "LineWidth", 1.2);
end

function [displayValues, displayUnit] = displayFactorValues(values, factor, unit)
displayValues = values;
displayUnit = unit;
if factor == "beamWidth" && unit == "rad"
    displayValues = 1e3 .* values;
    displayUnit = "mrad";
end
end

function label = axisLabel(name, unit)
if strlength(unit) == 0 || unit == "1" || unit == "index"
    label = name;
else
    label = name + " (" + unit + ")";
end
end

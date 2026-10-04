function figures = ea_tornado_plot(outputDirectory, visible)
%EA_TORNADO_PLOT Draw and export the publication tornado chart.

arguments
    outputDirectory (1,1) string
    visible (1,1) string {mustBeMember(visible, ["on", "off"])} = "off"
end

analysisFile = fullfile(outputDirectory, "ea_tornado_analysis.mat");
assert(isfile(analysisFile), "EATornado:MissingAnalysis", ...
    "Run analysis before plotting: %s", analysisFile);
saved = load(analysisFile, "analysis");
analysis = saved.analysis;
effects = analysis.effects;
levels = analysis.levels;
cfg = analysis.design.config;

figureDirectory = fullfile(outputDirectory, "figures");
if ~isfolder(figureDirectory)
    mkdir(figureDirectory);
end

fig = figure("Color", "w", "Units", "centimeters", ...
    "Position", [2, 2, 19, 12], "Visible", visible);
ax = axes(fig, "Position", [0.25, 0.19, 0.69, 0.69]);
hold(ax, "on");
nFactors = height(effects);
y = (1:nFactors).';

values = [effects.TornadoLow; effects.TornadoHigh; ...
    effects.LowKappaCI_Low; effects.LowKappaCI_High; ...
    effects.HighKappaCI_Low; effects.HighKappaCI_High; ...
    analysis.baselineKappaCI(:); analysis.baselineKappa; 1];
values = values(isfinite(values));
if isempty(values)
    xLimits = [0.7, 1.05];
else
    xMinimum = min(values);
    xMaximum = max(values);
    span = max(xMaximum - xMinimum, 0.08);
    xLimits = [xMinimum - 0.16 * span, xMaximum + 0.20 * span];
end
xlim(ax, xLimits);

baselineCI = analysis.baselineKappaCI;
if all(isfinite(baselineCI))
    patch(ax, [baselineCI(1), baselineCI(2), baselineCI(2), baselineCI(1)], ...
        [0.45, 0.45, nFactors + 0.55, nFactors + 0.55], ...
        [0.88, 0.88, 0.88], "EdgeColor", "none", ...
        "FaceAlpha", 0.55, "HandleVisibility", "off");
end
xline(ax, 1, "--", "Color", [0.35, 0.35, 0.35], ...
    "LineWidth", 1.0, "DisplayName", "\kappa = 1");
if isfinite(analysis.baselineKappa)
    xline(ax, analysis.baselineKappa, "k-", "LineWidth", 1.2, ...
        "DisplayName", "Baseline \kappa");
end

lowColor = [0.12, 0.47, 0.71];
highColor = [0.89, 0.36, 0.12];
barColor = [0.76, 0.79, 0.82];
middleColor = [0.33, 0.33, 0.33];
for rowIndex = 1:nFactors
    row = effects(rowIndex, :);
    if all(isfinite([row.TornadoLow, row.TornadoHigh]))
        patch(ax, [row.TornadoLow, row.TornadoHigh, ...
            row.TornadoHigh, row.TornadoLow], ...
            [y(rowIndex) - 0.22, y(rowIndex) - 0.22, ...
            y(rowIndex) + 0.22, y(rowIndex) + 0.22], ...
            barColor, "EdgeColor", [0.58, 0.61, 0.64], ...
            "LineWidth", 0.7, "HandleVisibility", "off");
    else
        text(ax, mean(xLimits), y(rowIndex), "NA / unattained", ...
            "HorizontalAlignment", "center", "Color", [0.45, 0.45, 0.45], ...
            "FontSize", 8);
    end

    drawEndpoint(ax, row.LowKappa, row.LowKappaCI_Low, ...
        row.LowKappaCI_High, y(rowIndex) - 0.11, lowColor, "o");
    drawEndpoint(ax, row.HighKappa, row.HighKappaCI_Low, ...
        row.HighKappaCI_High, y(rowIndex) + 0.11, highColor, "d");

    middle = levels(levels.Factor == row.Factor & ...
        (abs(levels.Scale - 0.9) < 1e-12 | ...
        abs(levels.Scale - 1.1) < 1e-12), :);
    middleY = y(rowIndex) + [-0.06; 0.06];
    finiteMiddle = isfinite(middle.Kappa);
    scatter(ax, middle.Kappa(finiteMiddle), middleY(finiteMiddle), ...
        18, middleColor, "+", "LineWidth", 0.8, ...
        "HandleVisibility", "off");

    annotateEndpoint(ax, row.LowKappa, y(rowIndex) - 0.11, ...
        formatValue(row.LowValue, row.Unit), lowColor, "right");
    annotateEndpoint(ax, row.HighKappa, y(rowIndex) + 0.11, ...
        formatValue(row.HighValue, row.Unit), highColor, "left");
end

lowLegend = scatter(ax, nan, nan, 30, lowColor, "o", "filled", ...
    "DisplayName", "-20% endpoint");
highLegend = scatter(ax, nan, nan, 30, highColor, "d", "filled", ...
    "DisplayName", "+20% endpoint");
middleLegend = scatter(ax, nan, nan, 22, middleColor, "+", ...
    "LineWidth", 0.9, "DisplayName", "\pm10% audit levels");

labels = effects.FactorLabel;
labels(effects.IntermediateOutsideEndpointRange) = ...
    labels(effects.IntermediateOutsideEndpointRange) + " *";
yticks(ax, y);
yticklabels(ax, labels);
set(ax, "YDir", "reverse", "YLim", [0.45, nFactors + 0.55], ...
    "XLim", xLimits, "FontName", "Arial", "FontSize", 8.5, ...
    "LineWidth", 0.8, "TickDir", "out", "Box", "off", ...
    "XGrid", "on", "GridAlpha", 0.14, "Layer", "top");
xlabel(ax, "Angular-width ratio  \kappa = w_{D,annular}^{min} / w_{D,line}^{min}", ...
    "Interpreter", "tex");
ylabel(ax, "Parameter (ranked by endpoint/base span)");
title(ax, "EA-mode local parameter sensitivity", ...
    "FontWeight", "normal", "FontSize", 11);
subtitle(ax, sprintf( ...
    "Equal requirement: P_D(%.1f s) >= %.0f%%; smaller kappa means a stronger annular advantage", ...
    cfg.deadline_s, 100 * cfg.requiredProbability), ...
    "FontSize", 8.5);
legend(ax, [lowLegend, highLegend, middleLegend], ...
    "Location", "southoutside", "Orientation", "horizontal", ...
    "Box", "off", "FontSize", 7.5);
note = "Five OAT levels per factor; bars use -20%, baseline and +20%. " + ...
    "Whiskers are paired pointwise 95% intervals. N = " + cfg.N + ...
    ", bootstrap = " + analysis.options.Replicates + ...
    ". * intermediate level falls outside the endpoint/base range.";
annotation(fig, "textbox", [0.25, 0.025, 0.70, 0.055], ...
    "String", note, ...
    "EdgeColor", "none", "HorizontalAlignment", "center", ...
    "FontName", "Arial", "FontSize", 7.2);
hold(ax, "off");
drawnow;

baseFile = fullfile(figureDirectory, "ea_tornado_sensitivity");
pngFile = baseFile + ".png";
pdfFile = baseFile + ".pdf";
figFile = baseFile + ".fig";
exportgraphics(fig, pngFile, "Resolution", 400);
exportgraphics(fig, pdfFile, "ContentType", "vector");
savefig(fig, figFile);
figures = struct("png", pngFile, "pdf", pdfFile, "fig", figFile);
if visible == "off"
    close(fig);
end
end

function drawEndpoint(ax, value, low, high, y, color, marker)
if ~isfinite(value)
    return;
end
if all(isfinite([low, high]))
    line(ax, [low, high], [y, y], "Color", color, ...
        "LineWidth", 1.1, "HandleVisibility", "off");
    line(ax, [low, low], y + [-0.045, 0.045], "Color", color, ...
        "LineWidth", 0.8, "HandleVisibility", "off");
    line(ax, [high, high], y + [-0.045, 0.045], "Color", color, ...
        "LineWidth", 0.8, "HandleVisibility", "off");
end
plot(ax, value, y, marker, "Color", color, "MarkerFaceColor", color, ...
    "MarkerSize", 5, "HandleVisibility", "off");
end

function annotateEndpoint(ax, value, y, label, color, side)
if ~isfinite(value)
    return;
end
if side == "right"
    alignment = "right";
    offset = -0.006 * diff(xlim(ax));
else
    alignment = "left";
    offset = 0.006 * diff(xlim(ax));
end
text(ax, value + offset, y, label, "HorizontalAlignment", alignment, ...
    "VerticalAlignment", "middle", "Color", color, "FontSize", 7);
end

function label = formatValue(value, unit)
if unit == "uJ"
    label = sprintf("%g \\muJ", value);
elseif unit == "1e-5 /m"
    label = sprintf("%g \\times10^{-5}/m", value);
elseif unit == "m/s"
    label = sprintf("%g m/s", value);
else
    label = sprintf("%g", value);
end
end

%% SCANTIME_3BEAM Three-panel scan-time comparison for a two-column figure.
% This script reproduces the scan-time calculations in ROIScanTime_1par.m
% and EAScanTime_1par.m without modifying either original file.  It creates
% three aligned heat maps for a full-width Applied Optics-style figure.

close all;

% Figure-wide settings (dimensions include the shared colorbar).
figureWidthCm = 17.8;
figureHeightCm = 6.8;
fontName = 'Times New Roman';
fontSize = 8;
colorMapName = 'turbo';
colorLimits = [0.1, 200];
contourLevels = [2, 5, 10];
contourStyles = {'--', '-', ':'};

% Common optical parameters.
rangeM = 1000;
angularRate = 2 * pi;
beamDivergenceRad = 1e-3:0.2e-3:500e-3;
beamDivergenceMrad = beamDivergenceRad * 1e3;

% Reproduce the two ROI scan-time maps in ROIScanTime_1par.m.
thetaHalfDeg = 5:0.2:45;
thetaDeg = 2 * thetaHalfDeg;
roiTimes = computeRoiTimes(thetaHalfDeg, beamDivergenceRad, rangeM, angularRate);

% Reproduce the EA scan-time map in EAScanTime_1par.m.
phiRad = linspace(pi / 2, 0, 180);
eaTimes = computeEaTimes(phiRad, beamDivergenceRad, rangeM, angularRate);

% Create one full-width, three-panel figure.
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, figureWidthCm, figureHeightCm]);
layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax(1) = nexttile(layout);
drawScanMap(ax(1), beamDivergenceMrad, thetaDeg, roiTimes.discrete, ...
    colorMapName, colorLimits, contourLevels, contourStyles);
ylabel(ax(1), '\theta (deg)');
title(ax(1), '(a) ROI scan: discrete model', 'FontWeight', 'normal');

ax(2) = nexttile(layout);
drawScanMap(ax(2), beamDivergenceMrad, thetaDeg, roiTimes.continuous, ...
    colorMapName, colorLimits, contourLevels, contourStyles);
ylabel(ax(2), '\theta (deg)');
title(ax(2), '(b) ROI scan: continuous model', 'FontWeight', 'normal');

ax(3) = nexttile(layout);
drawScanMap(ax(3), beamDivergenceMrad, phiRad, eaTimes, ...
    colorMapName, colorLimits, contourLevels, contourStyles);
set(ax(3), 'YDir', 'normal', 'YTick', [0, pi / 6, pi / 3, pi / 2], ...
    'YTickLabel', {'0', '\pi/6', '\pi/3', '\pi/2'});
ylabel(ax(3), '\phi (rad)');
title(ax(3), '(c) EA scan', 'FontWeight', 'normal');

% A single colorbar keeps the time scale directly comparable across panels.
cb = colorbar(ax(3), 'southoutside');
cb.Layout.Tile = 'south';
cb.Label.String = 'Scan time (s)';
cb.Ticks = [0.1, 1, 10, 100];
cb.FontName = fontName;
cb.FontSize = fontSize;

% Export both a vector publication file and a high-resolution raster preview.
outputBase = fullfile(pwd, 'scantime_3beam');
exportgraphics(fig, outputBase + ".pdf", 'ContentType', 'vector');
exportgraphics(fig, outputBase + ".png", 'Resolution', 600);


function roiTimes = computeRoiTimes(thetaHalfDeg, divergenceRad, rangeM, angularRate)
%COMPUTEROITIMES Calculate discrete and continuous ROI scan-time maps.

    halfWidthM = rangeM * tand(thetaHalfDeg(:));
    beamDiameterM = 2 * rangeM * tan(divergenceRad(:).' / 2);
    halfWidthGrid = repmat(halfWidthM, 1, numel(beamDiameterM));
    beamDiameterGrid = repmat(beamDiameterM, numel(halfWidthM), 1);

    upperBound = 2 * pi * (halfWidthGrid + 0.5 * beamDiameterGrid) ./ beamDiameterGrid;
    integralTerm = upperBound .* sqrt(upperBound .^ 2 + 1) + ...
        log(abs(upperBound + sqrt(upperBound .^ 2 + 1)));
    continuousLengthM = beamDiameterGrid .* integralTerm / (4 * pi);

    fullCoverage = halfWidthGrid <= beamDiameterGrid;
    continuousLengthM(fullCoverage) = pi * beamDiameterGrid(fullCoverage);

    numberOfPasses = ceil(2 * halfWidthGrid ./ beamDiameterGrid);
    discreteLengthM = 2 * halfWidthGrid .* numberOfPasses + ...
        beamDiameterGrid .* (numberOfPasses - 1);
    singlePass = 2 * halfWidthGrid <= beamDiameterGrid;
    discreteLengthM(singlePass) = 2 * halfWidthGrid(singlePass);

    roiTimes.discrete = discreteLengthM / (angularRate * rangeM);
    roiTimes.continuous = continuousLengthM / (angularRate * rangeM);
end


function scanTimes = computeEaTimes(phiRad, divergenceRad, rangeM, angularRate)
%COMPUTEEATIMES Calculate the EA scan-time map.

    scanTimes = zeros(numel(phiRad), numel(divergenceRad));
    for phiIndex = 1:numel(phiRad)
        for divergenceIndex = 1:numel(divergenceRad)
            beamDiameterM = 2 * rangeM * tan(divergenceRad(divergenceIndex) / 2);
            k = beamDiameterM / (2 * pi * rangeM);
            lowerPhi = pi / 2 - phiRad(phiIndex) + asin(beamDiameterM / (2 * rangeM));
            upperPhi = pi / 2 - asin(beamDiameterM / (2 * rangeM));
            lowerLimit = lowerPhi / k;
            upperLimit = upperPhi / k;
            integrand = @(x) sqrt(k ^ 2 + sin(k * x) .^ 2);
            integralValue = integral(integrand, lowerLimit, upperLimit);

            scanTimes(phiIndex, divergenceIndex) = integralValue / angularRate + ...
                2 * pi * (sin(upperPhi) + sin(lowerPhi)) / angularRate;
        end
    end
end


function drawScanMap(ax, xMrad, yValues, scanTimes, colorMapName, colorLimits, contourLevels, contourStyles)
%DRAWSCANMAP Draw a formatted scan-time heat map with matching contours.

    imagesc(ax, xMrad, yValues, scanTimes);
    set(ax, 'ColorScale', 'log', 'CLim', colorLimits, 'YDir', 'normal', ...
        'FontName', 'Times New Roman', 'FontSize', 8, 'Box', 'on', ...
        'TickDir', 'out', 'Layer', 'top');
    colormap(ax, colorMapName);
    xlim(ax, [0, 500]);
    xticks(ax, 0:100:500);
    xlabel(ax, '\it{w}\rm_D (mrad)');

    hold(ax, 'on');
    for levelIndex = 1:numel(contourLevels)
        contour(ax, xMrad, yValues, scanTimes, ...
            [contourLevels(levelIndex), contourLevels(levelIndex)], ...
            'w', 'LineStyle', contourStyles{levelIndex}, 'LineWidth', 0.9);
    end
    hold(ax, 'off');
end

% Gaussian beam cross-sections at a propagation distance of 100 m.
wD = 5e-3;                   % Full angular width at 1/e^2 (rad)
wd = 1e-3;                    % Full short/radial width at 1/e^2 (rad)
propagationDistance = 100;    % m

% Use one physical coordinate scale for all three cross-sections.
characteristicRadius = propagationDistance * tan(wD / 2);
viewLimit = 1.5 * characteristicRadius;
coordinate = linspace(-viewLimit, viewLimit, 601);
[X, Y] = meshgrid(coordinate, coordinate);

% Angular coordinates in the observation plane.
radialAngle = atan2(hypot(X, Y), propagationDistance);
longAngle = atan2(Y, propagationDistance);
shortAngle = atan2(X, propagationDistance);

% Profiles used by MC_point_snr.m, MC_line_snr.m and MC_ring_snr.m.
pointHalfWidth = wD / 2;
lineHalfLength = wD / 2;
lineHalfWidth = wd / 2;
ringHalfWidth = wd / 2;
ringRadius = wD / 2 - ringHalfWidth;

pointIntensity = exp(-2 * (radialAngle / pointHalfWidth).^2);
lineIntensity = (abs(longAngle) <= lineHalfLength) .* ...
    exp(-2 * (shortAngle / lineHalfWidth).^2);
ringIntensity = exp(-2 * ((radialAngle - ringRadius) / ...
    ringHalfWidth).^2);

crossSectionFigure = figure("Color", "w", ...
    "Position", [100, 100, 1200, 400]);
crossSectionLayout = tiledlayout(crossSectionFigure, 1, 3, ...
    "Padding", "compact", "TileSpacing", "compact");
axPoint = nexttile(crossSectionLayout);
axLine = nexttile(crossSectionLayout);
axRing = nexttile(crossSectionLayout);

drawIntensity(axPoint, coordinate, pointIntensity);
drawIntensity(axLine, coordinate, lineIntensity);
drawIntensity(axRing, coordinate, ringIntensity);

% White denotes zero intensity; a darkened Fig. 2 green denotes the peak.
numberOfColors = 256;
lowColor = [1, 1, 1];
highColor = [69, 107, 43] / 255*1.5;
blend = linspace(0, 1, numberOfColors).';
intensityMap = lowColor .* (1 - blend) + highColor .* blend;
colormap(crossSectionFigure, intensityMap);

% Closed outlines at the corresponding 1/e^2 locations.
outlineColor = [68, 114, 196] / 255;
spotRadius = propagationDistance * tan(pointHalfWidth);
lineHalfLengthPhysical = propagationDistance * tan(lineHalfLength);
lineHalfWidthPhysical = propagationDistance * tan(lineHalfWidth);
ringInnerRadius = propagationDistance * tan(ringRadius - ringHalfWidth);
ringOuterRadius = propagationDistance * tan(ringRadius + ringHalfWidth);
outlineAngle = linspace(0, 2 * pi, 720);

plot(axPoint, spotRadius * cos(outlineAngle), ...
    spotRadius * sin(outlineAngle), "Color", outlineColor, "LineWidth", 2);

rectangle(axLine, "Position", ...
    [-lineHalfWidthPhysical, -lineHalfLengthPhysical, ...
    2 * lineHalfWidthPhysical, 2 * lineHalfLengthPhysical], ...
    "EdgeColor", outlineColor, "LineWidth", 2);

plot(axRing, ringInnerRadius * cos(outlineAngle), ...
    ringInnerRadius * sin(outlineAngle), "Color", outlineColor, "LineWidth", 2);
plot(axRing, ringOuterRadius * cos(outlineAngle), ...
    ringOuterRadius * sin(outlineAngle), "Color", outlineColor, "LineWidth", 2);

scriptFolder = fileparts(mfilename("fullpath"));
if strlength(scriptFolder) == 0
    scriptFolder = pwd;
end
exportgraphics(crossSectionFigure, fullfile(scriptFolder, "beam.png"), ...
    "Resolution", 300);

% Three-dimensional side views of the beam-divergence geometry.
beam3dFigure = figure("Color", "w", "Position", [120, 120, 1200, 400]);
beam3dLayout = tiledlayout(beam3dFigure, 1, 3, ...
    "Padding", "loose", "TileSpacing", "loose");
beamSurfaceColor = highColor;

axPoint3d = nexttile(beam3dLayout);
drawPointBeam3d(axPoint3d, propagationDistance, pointHalfWidth, ...
    beamSurfaceColor, outlineColor);

axLine3d = nexttile(beam3dLayout);
drawLineBeam3d(axLine3d, propagationDistance, lineHalfLength, ...
    lineHalfWidth, beamSurfaceColor, outlineColor);

axRing3d = nexttile(beam3dLayout);
drawRingBeam3d(axRing3d, propagationDistance, ringRadius, ...
    ringHalfWidth, beamSurfaceColor, outlineColor);

maximumRadius = propagationDistance * tan(ringRadius + ringHalfWidth);
formatBeamAxes(axPoint3d, propagationDistance, maximumRadius);
formatBeamAxes(axLine3d, propagationDistance, maximumRadius);
formatBeamAxes(axRing3d, propagationDistance, maximumRadius);

exportgraphics(beam3dFigure, fullfile(scriptFolder, "beam_3d.png"), ...
    "Resolution", 300);

function drawIntensity(ax, coordinate, intensity)
imagesc(ax, coordinate, coordinate, intensity);
set(ax, "YDir", "normal");
axis(ax, "image");
axis(ax, "off");
clim(ax, [0, 1]);
hold(ax, "on");
end

function drawPointBeam3d(ax, distance, halfAngle, beamColor, outlineColor)
axialCoordinate = linspace(0, distance, 90);
azimuth = linspace(0, 2 * pi, 180);
[axialGrid, azimuthGrid] = meshgrid(axialCoordinate, azimuth);
radiusGrid = axialGrid * tan(halfAngle);

surf(ax, axialGrid, radiusGrid .* cos(azimuthGrid), ...
    radiusGrid .* sin(azimuthGrid), ...
    "FaceColor", beamColor, "FaceAlpha", 0.18, "EdgeColor", "none");
hold(ax, "on");

endRadius = distance * tan(halfAngle);
plot3(ax, distance * ones(size(azimuth)), ...
    endRadius * cos(azimuth), endRadius * sin(azimuth), ...
    "Color", outlineColor, "LineWidth", 2);
drawConeRays(ax, distance, endRadius, outlineColor);
end

function drawLineBeam3d(ax, distance, longHalfAngle, shortHalfAngle, ...
        beamColor, outlineColor)
longHalfSize = distance * tan(longHalfAngle);
shortHalfSize = distance * tan(shortHalfAngle);
vertices = [0, 0, 0; ...
    distance, -longHalfSize, -shortHalfSize; ...
    distance, longHalfSize, -shortHalfSize; ...
    distance, longHalfSize, shortHalfSize; ...
    distance, -longHalfSize, shortHalfSize];
faces = [1, 2, 3; 1, 3, 4; 1, 4, 5; 1, 5, 2; 2, 3, 4; 2, 4, 5];

patch(ax, "Vertices", vertices, "Faces", faces, ...
    "FaceColor", beamColor, "FaceAlpha", 0.18, ...
    "EdgeColor", outlineColor, "LineWidth", 1.4);
hold(ax, "on");

endLoop = [2, 3, 4, 5, 2];
plot3(ax, vertices(endLoop, 1), vertices(endLoop, 2), ...
    vertices(endLoop, 3), "Color", outlineColor, "LineWidth", 2);
end

function drawRingBeam3d(ax, distance, ringAngle, halfWidth, ...
        beamColor, outlineColor)
axialCoordinate = linspace(0, distance, 90);
azimuth = linspace(0, 2 * pi, 180);
[axialGrid, azimuthGrid] = meshgrid(axialCoordinate, azimuth);

innerRadiusGrid = axialGrid * tan(ringAngle - halfWidth);
outerRadiusGrid = axialGrid * tan(ringAngle + halfWidth);
surf(ax, axialGrid, outerRadiusGrid .* cos(azimuthGrid), ...
    outerRadiusGrid .* sin(azimuthGrid), ...
    "FaceColor", beamColor, "FaceAlpha", 0.14, "EdgeColor", "none");
hold(ax, "on");
surf(ax, axialGrid, innerRadiusGrid .* cos(azimuthGrid), ...
    innerRadiusGrid .* sin(azimuthGrid), ...
    "FaceColor", beamColor, "FaceAlpha", 0.20, "EdgeColor", "none");

innerEndRadius = distance * tan(ringAngle - halfWidth);
outerEndRadius = distance * tan(ringAngle + halfWidth);
endRadiusGrid = [innerEndRadius * ones(size(azimuth)); ...
    outerEndRadius * ones(size(azimuth))];
endAzimuthGrid = [azimuth; azimuth];
surf(ax, distance * ones(size(endRadiusGrid)), ...
    endRadiusGrid .* cos(endAzimuthGrid), ...
    endRadiusGrid .* sin(endAzimuthGrid), ...
    "FaceColor", beamColor, "FaceAlpha", 0.35, "EdgeColor", "none");

plot3(ax, distance * ones(size(azimuth)), ...
    innerEndRadius * cos(azimuth), innerEndRadius * sin(azimuth), ...
    "Color", outlineColor, "LineWidth", 2);
plot3(ax, distance * ones(size(azimuth)), ...
    outerEndRadius * cos(azimuth), outerEndRadius * sin(azimuth), ...
    "Color", outlineColor, "LineWidth", 2);
drawConeRays(ax, distance, innerEndRadius, outlineColor);
drawConeRays(ax, distance, outerEndRadius, outlineColor);
end

function drawConeRays(ax, distance, endRadius, color)
for azimuth = [0, pi / 2, pi, 3 * pi / 2]
    plot3(ax, [0, distance], [0, endRadius * cos(azimuth)], ...
        [0, endRadius * sin(azimuth)], ...
        "Color", color, "LineWidth", 1.4);
end
end

function formatBeamAxes(ax, distance, maximumRadius)
xlim(ax, [0, distance]);
ylim(ax, 1.18 * [-maximumRadius, maximumRadius]);
zlim(ax, 1.18 * [-maximumRadius, maximumRadius]);
pbaspect(ax, [2.8, 1, 1]);
view(ax, [-38, 18]);
axis(ax, "off");
set(ax, "Projection", "perspective");
end

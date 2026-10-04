function scan = multiuav_build_scan(cfg, beam, region, trajectory)
%MULTIUAV_BUILD_SCAN Build ROI or ROE scan directions and line frames.

arguments
    cfg (1,1) struct
    beam (1,1) struct
    region (1,1) string
    trajectory (1,1) string
end

region = upper(region);
trajectory = upper(trajectory);
if ~any(region == ["ROI", "ROE"])
    error('multiuav_build_scan:UnknownRegion', ...
        'region must be ROI or ROE.');
end
if region == "ROI" && ~any(trajectory == ["RA", "SP"])
    error('multiuav_build_scan:UnknownTrajectory', ...
        'ROI trajectory must be RA or SP.');
end

if region == "ROI"
    pathType = "raster";
    if trajectory == "SP"
        pathType = "spiral";
    end
    directions = build_ROI_scan_cycle( ...
        cfg.system.roiHalfAngle_deg, beam.wD_rad, ...
        cfg.system.prf_Hz, cfg.system.scanRate_rad_s, pathType);
else
    directions = buildRoeCycle(cfg.system, beam.wD_rad);
    trajectory = "SP";
end

directions = directions ./ vecnorm(directions, 2, 2);
scan = struct();
scan.region = region;
scan.trajectory = trajectory;
scan.B = directions;
scan.cycleLength = size(directions, 1);
scan.cycleTime_s = scan.cycleLength / cfg.system.prf_Hz;

if strcmpi(beam.name, 'line')
    frame = build_ROI_line_scan(directions);
    scan.longAxes = frame.LongAxes;
    scan.shortAxes = frame.ShortAxes;
else
    scan.longAxes = zeros(0, 3);
    scan.shortAxes = zeros(0, 3);
end
end

function cycle = buildRoeCycle(system, wD)
step = system.scanRate_rad_s / system.prf_Hz;
halfWidth = wD / 2;
phiMin = min(max(halfWidth, eps), pi / 2);
phiMax = max(phiMin, pi / 2 - halfWidth);
spiral = sphericalSpiral(phiMin, phiMax, wD, step);
startAzimuth = atan2(spiral(1, 2), spiral(1, 1));
endAzimuth = atan2(spiral(end, 2), spiral(end, 1));
lowerCircle = latitudeCircle(phiMin, step, startAzimuth);
upperCircle = latitudeCircle(phiMax, step, endAzimuth);
reverseSpiral = flipud(spiral);
if size(reverseSpiral, 1) > 2
    reverseInterior = reverseSpiral(2:end-1, :);
else
    reverseInterior = zeros(0, 3);
end
cycle = [lowerCircle; spiral(2:end, :); ...
    upperCircle(2:end, :); reverseInterior];
end

function B = latitudeCircle(phi, step, azimuth0)
circumference = 2 * pi * max(sin(phi), eps);
pointCount = max(4, ceil(circumference / step));
azimuth = azimuth0 + (0:pointCount-1).' .* (2 * pi / pointCount);
B = [sin(phi) .* cos(azimuth), sin(phi) .* sin(azimuth), ...
    cos(phi) .* ones(pointCount, 1)];
end

function B = sphericalSpiral(phiMin, phiMax, pitch, step)
if phiMax <= phiMin + eps
    B = [sin(phiMin), 0, cos(phiMin)];
    return;
end
b = pitch / (2 * pi);
thetaMax = (phiMax - phiMin) / b;
speed = @(theta) sqrt(b.^2 + sin(phiMin + b .* theta).^2);
arcLength = integral(speed, 0, thetaMax, ...
    'RelTol', 1e-8, 'AbsTol', 1e-10);
theta = zeros(max(2, ceil(arcLength / step) + 3), 1);
count = 1;
while theta(count) < thetaMax
    phi = phiMin + b * theta(count);
    dTheta = step / sqrt(b^2 + sin(phi)^2);
    midpoint = min(theta(count) + dTheta / 2, thetaMax);
    midpointPhi = phiMin + b * midpoint;
    dTheta = step / sqrt(b^2 + sin(midpointPhi)^2);
    count = count + 1;
    if count > numel(theta)
        theta = [theta; zeros(100, 1)]; %#ok<AGROW>
    end
    theta(count) = min(theta(count - 1) + dTheta, thetaMax);
end
theta = theta(1:count);
phi = phiMin + b .* theta;
B = [sin(phi) .* cos(theta), sin(phi) .* sin(theta), cos(phi)];
end

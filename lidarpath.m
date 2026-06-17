function beamVectors = lidarpath(f, D, Omega, R)
%LIDARPATH Generate unit pointing vectors on a spherical spiral.
%   beamVectors = LIDARPATH(f, D, Omega, R) returns an N-by-3 array.
%
%   Inputs:
%       f     - Laser pulse repetition frequency, Hz
%       D     - Spacing of the spherical spiral, m
%       Omega - Scan angular speed, rad/s
%       R     - Radius of the scanned hemisphere, m
%
%   The trajectory contains only the spherical spiral, sampled on one pulse
%   time grid with spacing 1/f. The surface scan speed is Omega*R.

arguments
    f (1, 1) double {mustBeFinite, mustBePositive}
    D (1, 1) double {mustBeFinite, mustBePositive}
    Omega (1, 1) double {mustBeFinite, mustBePositive}
    R (1, 1) double {mustBeFinite, mustBePositive}
end

if D > 2 * R
    error("lidarpath:InvalidSpacing", ...
        "D must satisfy D <= 2*R so that the endpoint formula is real.");
end

endpointValue = D * sqrt(R^2 - D^2 / 4) / R^2;
endpointValue = min(max(endpointValue, 0), 1);

% phi is the polar angle measured from the positive z-axis.
phiLower = acos(endpointValue);
phiLower = phiLower+0.5.*asin(endpointValue);
phiUpper = asin(endpointValue)./2;

if phiLower < phiUpper
    error("lidarpath:ReversedEndpoints", ...
        ["The specified D and R make the stated lower endpoint lie " ...
         "above the stated upper endpoint."]);
end

% For the spherical Archimedean spiral, phi = phiLower - k*theta.
% Adjacent turns have nominal surface spacing D.
k = D / (2 * pi * R);

% Arc length from phiLower to phi is
%   s(phi) = R/k * integral_phi^phiLower sqrt(k^2 + sin(u)^2) du.
arcIntegrand = @(phi) sqrt(k^2 + sin(phi).^2);
totalArcLength = R / k * integral(arcIntegrand, phiUpper, phiLower, ...
    "RelTol", 1e-11, "AbsTol", 1e-12);
scanDuration = totalArcLength / (Omega * R);

lastSpiralPulseIndex = floor(scanDuration * f);
spiralTimes = (0:lastSpiralPulseIndex).' / f;
sampleNormalizedArcLengths = Omega * spiralTimes;
totalNormalizedArcLength = totalArcLength / R;

% Invert the arc-length relation using q = s/R and its differential form:
%   dphi/dq = -k / sqrt(k^2 + sin(phi)^2).
% Solving once on an adaptive mesh is much faster than root-finding for
% every laser pulse.
inverseArcLength = @(~, phi) ...
    -k ./ sqrt(k^2 + sin(phi).^2);
solverOptions = odeset("RelTol", 1e-9, "AbsTol", 1e-11);
solution = ode45(inverseArcLength, [0, totalNormalizedArcLength], phiLower, ...
    solverOptions);
spiralPhi = deval(solution, sampleNormalizedArcLengths).';

% Clamp only roundoff-sized endpoint deviations.
spiralPhi = min(phiLower, max(phiUpper, spiralPhi));
spiralTheta = (phiLower - spiralPhi) / k;

beamVectors = [sin(spiralPhi) .* cos(spiralTheta), ...
               sin(spiralPhi) .* sin(spiralTheta), ...
               cos(spiralPhi)];

% Explicit normalization protects against accumulated floating-point error.
beamVectors = beamVectors ./ vecnorm(beamVectors, 2, 2);

% Display the sampling-point distribution on the unit hemisphere.
figure("Name", "LiDAR Spherical Spiral");
[sphereX, sphereY, sphereZ] = sphere(80);
surf(sphereX, sphereY, sphereZ, ...
    "FaceAlpha", 0.08, "EdgeAlpha", 0.12, "FaceColor", [0.3, 0.7, 1]);
hold on
plot3(beamVectors(:, 1), beamVectors(:, 2), beamVectors(:, 3), ...
    "-", "Color", [0.15, 0.35, 0.8], "LineWidth", 0.5);
scatter3(beamVectors(:, 1), beamVectors(:, 2), beamVectors(:, 3), ...
    8, beamVectors(:, 3), "filled");
hold off
axis equal
grid on
xlabel("x");
ylabel("y");
zlabel("z");
title("LiDAR spherical spiral");
view(3);
colorbar
end

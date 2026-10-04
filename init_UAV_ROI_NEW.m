function [P, V] = init_UAV_ROI_NEW(R, jiaodu, v_min, v_max, N)
%INIT_UAV_ROI_NEW Initialize targets on the side surface of a cone ROI.

if nargin < 5
    error('init_UAV_ROI_NEW:NotEnoughInputs', ...
        'Expected R, jiaodu, v_min, v_max, and N.');
end
if R <= 100
    error('init_UAV_ROI_NEW:InvalidRange', ...
        'R must be greater than the 100 m initial range lower bound.');
end
if jiaodu <= 0 || jiaodu >= 90
    error('init_UAV_ROI_NEW:InvalidAngle', ...
        'jiaodu must be between 0 and 90 degrees.');
end

if v_min < 0.2
    v_min = 0.2;
end
v_max = max(v_max, v_min);

alpha = deg2rad(jiaodu);
theta = 2 * pi * rand(N, 1);
minimum_slant_range = 100;
slant_range = sqrt(minimum_slant_range^2 + ...
    (R^2 - minimum_slant_range^2) * rand(N, 1));

P = [ ...
    slant_range .* sin(alpha) .* cos(theta), ...
    slant_range .* sin(alpha) .* sin(theta), ...
    slant_range .* cos(alpha)];

Vmag = v_max .* ones(N, 1);

theta_v = 2 * pi * rand(N, 1);
u_v = 2 * rand(N, 1) - 1;
phi_v = acos(u_v);
Vdir = [ ...
    sin(phi_v) .* cos(theta_v), ...
    sin(phi_v) .* sin(theta_v), ...
    cos(phi_v)];

% The cone side is g(x,y,z)=hypot(x,y)-z*tan(alpha).  Force dg/dt < 0
% at the initial side-surface point, so the first movement enters the cone.
normal = [cos(theta), sin(theta), -tan(alpha) * ones(N, 1)];
normal = normal ./ vecnorm(normal, 2, 2);
outward = sum(Vdir .* normal, 2) >= 0;
Vdir(outward, :) = Vdir(outward, :) - ...
    2 * sum(Vdir(outward, :) .* normal(outward, :), 2) .* normal(outward, :);

V = Vdir .* Vmag;
end

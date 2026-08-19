function [beam_cycle, one_way_path] = build_ROI_scan_cycle( ...
    jiaodu, fasan_D, f, omiga, path_type)
%BUILD_ROI_SCAN_CYCLE Build a constant-angular-step spherical ROI scan.

if jiaodu <= 0 || jiaodu >= 90
    error('build_ROI_scan_cycle:InvalidAngle', ...
        'jiaodu must be between 0 and 90 degrees.');
end
if fasan_D <= 0 || fasan_D >= pi
    error('build_ROI_scan_cycle:InvalidBeamWidth', ...
        'fasan_D must be between 0 and pi radians.');
end
if f <= 0
    error('build_ROI_scan_cycle:InvalidPulseRate', 'f must be positive.');
end
if omiga <= 0
    error('build_ROI_scan_cycle:InvalidAngularSpeed', ...
        'omiga must be positive.');
end

roi_half_angle = deg2rad(jiaodu);
coverage_margin = fasan_D / 2;
path_radius = roi_half_angle + coverage_margin;
if path_radius >= pi
    error('build_ROI_scan_cycle:InvalidPathRadius', ...
        'jiaodu plus fasan_D/2 must be less than pi radians.');
end

angular_step = omiga / f;
switch lower(path_type)
    case 'raster'
        one_way_path = generate_spherical_raster_path( ...
            path_radius, fasan_D, angular_step);
    case 'spiral'
        one_way_path = generate_spherical_spiral_path( ...
            path_radius, fasan_D, angular_step);
    otherwise
        error('build_ROI_scan_cycle:UnknownPathType', ...
            'Unknown path_type: %s. Use raster or spiral.', path_type);
end

beam_cycle = make_ping_pong_cycle(one_way_path);
end

function path = generate_spherical_spiral_path(path_radius, pitch, angular_step)
raw_step = angular_step / 4;
b = pitch / (2 * pi);
phi_max = path_radius / b;
planar_length = 0.5 * b * (phi_max * sqrt(phi_max^2 + 1) + ...
    log(phi_max + sqrt(phi_max^2 + 1)));
point_capacity = max(2, ceil(planar_length / raw_step) + 2);
angular_path = zeros(point_capacity, 2);
point_count = 1;
phi = 0;

while true
    rho = b * phi;
    spherical_speed = sqrt(b^2 + sin(rho)^2);
    dphi = raw_step / max(spherical_speed, eps);
    phi_new = phi + dphi;
    rho_new = b * phi_new;
    if rho_new >= path_radius
        phi_end = path_radius / b;
        point_count = point_count + 1;
        angular_path(point_count, :) = ...
            [path_radius * cos(phi_end), path_radius * sin(phi_end)];
        break;
    end

    point_count = point_count + 1;
    angular_path(point_count, :) = ...
        [rho_new * cos(phi_new), rho_new * sin(phi_new)];
    phi = phi_new;
end

angular_path = angular_path(1:point_count, :);
raw_directions = angular_plane_to_unit(angular_path);
path = resample_spherical_path(raw_directions, angular_step);
end

function path = generate_spherical_raster_path(path_radius, pitch, angular_step)
raw_step = angular_step / 4;
if pitch >= 2 * path_radius
    x_list = 0;
else
    line_count = max(2, ceil(2 * path_radius / pitch));
    x_list = linspace(-path_radius + pitch / 2, ...
        path_radius - pitch / 2, line_count);
end

path_parts = cell(2 * numel(x_list) - 1, 1);
part_count = 0;
previous_end = [];
for i = 1:numel(x_list)
    x0 = x_list(i);
    y_extent = sqrt(max(path_radius^2 - x0^2, 0));
    if mod(i, 2) == 1
        y_start = -y_extent;
        y_end = y_extent;
    else
        y_start = y_extent;
        y_end = -y_extent;
    end

    line_sample_count = max(2, ceil(abs(y_end - y_start) / raw_step) + 1);
    y_line = linspace(y_start, y_end, line_sample_count).';
    line_part = [x0 .* ones(line_sample_count, 1), y_line];

    if ~isempty(previous_end)
        connector = generate_boundary_connector( ...
            previous_end, line_part(1, :), path_radius, raw_step);
        part_count = part_count + 1;
        path_parts{part_count} = connector(2:end, :);
    end

    part_count = part_count + 1;
    if i == 1
        path_parts{part_count} = line_part;
    else
        path_parts{part_count} = line_part(2:end, :);
    end
    previous_end = line_part(end, :);
end

angular_path = vertcat(path_parts{1:part_count});
raw_directions = angular_plane_to_unit(angular_path);
path = resample_spherical_path(raw_directions, angular_step);
end

function connector = generate_boundary_connector(p_start, p_end, radius, raw_step)
phi_start = atan2(p_start(2), p_start(1));
phi_end = atan2(p_end(2), p_end(1));
delta_phi = atan2(sin(phi_end - phi_start), cos(phi_end - phi_start));
arc_length = radius * abs(delta_phi);
sample_count = max(2, ceil(arc_length / raw_step) + 1);
phi = phi_start + linspace(0, delta_phi, sample_count).';
connector = radius .* [cos(phi), sin(phi)];
end

function directions = angular_plane_to_unit(angular_path)
theta = hypot(angular_path(:, 1), angular_path(:, 2));
phi = atan2(angular_path(:, 2), angular_path(:, 1));
sin_theta = sin(theta);
directions = [ ...
    sin_theta .* cos(phi), ...
    sin_theta .* sin(phi), ...
    cos(theta)];
directions = directions ./ vecnorm(directions, 2, 2);
end

function path = resample_spherical_path(raw_path, angular_step)
segment_dot = sum(raw_path(1:end-1, :) .* raw_path(2:end, :), 2);
segment_angle = acos(max(min(segment_dot, 1), -1));
keep = [true; segment_angle > 100 * eps];
raw_path = raw_path(keep, :);

if size(raw_path, 1) == 1
    path = raw_path;
    return;
end

segment_dot = sum(raw_path(1:end-1, :) .* raw_path(2:end, :), 2);
segment_angle = acos(max(min(segment_dot, 1), -1));
cumulative_angle = [0; cumsum(segment_angle)];
total_angle = cumulative_angle(end);
sample_angle = (0:angular_step:total_angle).';
if isempty(sample_angle) || total_angle - sample_angle(end) > 100 * eps(total_angle)
    sample_angle = [sample_angle; total_angle];
end

segment_index = discretize(sample_angle, cumulative_angle);
segment_index(isnan(segment_index)) = numel(segment_angle);
local_fraction = (sample_angle - cumulative_angle(segment_index)) ./ ...
    segment_angle(segment_index);

delta = segment_angle(segment_index);
sin_delta = sin(delta);
weight_start = sin((1 - local_fraction) .* delta) ./ sin_delta;
weight_end = sin(local_fraction .* delta) ./ sin_delta;
path = weight_start .* raw_path(segment_index, :) + ...
    weight_end .* raw_path(segment_index + 1, :);
path = path ./ vecnorm(path, 2, 2);
end

function cycle = make_ping_pong_cycle(path)
if size(path, 1) <= 2
    cycle = path;
else
    cycle = [path; flipud(path(2:end-1, :))];
end
end

function line_scan = build_ROI_line_scan(beam_cycle)
%BUILD_ROI_LINE_SCAN Build spherical long/short axes for a line-beam scan.

if size(beam_cycle, 2) ~= 3 || size(beam_cycle, 1) < 2
    error('build_ROI_line_scan:InvalidCycle', ...
        'beam_cycle must contain at least two N-by-3 direction vectors.');
end

directions = beam_cycle ./ vecnorm(beam_cycle, 2, 2);
cycle_length = size(directions, 1);
previous_index = [cycle_length, 1:cycle_length-1];
next_index = [2:cycle_length, 1];

previous_direction = directions(previous_index, :);
next_direction = directions(next_index, :);
scan_tangent = next_direction - previous_direction;
scan_tangent = scan_tangent - ...
    directions .* sum(directions .* scan_tangent, 2);

tangent_norm = vecnorm(scan_tangent, 2, 2);
degenerate = tangent_norm <= 100 * eps;
if any(degenerate)
    forward_tangent = next_direction - directions;
    forward_tangent = forward_tangent - ...
        directions .* sum(directions .* forward_tangent, 2);
    backward_tangent = directions - previous_direction;
    backward_tangent = backward_tangent - ...
        directions .* sum(directions .* backward_tangent, 2);

    use_forward = degenerate & vecnorm(forward_tangent, 2, 2) > 100 * eps;
    scan_tangent(use_forward, :) = forward_tangent(use_forward, :);
    use_backward = degenerate & ~use_forward;
    scan_tangent(use_backward, :) = backward_tangent(use_backward, :);
end

short_axes = scan_tangent ./ vecnorm(scan_tangent, 2, 2);
long_axes = cross(short_axes, directions, 2);
long_axes = long_axes ./ vecnorm(long_axes, 2, 2);

line_scan = struct( ...
    'Directions', directions, ...
    'LongAxes', long_axes, ...
    'ShortAxes', short_axes);
end

function pulse_count = post_detection_pulse_count( ...
    f, base_angular_span, omiga, detection_range, target_speed)
%POST_DETECTION_PULSE_COUNT Compute the causal post-detection capture length.
%   The target is assumed to move in the same angular direction as the scan.
%   Detection ranges below 100 m are evaluated at 100 m.

TARGET_W = 0.297;
MIN_CALCULATION_RANGE = 100;

calculation_range = max(detection_range, MIN_CALCULATION_RANGE);
target_half_diagonal_angle = atan( ...
    (sqrt(2) * TARGET_W / 2) / calculation_range);
relative_angular_speed = omiga - target_speed / calculation_range;
post_detection_angle = base_angular_span + ...
    2 * target_half_diagonal_angle;

pulse_count = ceil(f * post_detection_angle / relative_angular_speed);
end

function [resolved, ambiguous, rangeBin] = ...
    multiuav_apply_range_resolution(ranges_m, candidate, resolution_m)
%MULTIUAV_APPLY_RANGE_RESOLUTION Apply conservative per-pulse range gating.
% Candidate targets sharing a range bin are ambiguous and cannot be
% credited as individually detected on that pulse.

arguments
    ranges_m (:,:) double
    candidate (:,:) logical
    resolution_m (1,1) double {mustBePositive}
end

if ~isequal(size(ranges_m), size(candidate))
    error('multiuav_apply_range_resolution:SizeMismatch', ...
        'ranges_m and candidate must have the same size.');
end

rangeBin = floor(ranges_m ./ resolution_m);
rangeBin(~isfinite(ranges_m) | ranges_m < 0) = -1;
ambiguous = false(size(candidate));
targetCount = size(candidate, 2);

for firstTarget = 1:targetCount-1
    for secondTarget = firstTarget+1:targetCount
        collision = candidate(:, firstTarget) & ...
            candidate(:, secondTarget) & ...
            rangeBin(:, firstTarget) == rangeBin(:, secondTarget);
        ambiguous(:, firstTarget) = ambiguous(:, firstTarget) | collision;
        ambiguous(:, secondTarget) = ambiguous(:, secondTarget) | collision;
    end
end

resolved = candidate & ~ambiguous;
end

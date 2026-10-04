function [minimumWidth, gridWidth, bracket] = ...
        ea_tornado_min_width(widthMrad, probability, targetProbability)
%EA_TORNADO_MIN_WIDTH First raw probability crossing without extrapolation.
%   Each column of PROBABILITY is treated independently. The function does
%   not smooth the curve or impose a monotone envelope.

arguments
    widthMrad (:,1) double {mustBeFinite,mustBePositive}
    probability double {mustBeNonnegative}
    targetProbability (1,1) double {mustBePositive}
end

assert(issorted(widthMrad, "strictascend") && numel(widthMrad) >= 2, ...
    "EATornado:WidthGrid", ...
    "Width samples must be unique and strictly increasing.");
assert(size(probability, 1) == numel(widthMrad) && ...
    all(isfinite(probability), "all") && all(probability <= 1, "all") && ...
    targetProbability <= 1, ...
    "EATornado:Probability", "Invalid detection-probability array.");

nColumns = size(probability, 2);
minimumWidth = nan(1, nColumns);
gridWidth = nan(1, nColumns);
bracket = nan(nColumns, 2);

for columnIndex = 1:nColumns
    crossingIndex = find( ...
        probability(:, columnIndex) >= targetProbability, 1, "first");
    if isempty(crossingIndex)
        continue;
    end
    gridWidth(columnIndex) = widthMrad(crossingIndex);
    if crossingIndex == 1
        minimumWidth(columnIndex) = widthMrad(1);
        bracket(columnIndex, :) = widthMrad(1);
        continue;
    end

    lowIndex = crossingIndex - 1;
    lowProbability = probability(lowIndex, columnIndex);
    highProbability = probability(crossingIndex, columnIndex);
    bracket(columnIndex, :) = widthMrad([lowIndex, crossingIndex]);
    if highProbability <= lowProbability
        minimumWidth(columnIndex) = widthMrad(crossingIndex);
    else
        fraction = (targetProbability - lowProbability) / ...
            (highProbability - lowProbability);
        minimumWidth(columnIndex) = widthMrad(lowIndex) + fraction * ...
            (widthMrad(crossingIndex) - widthMrad(lowIndex));
    end
end
end

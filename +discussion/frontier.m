function summary = frontier(curve,timeLimit,probability)
%FRONTIER Minimum TESTED feasible width, without monotonicizing responses.
% Simultaneous Clopper-Pearson/Bonferroni bounds cover all widths and both beams
% within this one (timeLimit,probability) comparison. No interpolation.
arguments
    curve (1,1) struct
    timeLimit (1,1) double {mustBeNonnegative}
    probability (1,1) double {mustBeGreaterThan(probability,0),mustBeLessThan(probability,1)}
end
hit = curve.times <= timeLimit;
p = reshape(mean(hit,1),numel(curve.widths_mrad),2);
[lo,hi] = discussion.binomialBounds(p*curve.N,curve.N,0.05/(2*numel(curve.widths_mrad)));
w = curve.widths_mrad;
observed = nan(1,2);
lower = nan(1,2);
upper = nan(1,2);
left = false(1,2);
status = strings(1,2);
for k = 1:2
    observed(k) = firstFeasible(w,p(:,k),probability);
    lower(k) = firstFeasible(w,hi(:,k),probability);
    upper(k) = firstFeasible(w,lo(:,k),probability);
    left(k) = p(1,k) >= probability;
    if isnan(observed(k))
        status(k) = "unattained-on-grid";
    elseif left(k)
        status(k) = "left-boundary";
    else
        status(k) = "interior-grid-minimum";
    end
end
kappa = observed(2)/observed(1);
kappaLow = lower(2)/upper(1);
kappaHigh = upper(2)/lower(1);
% These bounds concern the finite tested grid. They are not uncertainty
% bounds for an unobserved continuous optimum outside the grid.
summary = struct('timeLimit_s',timeLimit,'targetProbability',probability, ...
    'lineWidth_mrad',observed(1),'ringWidth_mrad',observed(2), ...
    'kappa',kappa,'kappaLow',kappaLow,'kappaHigh',kappaHigh, ...
    'lineStatus',status(1),'ringStatus',status(2), ...
    'bothInterior',all(~left & isfinite(observed)), ...
    'certifiedGridAdvantage',isfinite(kappaHigh) && kappaHigh < 1, ...
    'lineCertifiedWidth_mrad',upper(1),'ringCertifiedWidth_mrad',upper(2), ...
    'probability',p,'probabilityLow',lo,'probabilityHigh',hi);
end

function width = firstFeasible(grid,p,target)
index = find(p >= target,1);
width = nan;
if ~isempty(index)
    width = grid(index);
end
end

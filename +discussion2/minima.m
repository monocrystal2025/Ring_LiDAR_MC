function result = minima(curve,required)
%MINIMA First TESTED width meeting an overall-probability requirement.
% Inf firstTime denotes a miss and must never be counted by Inf<=Inf.
arguments
    curve (1,1) struct
    required (1,1) double {mustBeGreaterThan(required,0),mustBeLessThan(required,1)}
end
p = curve.probability;
w = curve.widths_mrad;
[lo,hi] = discussion.binomialBounds(p*curve.N,curve.N,0.05/(2*numel(w)));
indices = nan(1,2);
values = nan(1,2);
lower = nan(1,2);
upper = nan(1,2);
status = strings(1,2);
for k = 1:2
    index = find(p(:,k)>=required,1);
    if isempty(index)
        status(k) = "unattained";
    else
        indices(k) = index;
        values(k) = w(index);
        status(k) = "interior";
        if index == 1
            status(k) = "left-boundary";
        end
    end
    index = find(hi(:,k)>=required,1);
    if ~isempty(index)
        lower(k) = w(index);
    end
    index = find(lo(:,k)>=required,1);
    if ~isempty(index)
        upper(k) = w(index);
    end
end
result = struct('required',required,'indices',indices,'widths_mrad',values, ...
    'status',status,'kappa',values(2)/values(1), ...
    'kappaLow',lower(2)/upper(1),'kappaHigh',upper(2)/lower(1), ...
    'widthLower_mrad',lower,'widthUpper_mrad',upper);
end

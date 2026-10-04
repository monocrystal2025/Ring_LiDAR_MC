function statistics = paired(line,ring,cfg)
%PAIRED Interval for a paired mean difference, retaining trajectory pairing.
% For differences in {-1,0,1}, use simultaneous exact binomial bounds on
% positive/negative discordance. This avoids zero-width bootstrap intervals
% when a small MC ensemble contains no discordant observations.
arguments
    line (:,1) double
    ring (:,1) double
    cfg (1,1) struct
end
assert(numel(line)==numel(ring) && all(isfinite(line)) && all(isfinite(ring)));
difference = ring-line;
n = numel(difference);
if all(ismember(difference,[-1,0,1]))
    [posLo,posHi] = discussion.binomialBounds(nnz(difference==1),n,0.025);
    [negLo,negHi] = discussion.binomialBounds(nnz(difference==-1),n,0.025);
    lo = posLo-negHi;
    hi = posHi-negLo;
    method = "Paired discordance; simultaneous exact binomial bounds (conservative 95%)";
else
    stream = RandStream('mt19937ar','Seed',cfg.seed+53);
    samples = zeros(cfg.bootstrapN,1);
    for j = 1:cfg.bootstrapN
        index = randi(stream,n,n,1);
        samples(j) = mean(difference(index));
    end
    samples = sort(samples);
    lo = samples(max(1,ceil(0.025*cfg.bootstrapN)));
    hi = samples(min(cfg.bootstrapN,ceil(0.975*cfg.bootstrapN)));
    method = "Paired trial bootstrap; pointwise percentile interval";
end
statistics = struct('lineMean',mean(line),'ringMean',mean(ring), ...
    'delta',mean(difference),'low',lo,'high',hi,'N',n,'method',method);
end

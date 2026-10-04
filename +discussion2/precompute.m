function keys = precompute(cfg,indices)
%PRECOMPUTE Warm independent late sensitivity cases in a second worker.
arguments
    cfg (1,1) struct
    indices (1,:) double {mustBeInteger,mustBePositive}
end
cases = discussion2.cases(cfg);
bank = discussion.bank(cfg.N,cfg.seed);
keys = strings(size(indices));
for j = 1:numel(indices)
    item = cases(indices(j));
    fprintf('Precompute case %d: %s x%g\n',indices(j),item.name,item.scale);
    output = discussion2.runCurve(item.config,bank,"dynamic");
    keys(j) = discussion2.hash(output.probability);
end
end

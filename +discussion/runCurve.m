function curve = runCurve(cfg,bank,tag)
%RUNCURVE Paired line/ring response over the predeclared divergence grid.
arguments
    cfg (1,1) struct
    bank (1,1) struct
    tag (1,1) string
end
widths = cfg.widthGrid_mrad;
times = inf(bank.N,numel(widths),2);
singleTimes = times;
uniformTimes = times;
for j = 1:numel(widths)
    for k = 1:2
        names = ["line","ring"];
        b = discussion.beam(names(k),widths(j),cfg.system.beamWidth_rad);
        r = discussion.runCase(cfg,b,bank,tag);
        times(:,j,k) = r.metrics.firstTime;
        singleTimes(:,j,k) = r.metrics.singleTime;
        uniformTimes(:,j,k) = r.metrics.uniformTime;
    end
end
curve = struct('widths_mrad',widths,'times',times, ...
    'singleTimes',singleTimes,'uniformTimes',uniformTimes,'N',bank.N, ...
    'tag',tag,'seed',bank.seed);
end

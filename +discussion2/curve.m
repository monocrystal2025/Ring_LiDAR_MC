function curve = curve(cfg,bank,kind)
%CURVE Width response without smoothing or forcing monotonicity.
arguments
    cfg (1,1) struct
    bank (1,1) struct
    kind (1,1) string {mustBeMember(kind,["dynamic","static"])} = "dynamic"
end
widths = cfg.widthGrid_mrad;
raw = cell(numel(widths),2);
probability = zeros(numel(widths),2);
times = inf(bank.N,numel(widths),2);
pulseMean = nan(numel(widths),2);
clusterMean = pulseMean;
multiProbability = pulseMean;
qualifiedMultiProbability = pulseMean;
for j = 1:numel(widths)
    for k = 1:2
        names = ["line","ring"];
        beam = discussion.beam(names(k),widths(j),cfg.system.beamWidth_rad);
        raw{j,k} = discussion2.evaluate(cfg,beam,bank,kind);
        m = raw{j,k}.metrics;
        detected = isfinite(m.firstTime);
        probability(j,k) = mean(detected);
        times(:,j,k) = m.firstTime;
        if any(detected)
            pulseMean(j,k) = mean(m.effectivePulses(detected));
            clusterMean(j,k) = mean(m.clusters(detected));
            multiProbability(j,k) = mean(m.clusters(detected)>=2);
            qualifiedMultiProbability(j,k) = mean(m.qualifiedClusters(detected)>=2);
        end
    end
    fprintf('Curve %s wD=%g mrad: L=%.4f R=%.4f\n', ...
        kind,widths(j),probability(j,1),probability(j,2));
end
curve = struct('widths_mrad',widths,'raw',{raw},'probability',probability, ...
    'times',times,'pulseMean',pulseMean,'clusterMean',clusterMean, ...
    'multiProbability',multiProbability,'qualifiedMultiProbability', ...
    qualifiedMultiProbability,'N',bank.N,'kind',kind);
end

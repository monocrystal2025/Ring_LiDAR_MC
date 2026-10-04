function [metrics,window] = scoreEvents(events,model,cfg,windowPulses)
%SCOREEVENTS Fixed-window detection and first-detection pulse/cluster counts.
% N_eff and C_eff use S>0, matching advantage_data.m. A second cluster count
% requires each cluster to contribute Q>=clusterQThreshold; it is a temporal
% information diagnostic, not a velocity estimate or an identifiability test.
arguments
    events (:,3) double
    model (1,1) struct
    cfg (1,1) struct
    windowPulses (1,1) double {mustBeInteger,mustBePositive}
end
metrics = struct('firstTime',Inf,'encounterTime',Inf,'effectivePulses',NaN, ...
    'clusters',NaN,'qualifiedClusters',NaN,'clusterSpan_ms',NaN,'maxQ',0);
window = struct('steps',zeros(0,1),'signal',zeros(0,1),'q',zeros(0,1), ...
    'cluster',zeros(0,1),'qualifiedClusterQ',zeros(0,1));
if isempty(events)
    return
end
step = events(:,1);
[signal,noise] = discussion.photons(model,events(:,2),events(:,3),cfg.returnScale);
q = signal.^2./max(signal+noise,realmin);
count = numel(step);
if count == 1
    left = 0;
else
    left = interp1(step,(1:count).',step-windowPulses,'previous',0);
end
cumulative = [0;cumsum(q)];
qWindow = max(0,cumulative(2:end)-cumulative(left+1));
metrics.encounterTime = step(1)/cfg.system.prf_Hz;
metrics.maxQ = max(qWindow);
tail = find(qWindow>=cfg.system.snrThreshold^2,1);
if isempty(tail)
    return
end
first = left(tail)+1;
selected = (first:tail).';
wstep = step(selected);
wq = q(selected);
cluster = cumsum([true;diff(wstep)>1]);
clusterNumber = max(cluster);
clusterQ = zeros(clusterNumber,1);
centroid = zeros(clusterNumber,1);
for j = 1:clusterNumber
    mask = cluster==j;
    clusterQ(j) = sum(wq(mask));
    centroid(j) = sum(wq(mask).*wstep(mask))/max(clusterQ(j),realmin);
end
qualified = clusterQ>=cfg.clusterQThreshold;
metrics.firstTime = step(tail)/cfg.system.prf_Hz;
metrics.effectivePulses = numel(wstep);
metrics.clusters = clusterNumber;
metrics.qualifiedClusters = nnz(qualified);
if nnz(qualified)>=2
    metrics.clusterSpan_ms = (max(centroid(qualified))-min(centroid(qualified))) ...
        /cfg.system.prf_Hz*1000;
end
window = struct('steps',wstep,'signal',signal(selected),'q',wq, ...
    'cluster',cluster,'qualifiedClusterQ',clusterQ(qualified));
end

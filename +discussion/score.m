function [metrics,trace] = score(signal,noise,time,window,threshold,illumination)
%SCORE Causal fixed-window rules and all-trial mechanism diagnostics.
% No selection on successful detection. Misses have firstTime=Inf.
arguments
    signal (:,1) double
    noise (:,1) double
    time (:,1) double
    window (1,1) double {mustBeInteger,mustBePositive}
    threshold (1,1) double {mustBePositive}
    illumination (:,1) double
end
q = signal.^2./max(signal+noise,realmin);
qWindow = causalSum(q,window);
signalWindow = causalSum(signal,window);
noiseWindow = causalSum(noise,window);
uniformQ = signalWindow.^2./max(signalWindow+noiseWindow,realmin);
pulseSquared = causalSum(q.^2,window);
[maxQ,peakIndex] = max(qWindow);
if isempty(peakIndex)
    peakIndex = 1;
    maxQ = 0;
end
% A physical 1/e^2 intensity contour avoids counting tiny Gaussian tails
% as separate opportunities. Geometric prefilter encounters are also saved.
core = illumination >= exp(-2);
starts = core & ~[false;core(1:end-1)];
hitCount = causalSum(double(core),window);
participation = qWindow.^2./max(pulseSquared,realmin);
metrics = struct('firstTime',firstTime(qWindow >= threshold^2,time), ...
    'singleTime',firstTime(q >= threshold^2,time), ...
    'uniformTime',firstTime(uniformQ >= threshold^2,time), ...
    'encounterTime',firstTime(illumination > 0,time), ...
    'coreEncounterTime',firstTime(core,time), ...
    'maxQ',maxQ,'maxSingleQ',max(q),'totalQ',sum(q), ...
    'signalPhotons',sum(signal),'corePulses',nnz(core), ...
    'passages',nnz(starts),'participation',participation(peakIndex), ...
    'peakWindowCorePulses',hitCount(peakIndex), ...
    'exposure_s',time(end)-time(1),'pulses',numel(time));
if nargout > 1
    trace = struct('time',time,'signal',signal,'noise',noise, ...
        'q',q,'qWindow',qWindow,'uniformQ',uniformQ,'core',core);
end
end

function y = causalSum(x,window)
cumulative = [0;cumsum(x)];
ends = (1:numel(x)).';
starts = max(0,ends-window);
y = cumulative(ends+1)-cumulative(starts+1);
% Floating-point subtraction can produce tiny negative sums.
y = max(0,y);
end

function value = firstTime(condition,time)
index = find(condition,1);
value = inf;
if ~isempty(index)
    value = time(index);
end
end

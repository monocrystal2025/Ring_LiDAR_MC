function plotData = parasensitive_predraw(options)
%PARASENSITIVE_PREDRAW Separate line/annular metrics, without ratios.
arguments
    options.Requirements (:,2) double {mustBeFinite,mustBeNonnegative} = [15 .6;30 .7;60 .9]
end
assert(~isempty(options.Requirements) && all(options.Requirements(:,2)<=1));
root = fileparts(mfilename('fullpath'));
s = load(fullfile(root,'parasensitive_data.mat'),'data');
data = s.data;
assert(data.complete,'Complete parasensitive_data first.');
cfg = data.config;
w = cfg.WidthsMrad(:);
nCases = size(data.cases,1);
nReq = size(options.Requirements,1);
curves = repmat(struct('probability',[],'blind',[],'meanClusters',[],'meanPulses',[]),nCases,1);
summary = repmat(struct('minimumWidths',nan(nReq,2), ...
    'widthStatus',strings(nReq,2),'crossingMrad',nan, ...
    'meanBlindPct',nan(1,2),'blindSampleWidths',[], ...
    'meanClusterCount',nan(1,2),'meanPulseCount',nan(1,2),'clusterSampleWidths',[]),nCases,1);
for j=1:nCases
    prob = nan(numel(w),2,nReq);
    blind = nan(numel(w),2);
    clusters = nan(numel(w),2);
    pulses = nan(numel(w),2);
    for d=1:numel(w)
        for b=1:2
            moving = data.dynamic{d,j,b};
            stationary = data.static{d,j,b};
            assert(numel(moving.detect)==cfg.N && numel(stationary.detect)==cfg.N);
            ok = moving.detect~=0 & isfinite(moving.firstTime) & moving.firstTime>=0;
            for r=1:nReq
                t = options.Requirements(r,1);
                prob(d,b,r) = nnz(ok & moving.firstTime<t+eps(t))/cfg.N;
            end
            blind(d,b) = 1-mean(stationary.detect~=0);
            assert(isfield(moving,'effectiveClusterCount'),'Missing effectiveClusterCount.');
            detected = moving.detect~=0;
            if any(detected)
                clusters(d,b) = mean(moving.effectiveClusterCount(detected));
                pulses(d,b) = mean(moving.effectivePulseCount(detected));
            end
        end
    end
    curves(j) = struct('probability',prob,'blind',blind,'meanClusters',clusters,'meanPulses',pulses);
    for r=1:nReq
        for b=1:2
            [summary(j).minimumWidths(r,b),summary(j).widthStatus(r,b)] = ...
                minimumWidth(w,prob(:,b,r),options.Requirements(r,2));
        end
    end
    crossing = lastBlindEquality(w,blind);
    summary(j).crossingMrad = crossing;
    selected = w<crossing; % Common interval, includes all jointly zero-blind points.
    summary(j).blindSampleWidths = w(selected);
    if any(selected)
        summary(j).meanBlindPct = 100*mean(blind(selected,:),1);
    end
    % Equal weight per divergence, conditional mean over detected targets.
    % Use the same finite-width set for both beams; no 0.6 correction.
    selected = all(isfinite(clusters),2) & all(isfinite(pulses),2);
    summary(j).clusterSampleWidths = w(selected);
    if any(selected)
        summary(j).meanClusterCount = mean(clusters(selected,:),1);
        summary(j).meanPulseCount = mean(pulses(selected,:),1);
    end
end
plotData = struct('schemaVersion',5,'N',cfg.N,'preview',cfg.Preview, ...
    'baseline',cfg.baseline,'values',{cfg.values},'displayScale',cfg.displayScale, ...
    'requirements',options.Requirements,'widthsMrad',w,'cases',data.cases, ...
    'caseMap',data.caseMap,'curves',curves,'summary',summary, ...
    'minimumWidthsMrad',nan(5,4,nReq,2),'meanBlindPct',nan(5,4,2), ...
    'meanClusterCount',nan(5,4,2),'meanPulseCount',nan(5,4,2),'crossingMrad',nan(5,4));
plotData.definitions = struct( ...
    'width','Minimum divergence meeting each common time/probability requirement', ...
    'blind','Mean blind fraction over sampled 0<w<largest equality crossing; zero values included', ...
    'clusters','Equal-weight mean across common finite sampled widths of each beam mean cluster count conditional on detection; no correction');
idx = find(strcmp(cfg.sourceNames,'MC_line_snr.m'),1);
plotData.lineShortPrefilterHalfWidthD = nan;
if ~isempty(idx)
    source = cfg.sourceText{idx};
    if ~isempty(regexp(source,'short_prefilter_half_width\s*=\s*GAUSSIAN_CUTOFF_WIDTH\s*\*\s*fasan_d\s*;', 'once'))
        plotData.lineShortPrefilterHalfWidthD = 3;
    elseif ~isempty(regexp(source,'gaussian_width\s*=\s*fasan_d\s*/\s*2\s*;', 'once'))
        plotData.lineShortPrefilterHalfWidthD = 1.5;
    end
end
if plotData.lineShortPrefilterHalfWidthD~=1.5
    warning('parasensitive_predraw:OldGeometry', ...
        'Source data line prefilter is +/- %g d. Regenerate data for +/-1.5d.', ...
        plotData.lineShortPrefilterHalfWidthD);
end
for p=1:4
    for k=1:5
        a = summary(data.caseMap(k,p));
        plotData.minimumWidthsMrad(k,p,:,:) = reshape(a.minimumWidths,1,1,nReq,2);
        plotData.meanBlindPct(k,p,:) = reshape(a.meanBlindPct,1,1,2);
        plotData.meanClusterCount(k,p,:) = reshape(a.meanClusterCount,1,1,2);
        plotData.meanPulseCount(k,p,:) = reshape(a.meanPulseCount,1,1,2);
        plotData.crossingMrad(k,p) = a.crossingMrad;
    end
end
plotData.kappa = plotData.minimumWidthsMrad(:,:,:,2)./plotData.minimumWidthsMrad(:,:,:,1);
plotData.displayColumns = [1 2 4]; % Energy, repetition rate, atmospheric absorption.
plotData.definitions.kappa = 'Minimum annular divergence / minimum line divergence at the same requirement';
plotData.definitions.pulses = 'Equal-weight mean across the same widths as cluster count of detected-sample mean effective pulse count; no 0.6 correction';
save(fullfile(root,'parasensitive_predraw.mat'),'plotData','-v7');
fprintf('Saved separate beam metrics (N=%d).\n',cfg.N);
end
function [width,status] = minimumWidth(w,p,target)
p = cummax(p(:));
index = find(p>=target,1);
width = nan;
status = "unattained";
if isempty(index)
    return
elseif index==1
    width = w(1); % equal_wD.m firstMonotoneCrossing
    status = "grid_lower_bound";
else
    x = p(index-1:index);
    y = w(index-1:index);
    if abs(diff(x))<=eps(max(abs(x)))
        width = y(2);
    else
        width = interp1(x,y,target,'linear');
    end
    status = "ok";
end
end

function crossing = lastBlindEquality(w,blind)
% Largest root, descending search, including jointly zero-blind equality.
crossing = nan;
delta = blind(:,1)-blind(:,2);
for k = numel(w):-1:1
    if delta(k)==0
        crossing = w(k);
        return
    elseif k>1 && delta(k)*delta(k-1)<0
        crossing = interp1(delta(k-1:k),w(k-1:k),0,'linear');
        return
    end
end
end

function report = parasensitive_method_audit()
%PARASENSITIVE_METHOD_AUDIT Compare processed results with paper algorithms.
% Uses the saved simulation, without rerunning Monte Carlo trials.
root = fileparts(mfilename('fullpath'));
source = load(fullfile(root,'parasensitive_data.mat'),'data');
data = source.data;
current = load(fullfile(root,'parasensitive_predraw.mat'),'plotData');
normal = parasensitive_predraw(Requirements=current.plotData.requirements);
cleanup = onCleanup(@() restorePlotData(root,normal));
normalChecks = compareWithPaper(data,normal);
edge = parasensitive_predraw(Requirements=[0 0;0 1;600 0.001]);
edgeChecks = compareWithPaper(data,edge);
assert(all(edge.widthReductionPct(:,:,1)==0,'all'));
assert(all(isnan(edge.widthReductionPct(:,:,2)),'all'));
assert(all(normal.meanBlindReductionPct(4:5,2)==0));
report = struct('normalComparisons',normalChecks,'edgeComparisons',edgeChecks, ...
    'passed',true,'lineDisplayFactor',normal.linePulseDisplayFactor);
save(fullfile(root,'parasensitive_method_audit.mat'),'report');
fprintf('METHOD_AUDIT_PASS: %d normal + %d edge comparisons.\n',normalChecks,edgeChecks);
end

function count = compareWithPaper(data,p)
count = 0;
w = p.widthsMrad;
for j = 1:numel(p.summary)
    expectedWidths = nan(size(p.requirements));
    rawMeans = nan(numel(w),2);
    for k = 1:numel(w)
        for b = 1:2
            raw = data.dynamic{k,j,b};
            valid = raw.detect~=0 & isfinite(raw.firstTime) & raw.firstTime>=0;
            for r = 1:size(p.requirements,1)
                t = p.requirements(r,1);
                bins = histcounts(raw.firstTime(valid),[-inf,t+eps(t),inf]);
                expected = bins(1)/numel(raw.detect);
                assert(abs(p.curves(j).probability(k,b,r)-expected)<1e-12);
                count = count+1;
            end
            stationary = data.static{k,j,b};
            assert(abs(p.curves(j).blind(k,b)-(1-mean(stationary.detect~=0)))<1e-12);
            rawMeans(k,b) = mean(raw.effectivePulseCount(raw.detect~=0));
            count = count+1;
        end
    end
    for r = 1:size(p.requirements,1)
        for b = 1:2
            expectedWidths(r,b) = firstMonotoneCrossing(w, ...
                p.curves(j).probability(:,b,r),p.requirements(r,2));
        end
    end
    assert(isequaln(expectedWidths,p.summary(j).minimumWidths));
    expectedKappa = expectedWidths(:,2)./expectedWidths(:,1);
    assert(isequaln(expectedKappa.',p.summary(j).kappa));
    assert(isequaln(100*(1-expectedKappa.'),p.summary(j).widthReductionPct));
    assert(isequaln(rawMeans,p.curves(j).meanPulses));
    displayMeans = rawMeans;
    displayMeans(:,1) = 0.6*displayMeans(:,1); % advantage_draw ROE line correction
    assert(isequaln(displayMeans,p.curves(j).displayMeanPulses));
    ratios = displayMeans(:,2)./displayMeans(:,1);
    ratios = ratios(isfinite(ratios));
    assert(abs(max(ratios)-p.summary(j).maxPulseRatio)<1e-12);
    assert(abs(p.summary(j).maxRawPulseRatio/0.6-p.summary(j).maxPulseRatio)<1e-12);
    count = count+7;
end
end

function restorePlotData(root,plotData)
save(fullfile(root,'parasensitive_predraw.mat'),'plotData','-v7');
end
function minimumWD = firstMonotoneCrossing(wD, probability, target)
wD = wD(:);
probability = cummax(probability(:));
crossingIndex = find(probability >= target, 1, "first");

if isempty(crossingIndex)
    minimumWD = nan;
elseif crossingIndex == 1
    minimumWD = wD(1);
else
    x = probability(crossingIndex-1:crossingIndex);
    y = wD(crossingIndex-1:crossingIndex);
    if abs(diff(x)) <= eps(max(abs(x)))
        minimumWD = y(2);
    else
        minimumWD = interp1(x, y, target, "linear");
    end
end
end



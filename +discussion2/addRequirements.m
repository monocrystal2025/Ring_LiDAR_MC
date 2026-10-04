function result = addRequirements(result)
%ADDREQUIREMENTS Report a range of common requirements, including null and failed cases.
arguments
    result (1,1) struct
end
levels = [0.70,0.75,0.80,0.85,0.90,0.95];
requirements = cell(numel(result.curves),numel(levels));
rows = cell(numel(requirements),1);
for a = 1:numel(result.curves)
    for b = 1:numel(levels)
        f = discussion2.minima(result.curves{a},levels(b));
        requirements{a,b} = f;
        rows{b+(a-1)*numel(levels)} = struct('Case',a,'Name',result.cases(a).name, ...
            'Scale',result.cases(a).scale,'RequiredPd',levels(b), ...
            'LineWidth_mrad',f.widths_mrad(1),'RingWidth_mrad',f.widths_mrad(2), ...
            'Kappa',f.kappa,'KappaLow95',f.kappaLow,'KappaHigh95',f.kappaHigh, ...
            'LineStatus',f.status(1),'RingStatus',f.status(2));
    end
end
result.requirementLevels = levels;
result.requirements = requirements;
result.requirementSummary = struct2table(vertcat(rows{:}));
result.paperPairWidths_mrad = [195,120]; % Original manuscript representative pair.
result.paperPair = summarizePaperPair(result);
save(fullfile(result.config.outputRoot,"tmp_1.mat"),'result','-v7.3');
writetable(result.requirementSummary,fullfile(result.config.outputRoot,"tmp_1_requirements.csv"));
end

function pair = summarizePaperPair(result)
curve = result.curves{1};
indexLine = find(curve.widths_mrad==195,1);
indexRing = find(curve.widths_mrad==120,1);
if isempty(indexLine) || isempty(indexRing)
    pair = struct('available',false);
    return
end
line = curve.raw{indexLine,1}.metrics;
ring = curve.raw{indexRing,2}.metrics;
hitLine = isfinite(line.firstTime);
hitRing = isfinite(ring.firstTime);
common = hitLine & hitRing;
cfg = result.config;
clusterDifference = discussion.paired(line.clusters(common),ring.clusters(common),cfg);
pulseDifference = discussion.paired(line.effectivePulses(common),ring.effectivePulses(common),cfg);
pDifference = discussion.paired(double(hitLine),double(hitRing),cfg);
ringCounts = ring.clusters(hitRing);
stream = RandStream('mt19937ar','Seed',cfg.seed+781);
sampleMeans = zeros(cfg.bootstrapN,1);
for j = 1:cfg.bootstrapN
    index = randi(stream,numel(ringCounts),numel(ringCounts),1);
    sampleMeans(j) = mean(ringCounts(index));
end
sampleMeans = sort(sampleMeans);
ci = sampleMeans([ceil(0.025*cfg.bootstrapN),ceil(0.975*cfg.bootstrapN)]);
pair = struct('available',true,'widths_mrad',[195,120], ...
    'Pd',pDifference,'commonDetectedN',nnz(common), ...
    'meanPulses',[mean(line.effectivePulses(hitLine)),mean(ring.effectivePulses(hitRing))], ...
    'meanClusters',[mean(line.clusters(hitLine)),mean(ringCounts)], ...
    'ringMeanClusterCI95',ci,'commonClusterDifference',clusterDifference, ...
    'commonPulseDifference',pulseDifference, ...
    'multiClusterConditional',[mean(line.clusters(hitLine)>=2),mean(ringCounts>=2)], ...
    'qualifiedMultiConditional',[mean(line.qualifiedClusters(hitLine)>=2), ...
    mean(ring.qualifiedClusters(hitRing)>=2)], ...
    'qualifiedMultiJoint',[mean(hitLine & line.qualifiedClusters>=2), ...
    mean(hitRing & ring.qualifiedClusters>=2)], ...
    'ringMedianClusterSpan_ms',median(ring.clusterSpan_ms,'omitnan'));
end

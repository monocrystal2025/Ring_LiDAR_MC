function result = computeStudy(cfg)
%COMPUTESTUDY Three manuscript claims across physical parameter changes.
% Overall Pd: dynamic target until its actual upper-half-ball exit.
% Blind fraction: uniform-volume static probes, exactly one complete scan.
% Pulse/cluster counts: first successful fixed window, conditional on Pd.
arguments
    cfg (1,1) struct
end
if ~isfolder(cfg.outputRoot)
    mkdir(cfg.outputRoot);
end
cases = discussion2.cases(cfg);
dynamicBank = discussion.bank(cfg.N,cfg.seed);
staticBank = discussion.bank(cfg.staticN,cfg.seed+1701);
curves = cell(numel(cases),1);
minima = cell(numel(cases),1);
selectedStatic = cell(numel(cases),2);
summaryRows = cell(numel(cases),1);
for j = 1:numel(cases)
    fprintf('tmp_1 case %d/%d: %s x%g\n',j,numel(cases),cases(j).name,cases(j).scale);
    caseCfg = cases(j).config;
    curves{j} = discussion2.runCurve(caseCfg,dynamicBank,"dynamic");
    minima{j} = discussion2.minima(curves{j},cfg.targetProbability);
    iso = minima{j};
    row = struct('Case',j,'Factor',cases(j).factor,'Name',cases(j).name, ...
        'Level',cases(j).level,'Scale',cases(j).scale, ...
        'LineWidth_mrad',iso.widths_mrad(1),'RingWidth_mrad',iso.widths_mrad(2), ...
        'Kappa',iso.kappa,'KappaLow95',iso.kappaLow,'KappaHigh95',iso.kappaHigh, ...
        'LineStatus',iso.status(1),'RingStatus',iso.status(2), ...
        'LineP',NaN,'RingP',NaN,'LineBlind',NaN,'RingBlind',NaN, ...
        'LinePulses',NaN,'RingPulses',NaN,'LineClusters',NaN,'RingClusters',NaN, ...
        'LineMulti',NaN,'RingMulti',NaN,'LineQualifiedMulti',NaN,'RingQualifiedMulti',NaN, ...
        'LineDetectedN',0,'RingDetectedN',0,'DynamicN',cfg.N,'StaticN',cfg.staticN);
    for k = 1:2
        if ~isfinite(iso.indices(k))
            continue
        end
        names = ["line","ring"];
        labels = ["Line","Ring"];
        index = iso.indices(k);
        beam = discussion.beam(names(k),iso.widths_mrad(k),cfg.system.beamWidth_rad);
        selectedStatic{j,k} = discussion2.evaluate(caseCfg,beam,staticBank,"static");
        row.(labels(k)+"P") = curves{j}.probability(index,k);
        row.(labels(k)+"Blind") = mean(~isfinite(selectedStatic{j,k}.metrics.firstTime));
        row.(labels(k)+"Pulses") = curves{j}.pulseMean(index,k);
        row.(labels(k)+"Clusters") = curves{j}.clusterMean(index,k);
        row.(labels(k)+"Multi") = curves{j}.multiProbability(index,k);
        row.(labels(k)+"QualifiedMulti") = curves{j}.qualifiedMultiProbability(index,k);
        row.(labels(k)+"DetectedN") = nnz(isfinite(curves{j}.raw{index,k}.metrics.firstTime));
    end
    summaryRows{j} = row;
    progress = struct('caseIndex',j,'caseCount',numel(cases),'name',cases(j).name, ...
        'scale',cases(j).scale,'time',string(datetime('now')),'state',"running");
    save(fullfile(cfg.outputRoot,"tmp_1_progress.mat"),'progress');
end
% A priori examples: half, nominal, and double pulse energy. These do not
% select cases by favorable outcomes. All three use identical geometry.
exampleIndices = [2,1,3];
staticCurves = cell(1,3);
for j = 1:3
    staticCurves{j} = discussion2.runCurve(cases(exampleIndices(j)).config,staticBank,"static");
end
summary = struct2table(vertcat(summaryRows{:}));
result = struct('config',cfg,'cases',cases,'curves',{curves}, ...
    'minima',{minima},'selectedStatic',{selectedStatic}, ...
    'exampleIndices',exampleIndices,'staticCurves',{staticCurves},'summary',summary, ...
    'definitions',struct('Pd',"Overall first detection before actual ROE exit", ...
    'blind',"Undetected volume fraction after one scan; static uniform-volume probes", ...
    'pulses',"S>0 pulse count in the first successful fixed window, conditional on detection", ...
    'clusters',"Contiguous S>0 runs in that window; separate if at least one zero-signal pulse", ...
    'qualifiedClusters',"Clusters whose sum S^2/(S+B) is at least 1"));
save(fullfile(cfg.outputRoot,"tmp_1.mat"),'result','-v7.3');
writetable(summary,fullfile(cfg.outputRoot,"tmp_1_summary.csv"));
progress.state = "complete";
save(fullfile(cfg.outputRoot,"tmp_1_progress.mat"),'progress');
end

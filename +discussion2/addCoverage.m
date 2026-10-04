function result = addCoverage(result)
%ADDCOVERAGE Attach both blind-zone definitions and their selected-width values.
arguments
    result (1,1) struct
end
file = fullfile(result.config.outputRoot,"tmp_1_geometric_coverage.mat");
if isfile(file)
    loaded = load(file,'result');
    assert(isequaln(loaded.result.config,result.config),'discussion2:CoverageConfig', ...
        'Coverage and sensitivity configurations do not match.');
    geometry = loaded.result;
else
    geometry = discussion2.computeCoverage(result.config);
end
result.geometricCoverage = geometry;
result.summary.LineGeometricBlind = nan(height(result.summary),1);
result.summary.RingGeometricBlind = nan(height(result.summary),1);
for j = 1:numel(result.cases)
    index = find(geometry.caseIndices==j,1);
    if isempty(index)
        index = 1; % Energy, range, speed, target size and receiver do not change point-centre geometry.
    end
    curve = geometry.curves{index};
    iso = result.minima{j};
    for k = 1:2
        if isfinite(iso.indices(k))
            names = ["LineGeometricBlind","RingGeometricBlind"];
            result.summary.(names(k))(j) = curve.blind(iso.indices(k),k);
        end
    end
end
result.definitions.geometricBlind = ...
    "One-cycle 1/e^2 point-centre footprint coverage, as static_MC_EA_NEW.m; no SNR or finite-target criterion.";
result.definitions.snrBlind = result.definitions.blind;
result.definitions.blind = result.definitions.geometricBlind;
save(fullfile(result.config.outputRoot,"tmp_1.mat"),'result','-v7.3');
writetable(result.summary,fullfile(result.config.outputRoot,"tmp_1_summary.csv"));
end

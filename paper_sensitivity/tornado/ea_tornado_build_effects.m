function effects = ea_tornado_build_effects(levels)
%EA_TORNADO_BUILD_EFFECTS Build endpoint tornado ranking from level results.

arguments
    levels table
end

required = ["Factor", "FactorLabel", "Scale", "Value", "Unit", ...
    "Kappa", "KappaCI_Low", "KappaCI_High", ...
    "DeltaKappa", "DeltaKappaCI_Low", "DeltaKappaCI_High", ...
    "BothFeasibleBootstrap"];
assert(all(ismember(required, string(levels.Properties.VariableNames))), ...
    "EATornado:LevelSchema", "Level table is missing required columns.");

factorNames = unique(levels.Factor, "stable");
records = repmat(emptyRecord(), numel(factorNames), 1);
for factorIndex = 1:numel(factorNames)
    rows = levels(levels.Factor == factorNames(factorIndex), :);
    low = rows(abs(rows.Scale - 0.8) < 1e-12, :);
    midLow = rows(abs(rows.Scale - 0.9) < 1e-12, :);
    baseline = rows(abs(rows.Scale - 1.0) < 1e-12, :);
    midHigh = rows(abs(rows.Scale - 1.1) < 1e-12, :);
    high = rows(abs(rows.Scale - 1.2) < 1e-12, :);
    assert(all([height(low), height(midLow), height(baseline), ...
        height(midHigh), height(high)] == 1), ...
        "EATornado:LevelDesign", ...
        "Every factor must contain the five configured levels.");

    endpointValues = [low.Kappa, baseline.Kappa, high.Kappa];
    finiteEndpoints = endpointValues(isfinite(endpointValues));
    if isempty(finiteEndpoints)
        rangeLow = nan;
        rangeHigh = nan;
        span = nan;
    else
        rangeLow = min(finiteEndpoints);
        rangeHigh = max(finiteEndpoints);
        span = rangeHigh - rangeLow;
    end
    intermediateValues = [midLow.Kappa, midHigh.Kappa];
    finiteIntermediate = intermediateValues(isfinite(intermediateValues));
    nonlinearOutside = false;
    if isfinite(rangeLow) && ~isempty(finiteIntermediate)
        tolerance = 1e-12 * max(1, max(abs(finiteEndpoints)));
        nonlinearOutside = any(finiteIntermediate < rangeLow - tolerance | ...
            finiteIntermediate > rangeHigh + tolerance);
    end

    record = emptyRecord();
    record.Factor = low.Factor;
    record.FactorLabel = low.FactorLabel;
    record.LowValue = low.Value;
    record.HighValue = high.Value;
    record.Unit = low.Unit;
    record.BaselineKappa = baseline.Kappa;
    record.LowKappa = low.Kappa;
    record.LowKappaCI_Low = low.KappaCI_Low;
    record.LowKappaCI_High = low.KappaCI_High;
    record.HighKappa = high.Kappa;
    record.HighKappaCI_Low = high.KappaCI_Low;
    record.HighKappaCI_High = high.KappaCI_High;
    record.LowDeltaKappa = low.DeltaKappa;
    record.LowDeltaCI_Low = low.DeltaKappaCI_Low;
    record.LowDeltaCI_High = low.DeltaKappaCI_High;
    record.HighDeltaKappa = high.DeltaKappa;
    record.HighDeltaCI_Low = high.DeltaKappaCI_Low;
    record.HighDeltaCI_High = high.DeltaKappaCI_High;
    record.LowBothFeasible = low.BothFeasibleBootstrap;
    record.HighBothFeasible = high.BothFeasibleBootstrap;
    record.TornadoLow = rangeLow;
    record.TornadoHigh = rangeHigh;
    record.TornadoSpan = span;
    record.IntermediateOutsideEndpointRange = nonlinearOutside;
    record.RobustAdvantageAtEndpoints = all(isfinite([low.Kappa, high.Kappa])) && ...
        low.Kappa < 1 && high.Kappa < 1;
    records(factorIndex) = record;
end

effects = struct2table(records);
sortSpan = effects.TornadoSpan;
sortSpan(~isfinite(sortSpan)) = -inf;
[~, order] = sort(sortSpan, "descend");
effects = effects(order, :);
effects.Rank = (1:height(effects)).';
effects = movevars(effects, "Rank", "Before", 1);
end

function record = emptyRecord()
record = struct( ...
    "Factor", "", "FactorLabel", "", ...
    "LowValue", nan, "HighValue", nan, "Unit", "", ...
    "BaselineKappa", nan, ...
    "LowKappa", nan, "LowKappaCI_Low", nan, ...
    "LowKappaCI_High", nan, ...
    "HighKappa", nan, "HighKappaCI_Low", nan, ...
    "HighKappaCI_High", nan, ...
    "LowDeltaKappa", nan, "LowDeltaCI_Low", nan, ...
    "LowDeltaCI_High", nan, ...
    "HighDeltaKappa", nan, "HighDeltaCI_Low", nan, ...
    "HighDeltaCI_High", nan, ...
    "LowBothFeasible", nan, "HighBothFeasible", nan, ...
    "TornadoLow", nan, "TornadoHigh", nan, "TornadoSpan", nan, ...
    "IntermediateOutsideEndpointRange", false, ...
    "RobustAdvantageAtEndpoints", false);
end

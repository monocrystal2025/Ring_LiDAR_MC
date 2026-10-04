function audit = ea_tornado_verify(outputDirectory)
%EA_TORNADO_VERIFY Verify completeness and cross-file consistency.

arguments
    outputDirectory (1,1) string
end

designSaved = load(fullfile(outputDirectory, "design.mat"), "design");
analysisSaved = load(fullfile(outputDirectory, ...
    "ea_tornado_analysis.mat"), "analysis");
design = designSaved.design;
analysis = analysisSaved.analysis;

expectedFiles = [ ...
    "ea_tornado_level_results.csv", ...
    "ea_tornado_effects.csv", ...
    "ea_tornado_probability_curves.csv", ...
    "physical_direction_checks.csv", ...
    "ea_tornado_analysis.mat", ...
    "methodology_and_results.md", ...
    fullfile("figures", "ea_tornado_sensitivity.png"), ...
    fullfile("figures", "ea_tornado_sensitivity.pdf"), ...
    fullfile("figures", "ea_tornado_sensitivity.fig")];
filesPresent = arrayfun(@(name) isfile(fullfile(outputDirectory, name)), ...
    expectedFiles);

levelCsv = readtable(fullfile(outputDirectory, ...
    "ea_tornado_level_results.csv"), TextType="string");
effectCsv = readtable(fullfile(outputDirectory, ...
    "ea_tornado_effects.csv"), TextType="string");
curveCsv = readtable(fullfile(outputDirectory, ...
    "ea_tornado_probability_curves.csv"), TextType="string");

audit = struct();
audit.SchemaVersion = 1;
audit.CreatedAt = string(datetime("now", ...
    "Format", "yyyy-MM-dd HH:mm:ss Z"));
audit.RequiredCheckpoints = nnz(design.required);
audit.CompletedCheckpoints = nnz(design.complete & design.required);
audit.AllRequiredComplete = all(design.complete(design.required));
audit.CaseCount = numel(design.cases);
audit.LevelRowCount = height(levelCsv);
audit.EffectRowCount = height(effectCsv);
audit.CurveRowCount = height(curveCsv);
audit.FilesPresent = all(filesPresent);
audit.MissingFiles = cellstr(expectedFiles(~filesPresent));
audit.LevelsMatchMat = tablesEquivalent(levelCsv, analysis.levels);
audit.EffectsMatchMat = tablesEquivalent(effectCsv, analysis.effects);
audit.ExpectedCaseCount = audit.CaseCount == 17;
audit.ExpectedLevelRows = audit.LevelRowCount == 20;
audit.ExpectedEffectRows = audit.EffectRowCount == 4;
audit.FiniteEndpointFactors = nnz(isfinite(analysis.effects.LowKappa) & ...
    isfinite(analysis.effects.HighKappa));
audit.LowFeasibilityFactors = nnz( ...
    analysis.effects.LowBothFeasible < ...
    design.config.feasibilityThreshold | ...
    analysis.effects.HighBothFeasible < ...
    design.config.feasibilityThreshold);
audit.UnattainedLevelRows = nnz(~isfinite(analysis.levels.Kappa));
audit.NonlinearAuditFlags = nnz( ...
    analysis.effects.IntermediateOutsideEndpointRange);
audit.PhysicalDirectionViolations = sum( ...
    analysis.physicalChecks.Violations);
audit.Passed = audit.AllRequiredComplete && audit.FilesPresent && ...
    audit.ExpectedCaseCount && audit.ExpectedLevelRows && ...
    audit.ExpectedEffectRows && audit.LevelsMatchMat && ...
    audit.EffectsMatchMat;

jsonFile = fullfile(outputDirectory, "verification.json");
fileId = fopen(jsonFile, "w", "n", "UTF-8");
assert(fileId >= 0, "EATornado:VerificationFile", ...
    "Unable to create verification report.");
cleanup = onCleanup(@() fclose(fileId)); %#ok<NASGU>
fprintf(fileId, "%s", jsonencode(audit, PrettyPrint=true));
assert(audit.Passed, "EATornado:VerificationFailed", ...
    "EA tornado output verification failed.");
end

function tf = tablesEquivalent(first, second)
firstNames = string(first.Properties.VariableNames);
secondNames = string(second.Properties.VariableNames);
if height(first) ~= height(second) || ~isequal(firstNames, secondNames)
    tf = false;
    return;
end
tf = true;
for variableIndex = 1:numel(firstNames)
    firstValue = first.(firstNames(variableIndex));
    secondValue = second.(secondNames(variableIndex));
    if (isnumeric(firstValue) || islogical(firstValue)) && ...
            (isnumeric(secondValue) || islogical(secondValue))
        a = double(firstValue);
        b = double(secondValue);
        sameNaN = isnan(a) & isnan(b);
        finite = isfinite(a) & isfinite(b);
        tolerance = 1e-10 .* max(1, max(abs(a), abs(b)));
        equal = sameNaN | (finite & abs(a - b) <= tolerance) | ...
            (isinf(a) & a == b);
        tf = tf && all(equal, "all");
    else
        tf = tf && isequal(string(firstValue), string(secondValue));
    end
    if ~tf
        return;
    end
end
end

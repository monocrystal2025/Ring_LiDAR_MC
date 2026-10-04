clc;
clearvars;

% Read the raw EA and ROI Monte Carlo results and cache the statistics.
% Every saved variable is an N-by-1 double vector. Its kth element is the
% mean over all successful detections at w_D = wDMrad(k).

scriptDir = fileparts(mfilename("fullpath"));
if strlength(scriptDir) == 0
    scriptDir = pwd;
end

eaResultDir = "D:\lzx\MC_RESULTS";
roiResultDir = fullfile(scriptDir, "MC_ROI_RESULTS_SNR_NEW");
outputFile = fullfile(scriptDir, "advantage.mat");
if ~isfolder(eaResultDir)
    error("EA result folder not found: %s", eaResultDir);
end
if ~isfolder(roiResultDir)
    error("ROI result folder not found: %s", roiResultDir);
end

wDMrad = (5:5:400).';
innerDiameterMrad = 1;
resultColumnIndex = 1;
linePulseCountFactor = 0.6; % Apply to line pulse means saved in advantage.mat.
beamTypes = ["POINT", "LINE", "RING"];
trackTypes = ["RA", "SP"];
numberOfDiameters = numel(wDMrad);
numberOfBeams = numel(beamTypes);
numberOfTracks = numel(trackTypes);

eaMeanPulseCount = nan(numberOfDiameters, numberOfBeams);
eaMeanDetectionRange = nan(numberOfDiameters, numberOfBeams);
roiMeanDetectionRange = nan(numberOfDiameters, numberOfBeams, numberOfTracks);
eaMeanPulseClusterCount = nan(numberOfDiameters, numberOfBeams);
roiMeanPulseCount = nan(numberOfDiameters, numberOfBeams, numberOfTracks);
roiMeanPulseClusterCount = ...
    nan(numberOfDiameters, numberOfBeams, numberOfTracks);

for beamIndex = 1:numberOfBeams
    beamType = beamTypes(beamIndex);
    beamCode = beamCodeForType(beamType);
    detectVariable = "detect_" + beamCode;
    pulseVariable = "effective_pulses_" + beamCode;
    fprintf("Processing %s EA results...\n", beamType);

    for diameterIndex = 1:numberOfDiameters
        resultFile = eaResultFileName(eaResultDir, beamType, ...
            wDMrad(diameterIndex), innerDiameterMrad);
        [eaMeanPulseCount(diameterIndex, beamIndex), ...
            eaMeanPulseClusterCount(diameterIndex, beamIndex), ...
            eaMeanDetectionRange(diameterIndex, beamIndex)] = ...
            processResultFile(resultFile, detectVariable, pulseVariable, ...
            resultColumnIndex);
    end
end

for trackIndex = 1:numberOfTracks
    trackType = trackTypes(trackIndex);
    for beamIndex = 1:numberOfBeams
        beamType = beamTypes(beamIndex);
        beamCode = beamCodeForType(beamType);
        variableSuffix = beamCode + "_" + trackType;
        detectVariable = "detect_" + variableSuffix;
        pulseVariable = "effective_pulses_" + variableSuffix;
        fprintf("Processing %s ROI-%s results...\n", beamType, trackType);

        for diameterIndex = 1:numberOfDiameters
            resultFile = roiResultFileName(roiResultDir, beamType, ...
                trackType, wDMrad(diameterIndex), innerDiameterMrad);
            [roiMeanPulseCount(diameterIndex, beamIndex, trackIndex), ...
                roiMeanPulseClusterCount( ...
                diameterIndex, beamIndex, trackIndex), ...
                roiMeanDetectionRange(diameterIndex, beamIndex, trackIndex)] = ...
                processResultFile(resultFile, detectVariable, ...
                pulseVariable, resultColumnIndex);
        end
    end
end

spotEffectivePulseCountMean = double(eaMeanPulseCount(:, 1));
spotEffectivePulseClusterCountMean = ...
    double(eaMeanPulseClusterCount(:, 1));
lineEffectivePulseCountMean = linePulseCountFactor .* double(eaMeanPulseCount(:, 2));
lineEffectivePulseClusterCountMean = ...
    double(eaMeanPulseClusterCount(:, 2));
annularEffectivePulseCountMean = double(eaMeanPulseCount(:, 3));
annularEffectivePulseClusterCountMean = ...
    double(eaMeanPulseClusterCount(:, 3));

roiRaSpotEffectivePulseCountMean = double(roiMeanPulseCount(:, 1, 1));
roiRaSpotEffectivePulseClusterCountMean = ...
    double(roiMeanPulseClusterCount(:, 1, 1));
roiRaLineEffectivePulseCountMean = ...
    linePulseCountFactor .* double(roiMeanPulseCount(:, 2, 1));
roiRaLineEffectivePulseClusterCountMean = ...
    double(roiMeanPulseClusterCount(:, 2, 1));
roiRaAnnularEffectivePulseCountMean = double(roiMeanPulseCount(:, 3, 1));
roiRaAnnularEffectivePulseClusterCountMean = ...
    double(roiMeanPulseClusterCount(:, 3, 1));

roiSpSpotEffectivePulseCountMean = double(roiMeanPulseCount(:, 1, 2));
roiSpSpotEffectivePulseClusterCountMean = ...
    double(roiMeanPulseClusterCount(:, 1, 2));
roiSpLineEffectivePulseCountMean = ...
    linePulseCountFactor .* double(roiMeanPulseCount(:, 2, 2));
roiSpLineEffectivePulseClusterCountMean = ...
    double(roiMeanPulseClusterCount(:, 2, 2));
roiSpAnnularEffectivePulseCountMean = double(roiMeanPulseCount(:, 3, 2));
roiSpAnnularEffectivePulseClusterCountMean = ...
    double(roiMeanPulseClusterCount(:, 3, 2));

spotDetectionRangeMean = eaMeanDetectionRange(:, 1);
lineDetectionRangeMean = eaMeanDetectionRange(:, 2);
annularDetectionRangeMean = eaMeanDetectionRange(:, 3);
roiRaSpotDetectionRangeMean = roiMeanDetectionRange(:, 1, 1);
roiRaLineDetectionRangeMean = roiMeanDetectionRange(:, 2, 1);
roiRaAnnularDetectionRangeMean = roiMeanDetectionRange(:, 3, 1);
roiSpSpotDetectionRangeMean = roiMeanDetectionRange(:, 1, 2);
roiSpLineDetectionRangeMean = roiMeanDetectionRange(:, 2, 2);
roiSpAnnularDetectionRangeMean = roiMeanDetectionRange(:, 3, 2);

save(outputFile, ...
    "spotDetectionRangeMean", "lineDetectionRangeMean", ...
    "annularDetectionRangeMean", "roiRaSpotDetectionRangeMean", ...
    "roiRaLineDetectionRangeMean", "roiRaAnnularDetectionRangeMean", ...
    "roiSpSpotDetectionRangeMean", "roiSpLineDetectionRangeMean", ...
    "roiSpAnnularDetectionRangeMean", ...
    "spotEffectivePulseCountMean", ...
    "spotEffectivePulseClusterCountMean", ...
    "lineEffectivePulseCountMean", ...
    "lineEffectivePulseClusterCountMean", ...
    "annularEffectivePulseCountMean", ...
    "annularEffectivePulseClusterCountMean", ...
    "roiRaSpotEffectivePulseCountMean", ...
    "roiRaSpotEffectivePulseClusterCountMean", ...
    "roiRaLineEffectivePulseCountMean", ...
    "roiRaLineEffectivePulseClusterCountMean", ...
    "roiRaAnnularEffectivePulseCountMean", ...
    "roiRaAnnularEffectivePulseClusterCountMean", ...
    "roiSpSpotEffectivePulseCountMean", ...
    "roiSpSpotEffectivePulseClusterCountMean", ...
    "roiSpLineEffectivePulseCountMean", ...
    "roiSpLineEffectivePulseClusterCountMean", ...
    "roiSpAnnularEffectivePulseCountMean", ...
    "roiSpAnnularEffectivePulseClusterCountMean");
fprintf("Saved processed data: %s\n", outputFile);

function [meanPulseCount, meanPulseClusterCount, meanDetectionRange] = processResultFile( ...
        resultFile, detectVariable, pulseVariable, columnIndex)
    meanDetectionRange = nan;
    if ~isfile(resultFile)
        warning("Missing result: %s", resultFile);
        meanPulseCount = nan;
        meanPulseClusterCount = nan;
        return;
    end

    result = load(resultFile, detectVariable, pulseVariable);
    if ~isfield(result, detectVariable) || ~isfield(result, pulseVariable)
        warning("Missing expected variables in: %s", resultFile);
        meanPulseCount = nan;
        meanPulseClusterCount = nan;
        return;
    end

    detect = selectResultColumn(result.(detectVariable), ...
        columnIndex, detectVariable, resultFile);
    effectivePulses = selectResultColumn(result.(pulseVariable), ...
        columnIndex, pulseVariable, resultFile);
    [meanPulseCount, meanPulseClusterCount] = ...
        summarizeSuccessfulWindows(detect, effectivePulses, resultFile);
    successfulWindows = effectivePulses(detect ~= 0);
    if ~isempty(successfulWindows)
        meanDetectionRange = mean(recoverDetectionRanges(successfulWindows, resultFile));
    end
end

function ranges = recoverDetectionRanges(windows, resultFile)
    % Raw files do not retain target positions. Recover the 0.5 m range bin
    % from deterministic backscatter at the reporting pulse (window end),
    % NOT the average distance across pulses in the window. Physical
    % constants and quadrature match all nine MC_*_snr beam implementations.
    % The original continuous distance is only recoverable to +/-0.25 m.
    persistent rangeGrid backscatterGrid
    if isempty(rangeGrid)
        rangeGrid = (0:0.5:2000).';
        h = 6.62607015e-34;
        c = 299792458;
        transmittedPhotons = 150e-6 / (h * c / 1550e-9);
        efficiency = 0.80 * 0.80;
        receiverArea = pi * (0.0508 / 2)^2;
        gateRange = c * 50e-9 / 2;
        offsets = linspace(-gateRange/2, gateRange/2, 101);
        r = rangeGrid + offsets;
        r(r <= 0) = 1e-6;
        integrand = 0.3e-6 .* exp(-2 .* 1.5e-5 .* r) ./ r.^2;
        backscatterGrid = transmittedPhotons .* efficiency .* receiverArea .* ...
            trapz(offsets, integrand, 2);
    end
    if ~isfield(windows, "backscatter_photons") || ...
            any(arrayfun(@(w) isempty(w.backscatter_photons), windows))
        error("Missing reporting-pulse backscatter in: %s", resultFile);
    end
    values = arrayfun(@(w) w.backscatter_photons(end), windows);
    indices = interp1(flipud(backscatterGrid), ...
        flipud((1:numel(rangeGrid)).'), values, "nearest", nan);
    if any(~isfinite(indices))
        error("Backscatter outside the distance lookup table: %s", resultFile);
    end
    expected = backscatterGrid(indices);
    if any(abs(values(:) - expected(:)) > 1e-10 .* abs(expected(:)))
        error("Backscatter does not match the simulation distance table: %s", resultFile);
    end
    ranges = rangeGrid(indices);
end

function resultFile = eaResultFileName( ...
        resultDir, beamType, diameterMrad, innerDiameterMrad)
    diameterText = num2str(diameterMrad, "%g");
    if beamType == "POINT"
        fileName = "MC_1par_EA_POINT_D" + diameterText + "mrad.mat";
    else
        innerDiameterText = num2str(innerDiameterMrad, "%g");
        fileName = "MC_1par_EA_" + beamType + "_D" + ...
            diameterText + "d" + innerDiameterText + "mrad.mat";
    end
    resultFile = fullfile(resultDir, fileName);
end

function resultFile = roiResultFileName( ...
        resultDir, beamType, trackType, diameterMrad, innerDiameterMrad)
    diameterText = num2str(diameterMrad, "%g");
    if beamType == "POINT"
        fileName = "MC_1par_ROI_NEW_POINT_" + trackType + ...
            "_D" + diameterText + "mrad.mat";
    else
        innerDiameterText = num2str(innerDiameterMrad, "%g");
        fileName = "MC_1par_ROI_NEW_" + beamType + "_" + trackType + ...
            "_D" + diameterText + "d" + innerDiameterText + "mrad.mat";
    end
    resultFile = fullfile(resultDir, fileName);
end

function values = selectResultColumn( ...
        values, columnIndex, variableName, resultFile)
    if ~ismatrix(values) || size(values, 2) < columnIndex
        error("Variable %s in %s does not contain column %d.", ...
            variableName, resultFile, columnIndex);
    end
    values = values(:, columnIndex);
end

function [meanPulseCount, meanPulseClusterCount] = ...
        summarizeSuccessfulWindows(detect, effectivePulses, resultFile)
    detect = detect(:);
    effectivePulses = effectivePulses(:);
    if numel(detect) ~= numel(effectivePulses)
        error("Detection and effective-pulse arrays differ in size: %s", ...
            resultFile);
    end
    if ~isstruct(effectivePulses) || ...
            ~isfield(effectivePulses, "signal_photons")
        error("Invalid effective-pulse structure in: %s", resultFile);
    end

    successfulWindows = effectivePulses(detect ~= 0);
    if isempty(successfulWindows)
        meanPulseCount = nan;
        meanPulseClusterCount = nan;
        return;
    end

    signalSeries = {successfulWindows.signal_photons};
    pulseCounts = cellfun(@countPositivePulses, signalSeries);
    pulseClusterCounts = cellfun(@countPositiveClusters, signalSeries);
    meanPulseCount = mean(pulseCounts);
    meanPulseClusterCount = mean(pulseClusterCounts);
end

function count = countPositivePulses(signalPhotons)
    hasTargetReturn = isfinite(signalPhotons(:)) & signalPhotons(:) > 0;
    count = nnz(hasTargetReturn);
end

function count = countPositiveClusters(signalPhotons)
    hasTargetReturn = isfinite(signalPhotons(:)) & signalPhotons(:) > 0;
    clusterStarts = hasTargetReturn & [true; ~hasTargetReturn(1:end-1)];
    count = nnz(clusterStarts);
end

function beamCode = beamCodeForType(beamType)
    switch beamType
        case "POINT"
            beamCode = "P";
        case "LINE"
            beamCode = "L";
        case "RING"
            beamCode = "R";
        otherwise
            error("Unsupported beam type: %s", beamType);
    end
end

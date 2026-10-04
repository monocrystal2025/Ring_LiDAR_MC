function cfg = multiuav_default_config()
%MULTIUAV_DEFAULT_CONFIG Configuration for multi-UAV LiDAR simulations.

moduleDir = fileparts(mfilename('fullpath'));

cfg = struct();
cfg.version = "1.0";
cfg.seed = 20260818;
cfg.smokeN = 24;
cfg.previewN = 300;
cfg.screenN = 3000;
cfg.formalN = 10000;
cfg.bootstrapN = 1000;
cfg.outputRoot = fullfile(moduleDir, 'output');
cfg.resultsDir = fullfile(cfg.outputRoot, 'results');
cfg.figuresDir = fullfile(cfg.outputRoot, 'figures');
cfg.checkpointDir = fullfile(cfg.outputRoot, 'checkpoints');

cfg.execution = struct( ...
    'useParallel', false, ...
    'showProgress', true, ...
    'force', false, ...
    'parallelBatchSize', 12, ...
    'motionUpdate_s', 0.02);

cfg.system = struct();
cfg.system.maxRange_m = 2000;
cfg.system.minimumEntryRange_m = 100;
cfg.system.roiHalfAngle_deg = 30;
cfg.system.prf_Hz = 5000;
cfg.system.scanRate_rad_s = 2 * pi;
cfg.system.pulseEnergy_J = 150e-6;
cfg.system.wavelength_m = 1550e-9;
cfg.system.receiverDiameter_m = 0.0508;
cfg.system.opticalTransmission = 0.80;
cfg.system.quantumEfficiency = 0.80;
cfg.system.rangeGate_s = 50e-9;
cfg.system.filterBandwidth_nm = 10;
cfg.system.darkCountRate_Hz = 400;
cfg.system.extinction_m_inv = 1.5e-5;
cfg.system.backscatter_m_inv_sr_inv = 0.3e-6;
cfg.system.skyRadiance_W_m2_sr_nm = 5e-9;
cfg.system.targetWidth_m = 0.30;
cfg.system.targetReflectivity = 0.50;
cfg.system.targetSpeed_m_s = 30;
cfg.system.snrThreshold = 2;
cfg.system.gaussianCutoffWidths = 3;
cfg.system.targetSampleGrid = 7;
cfg.system.noiseLutStep_m = 0.5;
cfg.system.maxObservation_s = 60;

cfg.scene = struct();
cfg.scene.maxTargets = 16;
cfg.scene.targetCounts = [1, 2, 4, 8, 16];
cfg.scene.entrySpan_s = 5;
cfg.scene.turnCorrelation_s = 2;
cfg.scene.turnRateRms_deg_s = 10;
cfg.scene.turnRateLimit_deg_s = 30;
cfg.scene.formationSpacing_m = 15;
cfg.scene.formationSpacingList_m = [5, 15, 30];
cfg.scene.scenarioNames = ["async_maneuver", ...
    "sync_dispersed", "sync_formation"];
cfg.scene.timelyThresholds_s = [1, 2, 5, 10, 30, 60];

cfg.beams = struct( ...
    'name', {"point", "line", "ring"}, ...
    'label', {"Spot", "Line", "Annular"}, ...
    'wd_rad', {NaN, 1e-3, 1e-3}, ...
    'referenceWD_mrad', {55, 195, 125});
cfg.trajectories = ["RA", "SP"];
cfg.screenWD_mrad = 25:25:250;
cfg.formalWD_mrad = [50, 75, 100, 125, 150, 175, 195, 225];
cfg.previewWD_mrad = [125, 195];
cfg.manuscriptWD_mrad = 5:5:500;
cfg.optimization = struct('cvFolds', 5, 'threshold_s', 5);

cfg.constants = struct();
cfg.constants.speedOfLight_m_s = 299792458;
cfg.constants.rangeResolution_m = ...
    cfg.constants.speedOfLight_m_s * cfg.system.rangeGate_s / 2;
end

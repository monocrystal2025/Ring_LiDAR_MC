function cfg = config(mode)
%CONFIG Shared, explicit parameters for the Discussion experiments.
arguments
    mode (1,1) string {mustBeMember(mode,["smoke","preview","paper"])} = "preview"
end
cfg = struct;
cfg.version = "discussion-1.1";
cfg.mode = mode;
cfg.seed = 20260908;
cfg.root = fileparts(fileparts(mfilename('fullpath')));
cfg.outputRoot = fullfile(cfg.root,"discussion_results",mode);
cfg.scanLibrary = "G:/BeamVEC_NEW";
cfg.scanSource = "library"; % Explicit portability alternative: "generated".
cfg.N = 128;
cfg.bootstrapN = 500;
cfg.widthGrid_mrad = 30:30:300;
cfg.maxTime_s = 15;
cfg.targetProbability = 0.80;
cfg.useParallel = false;
cfg.blockSize = 25000;
cfg.collectAll = true;
cfg.heading_deg = NaN; % NaN: original isotropic inward velocity distribution.
cfg.fixedWindow_s = NaN; % NaN: ceil(f*wD/Omega), as in the current MC code.
cfg.lineRotation_deg = 0;
cfg.jitterRms_rad = 0; % Independent pointing jitter per axis, per pulse.
cfg.ringGap_deg = 0;
cfg.gapRedistribute = false; % false: remove photons, do not renormalize.
cfg.returnScale = 1;
cfg.receiverMode = "matched"; % matched, envelope, or common.
cfg.commonFov_rad = 0.25;
cfg.phaseMode = "legacy"; % "cycle" samples the entire repeated scan cycle.
cfg.localRange_m = 2000;
cfg.localRelativeSpeed_rad_s = 2*pi;
cfg.localImpactScale = 1;
cfg.system = struct('maxRange_m',2000,'prf_Hz',5000, ...
    'scanRate_rad_s',2*pi,'targetSpeed_m_s',30, ...
    'beamWidth_rad',1e-3,'pulseEnergy_J',150e-6, ...
    'wavelength_m',1550e-9,'receiverDiameter_m',0.0508, ...
    'opticalTransmission',0.80,'quantumEfficiency',0.80, ...
    'rangeGate_s',50e-9,'filterBandwidth_nm',10, ...
    'darkCountRate_Hz',400,'extinction_m_inv',1.5e-5, ...
    'backscatter_m_inv_sr_inv',0.3e-6, ...
    'skyRadiance_W_m2_sr_nm',5e-9,'targetWidth_m',0.297, ...
    'targetReflectivity',0.50,'snrThreshold',2, ...
    'targetSampleGrid',7,'noiseLutStep_m',0.5);
switch mode
    case "smoke"
        cfg.N = 12;
        cfg.bootstrapN = 80;
        cfg.widthGrid_mrad = [90,150,210];
    case "paper"
        cfg.N = 3000;
        cfg.bootstrapN = 2000;
        cfg.widthGrid_mrad = 30:5:300;
    otherwise
        % Preview retains physical pulse rate, target grid, and horizon.
end
end

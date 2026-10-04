function manifest = para_data(options)
%PARA_DATA Generate the four one-at-a-time EA sensitivity sweeps.
%   para_data                       % N=800, 5:5:400 mrad, five values/parameter
%   para_data(N=1000)                % use a new OutputDirectory if N changes
%   para_data(DryRun=true)           % inspect configuration without simulation
%   para_predraw                    % preprocess after data generation
%   para_draw                       % draw the saved plotting snapshot
%
% Defaults: E=150 uJ, rho=0.5, v=30 m/s, alpha=1.5e-5 /m;
% Energy=[30 75 150 300 750] uJ; alpha=[0.5 1 1.5 2.25 4.5]*1e-5 /m.
% beta=alpha/50 throughout. Factors applies to reflectivity and speed only.
% Each raw MAT retains the MC_EA N-by-1 detect_L/R, first_time_L/R,
% first_encounter_time_L/R and effective_pulses_L/R variables. Static files
% use the SAME variable schema, in separately named PARA_STATIC_EA files.
% Extra metadata and first_detection_range_L/R (0.5 m bins) are included.
% Baseline cases and speed-independent STATIC cases share physical files.
%
% Algorithms: MC_func scan/phase convention and init_UAV for moving targets;
% blind_snr_EA uniform-volume stationary probes and one complete scan;
% the original MC_line_snr/MC_ring_snr kernels, with optional physics input.
% Moving windows are 'fixed' (MC_EA); static windows are 'variable'
% (blind_snr_EA). No finite observation cutoff is imposed by default.
% MaxStepIndex is ONLY a smoke-test limit; para_draw rejects truncated runs.
arguments
    options.N (1,1) double {mustBeInteger,mustBePositive} = 3000
    options.Factors (1,:) double {mustBePositive,mustBeFinite} = [0.5 0.75 1 1.25 1.5]
    options.EnergyMicroJ (1,:) double {mustBePositive,mustBeFinite} = [30 75 150 300 750]
    options.AlphaScaled (1,:) double {mustBePositive,mustBeFinite} = [0.5 1.0 1.5 2.25 4.5]
    options.WidthsMrad (1,:) double {mustBePositive,mustBeFinite} = 5:5:400
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 20260923
    options.OutputDirectory (1,1) string = ""
    options.BeamPathDirectory (1,1) string = "G:/BeamVEC_NEW"
    options.UseParallel (1,1) logical = true
    options.DryRun (1,1) logical = false
    options.MaxStepIndex (1,1) double {mustBeNonnegative} = Inf
end
root = fileparts(mfilename('fullpath'));
if options.OutputDirectory == ""
    options.OutputDirectory = fullfile(root,'para_sensative_N3000_range2');
end
assert(issorted(options.Factors) && numel(unique(options.Factors)) == numel(options.Factors) ...
    && any(options.Factors == 1) && max(options.Factors) <= 2, ...
    'para_data:Factors','Factors must be unique, increasing, include 1 and be <=2 (rho<=1).');
assert(issorted(options.EnergyMicroJ) && numel(unique(options.EnergyMicroJ)) == numel(options.EnergyMicroJ) ...
    && any(options.EnergyMicroJ == 150), ...
    'para_data:Energy','EnergyMicroJ must be unique, increasing and include 150 uJ.');
assert(issorted(options.AlphaScaled) && numel(unique(options.AlphaScaled)) == numel(options.AlphaScaled) ...
    && any(options.AlphaScaled == 1.5), ...
    'para_data:Alpha','AlphaScaled must be unique, increasing and include 1.5 (units: 1e-5 /m).');
assert(issorted(options.WidthsMrad) && numel(unique(options.WidthsMrad)) == numel(options.WidthsMrad), ...
    'para_data:Widths','WidthsMrad must be unique and increasing.');
assert(isinf(options.MaxStepIndex) || options.MaxStepIndex == fix(options.MaxStepIndex), ...
    'para_data:Horizon','MaxStepIndex must be an integer or Inf.');
assert(options.Seed <= 2^32-1,'para_data:Seed','Seed must fit uint32.');
manifest = makeManifest(options,root);
if options.DryRun
    return
end
% Fail before starting an expensive sweep if any required trajectory is absent.
assert(exist('readNPY','file') == 2,'para_data:ReadNPY','readNPY must be on the MATLAB path.');
for d = 1:numel(manifest.widthsMrad)
    assert(isfile(manifest.scanFiles(d)),'para_data:MissingScan', ...
        'Missing MC_EA scan trajectory: %s',manifest.scanFiles(d));
end
if ~isfolder(options.OutputDirectory)
    mkdir(options.OutputDirectory);
end
manifestFile = fullfile(options.OutputDirectory,'para_manifest.mat');
if isfile(manifestFile)
    saved = load(manifestFile,'manifest');
    assert(isequaln(saved.manifest.signature,manifest.signature), ...
        'para_data:DifferentRun', ...
        'Configuration or kernel source changed. Choose a new OutputDirectory.');
end
save(manifestFile,'manifest','-v7.3');
previousRng = rng;
rngCleanup = onCleanup(@() rng(previousRng));
rng(options.Seed,'twister');
[positions,velocity30] = init_UAV(manifest.R,1,30,options.N);
phase = rand(options.N,1);
azimuth = 2*pi*rand(options.N,1);
cosPolar = rand(options.N,1);
radius = manifest.R*nthroot(max(rand(options.N,1),realmin),3);
staticPositions = radius .* [sqrt(1-cosPolar.^2).*cos(azimuth), ...
    sqrt(1-cosPolar.^2).*sin(azimuth),cosPolar];
useParallel = options.UseParallel && license('test','Distrib_Computing_Toolbox') ...
    && ~isempty(ver('parallel'));
if useParallel && isempty(gcp('nocreate'))
    try
        parpool('threads');
    catch
        parpool('local');
    end
end
for d = numel(manifest.widthsMrad):-1:1
    wD = manifest.widthsMrad(d)*1e-3;
    scan = loadScan(manifest.scanFiles(d),manifest.f,manifest.omega);
    % Exactly MC_func's phase support (NOT the longer full forward/reverse cycle).
    initialIndices = floor(phase*scan.phaseLength)+1;
    for j = 1:numel(manifest.cases)
        oneCase = manifest.cases(j);
        for b = 1:2
            beam = manifest.beams(b);
            blockSize = mc_ea_block_size(wD,char(beam));
            for mode = 1:2
                isStatic = mode == 2;
                if isStatic
                    outputFile = oneCase.staticFiles(d,b);
                    P = staticPositions;
                    V = zeros(options.N,3);
                    indices = ones(options.N,1);
                    horizon = min(scan.cycleLength-1,options.MaxStepIndex);
                    windowMode = 'variable';
                else
                    outputFile = oneCase.dynamicFiles(d,b);
                    P = positions;
                    V = velocity30*(oneCase.speed/30);
                    indices = initialIndices;
                    horizon = options.MaxStepIndex;
                    windowMode = 'fixed';
                end
                meta = struct('schemaVersion',manifest.schemaVersion,'N',options.N, ...
                    'seed',options.Seed,'physics',oneCase.physics,'speed',oneCase.speed, ...
                    'widthMrad',manifest.widthsMrad(d),'beam',beam, ...
                    'isStatic',isStatic,'windowMode',windowMode,'maxStepIndex',horizon);
                if isStatic
                    meta.speed = 0;
                end
                if isfile(outputFile)
                    existing = load(outputFile,'meta');
                    assert(isfield(existing,'meta') && isequaln(existing.meta,meta), ...
                        'para_data:StaleFile','Incompatible saved result: %s',outputFile);
                    continue
                end
                fprintf('%s: %s, D=%g mrad, static=%d, N=%d\n', ...
                    oneCase.id,beam,manifest.widthsMrad(d),isStatic,options.N);
                raw = simulatePopulation(manifest,scan,P,V,indices,beam,wD, ...
                    blockSize,horizon,windowMode,oneCase.physics,useParallel);
                raw.meta = meta;
                if ~isfolder(fileparts(outputFile))
                    mkdir(fileparts(outputFile));
                end
                % Commit a file only after the full population completed.
                temporaryFile = outputFile + ".partial.mat";
                save(temporaryFile,'-struct','raw','-v7.3');
                movefile(temporaryFile,outputFile);
            end
        end
    end
end
manifest.complete = true;
save(manifestFile,'manifest','-v7.3');
fprintf('Saved sensitivity data: %s\nRun para_predraw, then para_draw.\n',manifestFile);
end

function manifest = makeManifest(options,root)
manifest = struct('schemaVersion',1,'N',options.N,'f',5000,'omega',2*pi, ...
    'R',2000,'smallWidthRad',1e-3,'baseline',[150e-6 0.5 30 1.5e-5], ...
    'alphaBetaRatio',50,'widthsMrad',options.WidthsMrad(:), ...
    'factors',options.Factors,'beams',["line","ring"],'complete',false, ...
    'truncated',isfinite(options.MaxStepIndex));
manifest.parameters = ["energy","reflectivity","speed","alpha"];
manifest.labels = ["Pulse energy (\muJ)","Reflectivity","Speed (m/s)","\alpha (10^{-5} m^{-1})"];
manifest.displayScale = [1e6 1 1 1e5];
% Scale from the baseline so shared baseline cases remain exactly identical.
manifest.parameterValues = {manifest.baseline(1)*(options.EnergyMicroJ/150), ...
    manifest.baseline(2)*options.Factors,manifest.baseline(3)*options.Factors, ...
    manifest.baseline(4)*(options.AlphaScaled/1.5)};
manifest.outputDirectory = options.OutputDirectory;
manifest.scanFiles = strings(numel(manifest.widthsMrad),1);
for d = 1:numel(manifest.widthsMrad)
    manifest.scanFiles(d) = fullfile(options.BeamPathDirectory, ...
        sprintf('beam_vector_fast_5000Hz_h3000m_D%gu.npy',manifest.widthsMrad(d)*1000));
end
sourceNames = ["para_data.m","MC_line_snr.m","MC_ring_snr.m", ...
    "init_UAV.m","append_and_check_photon_window.m", ...
    "empty_detection_window_photons.m","mc_ea_block_size.m"];
sourceText = cell(size(sourceNames));
for k = 1:numel(sourceNames)
    sourceText{k} = fileread(fullfile(root,sourceNames(k)));
end
signatureOptions = rmfield(options,{'DryRun','UseParallel'});
manifest.signature = struct('options',signatureOptions,'sourceNames',sourceNames, ...
    'sourceText',{sourceText});
prototype = struct('id',"",'parameterIndex',0,'factorIndex',0,'value',0, ...
    'physics',[],'speed',30,'dynamicFiles',[],'staticFiles',[]);
manifest.cases = repmat(prototype,sum(cellfun(@numel,manifest.parameterValues)),1);
j = 0;
for p = 1:4
    for k = 1:numel(manifest.parameterValues{p})
        j = j+1;
        values = manifest.baseline;
        values(p) = manifest.parameterValues{p}(k);
        id = manifest.parameters(p)+sprintf('_s%02d',k);
        dynamicId = id;
        if values(p) == manifest.baseline(p)
            dynamicId = "baseline";
        end
        staticId = dynamicId;
        if p == 3
            staticId = "baseline";
        end
        item = prototype;
        item.id = id;
        item.parameterIndex = p;
        item.factorIndex = k;
        item.value = values(p);
        item.physics = [values(1:2),values(4),values(4)/manifest.alphaBetaRatio];
        item.speed = values(3);
        item.dynamicFiles = strings(numel(manifest.widthsMrad),2);
        item.staticFiles = item.dynamicFiles;
        for d = 1:numel(manifest.widthsMrad)
            for b = 1:2
                suffix = sprintf('%s_D%gd1mrad.mat',upper(manifest.beams(b)),manifest.widthsMrad(d));
                item.dynamicFiles(d,b) = fullfile(options.OutputDirectory,'raw', ...
                    dynamicId,"PARA_MC_1par_EA_"+suffix);
                item.staticFiles(d,b) = fullfile(options.OutputDirectory,'raw', ...
                    staticId,"PARA_STATIC_MC_1par_EA_"+suffix);
            end
        end
        manifest.cases(j) = item;
    end
end
end

function scan = loadScan(file,f,omega)
scan.beam = readNPY(char(file));
validateattributes(scan.beam,{'numeric'},{'2d','ncols',3,'nonempty','finite'});
scan.path1 = boundaryCircle(scan.beam(1,:),f,omega);
scan.path2 = boundaryCircle(scan.beam(end,:),f,omega);
scan.phaseLength = size(scan.path1,1)+size(scan.beam,1)+size(scan.path2,1);
scan.cycleLength = scan.phaseLength+size(scan.beam,1);
end

function circle = boundaryCircle(point,f,omega)
theta0 = atan2(point(2),point(1));
r0 = hypot(point(1),point(2));
theta = (theta0:omega/f/r0:theta0+2*pi).';
circle = [r0*cos(theta),r0*sin(theta),ones(numel(theta),1)*point(3)];
end

function raw = simulatePopulation(cfg,scan,P,V,indices,beam,wD,blockSize,horizon,windowMode,physics,useParallel)
n = size(P,1);
detect = false(n,1);
firstPulse = nan(n,1);
encounterPulse = nan(n,1);
windows = repmat(empty_detection_window_photons(),n,1);
if beam == "line"
    kernel = @MC_line_snr;
    code = 'L';
else
    kernel = @MC_ring_snr;
    code = 'R';
end
workers = 0;
f = cfg.f;
omega = cfg.omega;
smallWidth = cfg.smallWidthRad;
regionRange = cfg.R;
if useParallel
    workers = Inf;
end
parfor (i = 1:n, workers)
    [detect(i),firstPulse(i),encounterPulse(i),windows(i)] = kernel( ...
        wD,smallWidth,scan.path1,scan.beam,scan.path2,f,omega, ...
        P(i,:),V(i,:),indices(i),regionRange,blockSize,horizon,windowMode,physics);
end
firstTime = firstPulse/cfg.f;
ranges = nan(n,1);
% Same bins recovered from reporting-pulse backscatter in advantage_data.m.
% Direct positions avoid incorrectly using its fixed-energy atmosphere LUT.
ranges(detect) = 0.5*round(vecnorm(P(detect,:)+V(detect,:).*firstTime(detect),2,2)/0.5);
raw = struct();
raw.(['detect_' code]) = logical(detect);
raw.(['first_time_' code]) = firstTime;
raw.(['first_encounter_time_' code]) = encounterPulse/cfg.f;
raw.(['effective_pulses_' code]) = windows;
raw.(['first_detection_range_' code]) = ranges;
end

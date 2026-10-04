function data = parasensitive_data(options)
%PARASENSITIVE_DATA Four one-at-a-time EA (paper: ROE) hemisphere sweeps.
%   parasensitive_data                         % generate fresh data, N=200
%   parasensitive_data(N=24,WidthsMrad=[5 25 50 100 150 200 300 400],Preview=true)
%   parasensitive_predraw; parasensitive_draw
% All samples and metadata are saved in parasensitive_data.mat beside this
% file. Every call computes fresh results and overwrites the output MAT file.
% Paper baseline: 150 uJ, 5 kHz, 30 m/s, alpha=1.5e-5 /m; beta=alpha/50.
% Moving targets: init_UAV, fixed detection window, run until region exit.
% Static blindness: uniform hemisphere volume, one entire scan cycle,
% variable detection window (same definition as blind_snr_EA/para_data).
% Count strictly positive signal pulses in the FIRST successful window.
% Zero divergence is singular; only positive divergences are simulated.
arguments
    options.N (1,1) double {mustBeInteger,mustBePositive} = 200
    options.WidthsMrad (1,:) double {mustBePositive,mustBeFinite} = 5:5:400
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 20260928
    options.UseParallel (1,1) logical = true
    options.Preview (1,1) logical = false
end
assert(numel(options.WidthsMrad)>=2 && all(diff(options.WidthsMrad)>0) ...
    && max(options.WidthsMrad)<1000, 'parasensitive_data:Grid', ...
    'Use at least two strictly increasing widths below 1000 mrad.');
assert(options.Seed<=2^32-1,'parasensitive_data:Seed','Seed must fit uint32.');
root = fileparts(mfilename('fullpath'));
file = fullfile(root,'parasensitive_data.mat');
config = rmfield(options,'UseParallel');
% MC_EA/blind_snr_EA use 5 kHz; MC_EA/init_UAV use 30 m/s.
% MC_line_snr/MC_ring_snr default physics: [150e-6 0.5 1.5e-5 0.3e-6].
config.baseline = [150e-6 5000 30 1.5e-5];
config.values = {[100 150 200 250 300]*1e-6, [2 4 6 8 10]*1e3, ...
    [20 25 30 35 40], [0.1 0.3 0.5 0.7 0.9]*1e-4};
config.parameters = ["energy","prf","speed","alpha"];
config.labels = ["Pulse energy (\muJ)","Repetition rate (kHz)", ...
    "Target speed (m/s)","\alpha (10^{-4} m^{-1})"];
config.displayScale = [1e6 1e-3 1 1e4];
config.R = 2000;
config.omega = 2*pi;
config.smallWidthRad = 1e-3;
config.lineShortPrefilterHalfWidthRad = 1.5*config.smallWidthRad;
config.lidarRatio = 50;
config.beams = ["line","ring"];
sourceNames = ["parasensitive_data.m","MC_line_snr.m","MC_ring_snr.m", ...
    "init_UAV.m","append_and_check_photon_window.m", ...
    "empty_detection_window_photons.m","mc_ea_block_size.m"];
config.sourceNames = sourceNames;
config.sourceText = cell(size(sourceNames));
for k = 1:numel(sourceNames)
    config.sourceText{k} = fileread(fullfile(root,sourceNames(k)));
end
% Use exactly the same kernel as MC_EA -> MC_func: 3*(d/2)=1.5*d.
lineSource = config.sourceText{2};
assert(~isempty(regexp(lineSource,'GAUSSIAN_CUTOFF_WIDTH\s*=\s*3\s*;', 'once')) ...
    && ~isempty(regexp(lineSource,'gaussian_width\s*=\s*fasan_d\s*/\s*2\s*;', 'once')) ...
    && ~isempty(regexp(lineSource,'short_prefilter_half_width\s*=\s*GAUSSIAN_CUTOFF_WIDTH\s*\*\s*gaussian_width\s*;', 'once')), ...
    'parasensitive_data:LineCutoff','MC_line_snr must use the paper short-side prefilter +/-1.5d.');
[cases,caseMap] = makeCases(config);
nWidths = numel(options.WidthsMrad);
nCases = size(cases,1);
data = struct('schemaVersion',1,'config',config,'cases',cases, ...
    'caseMap',caseMap,'dynamic',{cell(nWidths,nCases,2)}, ...
    'static',{cell(nWidths,nCases,2)},'complete',false, ...
    'created',datetime('now'),'elapsedSeconds',0);
fprintf('Simulation: N=%d, preview=%d, %d widths, %d cases.\n', ...
    options.N,options.Preview,nWidths,nCases);
fprintf('Fixed paper baseline: E=150 uJ, f=5 kHz, v=30 m/s, alpha=1.5e-5 /m.\n');
fprintf('w_D is beam divergence; scanning largest to smallest. Each sweep changes only its named parameter.\n');
previousRng = rng;
rngCleanup = onCleanup(@() rng(previousRng));
rng(options.Seed,'twister');
[positions,unitVelocity] = init_UAV(config.R,1,1,options.N);
phase = rand(options.N,1);
azimuth = 2*pi*rand(options.N,1);
cosPolar = rand(options.N,1);
radius = config.R*nthroot(max(rand(options.N,1),realmin),3);
staticPositions = radius.*[sqrt(1-cosPolar.^2).*cos(azimuth), ...
    sqrt(1-cosPolar.^2).*sin(azimuth),cosPolar];
data.population = struct('entryPositions',positions,'unitVelocity',unitVelocity, ...
    'phaseFraction',phase,'staticPositions',staticPositions);
workers = 0;
if options.UseParallel && license('test','Distrib_Computing_Toolbox')
    if isempty(gcp('nocreate'))
        parpool('threads');
    end
    workers = Inf;
end
timer = tic;
frequencies = unique(cases(:,2));
for d = nWidths:-1:1
    width = options.WidthsMrad(d)*1e-3;
    for q = 1:numel(frequencies)
        f = frequencies(q);
        members = find(cases(:,2)==f).';
        scan = makeScan(width,f,config.omega);
        indices = floor(phase*scan.phaseLength)+1;
        for j = members
            oneCase = cases(j,:);
            physics = [oneCase(1),0.5,oneCase(4),oneCase(4)/config.lidarRatio];
            for b = 1:2
                if isempty(data.dynamic{d,j,b})
                    data.dynamic{d,j,b} = simulate(config,scan,positions, ...
                        unitVelocity*oneCase(3),indices,b,width,physics,Inf,'fixed',workers);
                end
                if isempty(data.static{d,j,b})
                    % Static geometry and physics do not depend on speed.
                    donors = find(all(cases(:,[1 2 4])==oneCase([1 2 4]),2));
                    available = donors(~cellfun(@isempty,data.static(d,donors,b)));
                    if isempty(available)
                        data.static{d,j,b} = simulate(config,scan,staticPositions, ...
                            zeros(options.N,3),ones(options.N,1),b,width,physics, ...
                            scan.cycleLength-1,'variable',workers);
                    else
                        data.static{d,j,b} = data.static{d,available(1),b};
                    end
                end
            end
            [~,sweeps] = find(caseMap==j);
            sweepLabel = strjoin(config.parameters(unique(sweeps)), '/');
            fprintf('ROE [%s sweep]: w_D=%g mrad, E=%g uJ, f=%g kHz, v=%g m/s, alpha=%g /m\n', ...
                char(sweepLabel),width*1e3,oneCase(1)*1e6,f/1e3,oneCase(3),oneCase(4));
        end
    end
    data.elapsedSeconds = toc(timer);
    fprintf('Progress: %d/%d widths complete, %.1f s elapsed.\n', ...
        nWidths-d+1,nWidths,data.elapsedSeconds);
end
data.complete = true;
data.elapsedSeconds = toc(timer);
save(file,'data','-v7.3');
fprintf('Saved %s (N=%d, preview=%d).\n',file,options.N,options.Preview);
end

function [cases,caseMap] = makeCases(config)
allCases = repmat(config.baseline,20,1);
for p = 1:4
    allCases((p-1)*5+(1:5),p) = config.values{p}.';
end
[cases,~,mapping] = unique(allCases,'rows','stable');
caseMap = reshape(mapping,5,4);
end

function raw = simulate(cfg,scan,P,V,indices,b,width,physics,horizon,windowMode,workers)
n = size(P,1);
detect = false(n,1);
firstPulse = nan(n,1);
encounterPulse = nan(n,1);
pulseCount = nan(n,1);
clusterCount = nan(n,1);
if b==1
    kernel = @MC_line_snr;
else
    kernel = @MC_ring_snr;
end
blockSize = mc_ea_block_size(width,char(cfg.beams(b)));
smallWidth = cfg.smallWidthRad;
omega = cfg.omega;
regionRange = cfg.R;
f = scan.f;
path1 = scan.path1;
beam = scan.beam;
path2 = scan.path2;
parfor (i = 1:n,workers)
    [detect(i),firstPulse(i),encounterPulse(i),window] = kernel( ...
        width,smallWidth,path1,beam,path2, ...
        f,omega,P(i,:),V(i,:),indices(i),regionRange, ...
        blockSize,horizon,windowMode,physics);
    if detect(i)
        positive = isfinite(window.signal_photons(:)) & window.signal_photons(:)>0;
        pulseCount(i) = nnz(positive);
        clusterCount(i) = nnz(positive & [true;~positive(1:end-1)]);
    end
end
raw = struct('detect',detect,'firstTime',firstPulse/f, ...
    'firstEncounterTime',encounterPulse/f,'effectivePulseCount',pulseCount, ...
    'effectiveClusterCount',clusterCount);
end

function scan = makeScan(width,f,omega)
% Same spherical spiral and ordering as flipud(lidarpath(...)), without UI
% or external NPY files. Regenerate on the actual 1/f pulse-time grid.
spacing = 2*tan(width/2);
endpoint = min(1,spacing*sqrt(1-spacing^2/4));
phiLower = acos(endpoint)+0.5*asin(endpoint);
phiUpper = asin(endpoint)/2;
k = spacing/(2*pi);
arc = integral(@(phi) sqrt(k^2+sin(phi).^2),phiUpper,phiLower, ...
    'RelTol',1e-11,'AbsTol',1e-12)/k;
pulseTimes = (0:floor(arc/omega*f)).'/f;
solution = ode45(@(~,phi) -k./sqrt(k^2+sin(phi).^2),[0 arc], ...
    phiLower,odeset('RelTol',1e-9,'AbsTol',1e-11));
phi = deval(solution,omega*pulseTimes).';
phi = min(phiLower,max(phiUpper,phi));
theta = (phiLower-phi)/k;
beam = [sin(phi).*cos(theta),sin(phi).*sin(theta),cos(phi)];
scan.beam = flipud(beam./vecnorm(beam,2,2));
scan.path1 = boundaryCircle(scan.beam(1,:),f,omega);
scan.path2 = boundaryCircle(scan.beam(end,:),f,omega);
scan.phaseLength = size(scan.path1,1)+size(scan.beam,1)+size(scan.path2,1);
scan.cycleLength = scan.phaseLength+size(scan.beam,1);
scan.f = f;
end

function circle = boundaryCircle(point,f,omega)
theta0 = atan2(point(2),point(1));
radius = hypot(point(1),point(2));
theta = (theta0:omega/f/radius:theta0+2*pi).';
circle = [radius*cos(theta),radius*sin(theta),repmat(point(3),numel(theta),1)];
end

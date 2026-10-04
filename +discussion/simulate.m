function result = simulate(cfg,beam,bank)
%SIMULATE Paired finite-horizon ROE experiment using the current scan library.
% Photon statistics are EXPECTATIONS. This is an SNR-feasibility ensemble,
% not a measured or Poisson-realization detection probability at fixed Pfa.
arguments
    cfg (1,1) struct
    beam (1,1) struct
    bank (1,1) struct
end
assert(cfg.maxTime_s > 0 && cfg.system.prf_Hz > 0);
assert(cfg.ringGap_deg >= 0 && cfg.ringGap_deg < 360);
scan = discussion.scan(cfg,beam.wD_rad);
model = discussion.photonModel(cfg,beam);
[positions,velocities] = mapBank(cfg,bank);
window = max(1,ceil(cfg.system.prf_Hz*beam.wD_rad/cfg.system.scanRate_rad_s));
if isfinite(cfg.fixedWindow_s)
    window = max(1,ceil(cfg.system.prf_Hz*cfg.fixedWindow_s));
end
raw = cell(bank.N,1);
timer = tic;
if cfg.useParallel
    parfor i = 1:bank.N
        raw{i} = simulateOne(cfg,beam,bank,i,scan,model,positions(i,:),velocities(i,:),window);
    end
else
    for i = 1:bank.N
        raw{i} = simulateOne(cfg,beam,bank,i,scan,model,positions(i,:),velocities(i,:),window);
    end
end
result = struct('beam',beam,'metrics',struct2table(vertcat(raw{:})), ...
    'windowPulses',window,'scanSource',scan.source, ...
    'phaseMode',cfg.phaseMode,'cycleTime_s',scan.length/cfg.system.prf_Hz, ...
    'model',rmfield(model,'backscatter'),'elapsed_s',toc(timer), ...
    'bankSeed',bank.seed,'N',bank.N);
end

function metrics = simulateOne(cfg,beam,bank,i,scan,model,p,v,window)
s = cfg.system;
support = scan.initialSupport;
if cfg.phaseMode == "cycle"
    support = scan.length;
end
initial = floor(bank.u(i,5)*support)+1;
% Analytic exit time keeps misses and early exits in the common denominator.
exitSphere = -2*dot(p,v)/max(dot(v,v),realmin);
exitGround = inf;
if v(3) < 0
    exitGround = -p(3)/v(3);
end
horizon = max(0,min([cfg.maxTime_s,exitSphere,exitGround]));
steps = (0:floor(horizon*s.prf_Hz)).';
time = steps/s.prf_Hz;
target = p+time.*v;
ranges = max(vecnorm(target,2,2),1e-9);
u = target./ranges;
index = mod(initial-1+steps,scan.length)+1;
b = scan.B(index,:);
long = scan.long(index,:);
short = scan.short(index,:);
if cfg.jitterRms_rad > 0
    stream = RandStream('mt19937ar','Seed',mod(bank.seed+i*104729,2^32-1));
    jitter = cfg.jitterRms_rad*randn(stream,numel(steps),2);
    b = b+jitter(:,1).*long+jitter(:,2).*short;
    b = b./vecnorm(b,2,2);
    % Reorthogonalize receiver/transmitter frame after common pointing jitter.
    short = short-sum(short.*b,2).*b;
    short = short./vecnorm(short,2,2);
    long = cross(b,short,2);
end
factor = discussion.illumination(cfg,beam,ranges,sum(u.*b,2), ...
    sum(u.*long,2),sum(u.*short,2));
[signal,noise] = discussion.photons(model,ranges,factor,cfg.returnScale);
metrics = discussion.score(signal,noise,time,window,s.snrThreshold,factor);
end

function [p,v] = mapBank(cfg,bank)
u = bank.u;
azimuth = 2*pi*u(:,1);
cosPhi = u(:,2);
direction = [sqrt(1-cosPhi.^2).*cos(azimuth), ...
    sqrt(1-cosPhi.^2).*sin(azimuth),cosPhi];
p = cfg.system.maxRange_m*direction;
velocityAz = 2*pi*u(:,3);
velocityCos = 2*u(:,4)-1;
unitVelocity = [sqrt(1-velocityCos.^2).*cos(velocityAz), ...
    sqrt(1-velocityCos.^2).*sin(velocityAz),velocityCos];
outward = sum(unitVelocity.*direction,2) > 0;
unitVelocity(outward,:) = -unitVelocity(outward,:);
if isfinite(cfg.heading_deg)
    % Preserve radial incidence (and sphere dwell time), rotate the
    % transverse velocity at the ENTRY LOS. This is not the instantaneous
    % angle at a later target-beam encounter.
    radial = -sum(unitVelocity.*direction,2);
    short = [-sin(azimuth),cos(azimuth),zeros(bank.N,1)];
    long = cross(direction,short,2);
    angle = deg2rad(cfg.heading_deg);
    tangent = cos(angle)*long+sin(angle)*short;
    unitVelocity = -radial.*direction+sqrt(max(0,1-radial.^2)).*tangent;
end
v = cfg.system.targetSpeed_m_s*unitVelocity;
end

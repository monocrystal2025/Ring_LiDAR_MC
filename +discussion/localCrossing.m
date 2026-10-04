function result = localCrossing(cfg,beam,bank,angle_deg)
%LOCALCROSSING Controlled relative target-field crossing in a tangent plane.
% angle=0 is parallel to the line long axis. Relative angular speed is
% prescribed; it must not be mislabeled as UAV speed. Impact parameter and
% sub-pulse phase are paired across beams and angles.
arguments
    cfg (1,1) struct
    beam (1,1) struct
    bank (1,1) struct
    angle_deg (1,1) double
end
model = discussion.photonModel(cfg,beam);
omega = cfg.localRelativeSpeed_rad_s;
extent = beam.wD_rad/2+4*beam.wd_rad;
travelTime = 2*extent/omega;
window = max(1,ceil(cfg.system.prf_Hz*beam.wD_rad/cfg.system.scanRate_rad_s));
if isfinite(cfg.fixedWindow_s)
    window = max(1,ceil(cfg.system.prf_Hz*cfg.fixedWindow_s));
end
raw = cell(bank.N,1);
trace = struct;
for i = 1:bank.N
    time = (0:ceil(travelTime*cfg.system.prf_Hz)).'/cfg.system.prf_Hz;
    distance = omega*(time-travelTime/2+bank.u(i,6)/cfg.system.prf_Hz);
    impact = (2*bank.u(i,7)-1)*beam.wD_rad/2*cfg.localImpactScale;
    x = distance*cosd(angle_deg)-impact*sind(angle_deg);
    y = distance*sind(angle_deg)+impact*cosd(angle_deg);
    theta = hypot(x,y);
    scale = ones(size(theta));
    nz = theta > 0;
    scale(nz) = sin(theta(nz))./theta(nz);
    ranges = repmat(cfg.localRange_m,size(time));
    factor = discussion.illumination(cfg,beam,ranges,cos(theta),x.*scale,y.*scale);
    [signal,noise] = discussion.photons(model,ranges,factor,cfg.returnScale);
    if i == 1
        [raw{i},trace] = discussion.score(signal,noise,time,window, ...
            cfg.system.snrThreshold,factor);
    else
        raw{i} = discussion.score(signal,noise,time,window,cfg.system.snrThreshold,factor);
    end
end
result = struct('beam',beam,'metrics',struct2table(vertcat(raw{:})), ...
    'trace',trace,'angle_deg',angle_deg,'windowPulses',window, ...
    'model',rmfield(model,'backscatter'),'N',bank.N);
end

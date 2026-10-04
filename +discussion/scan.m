function scan = scan(cfg,width_rad)
%SCAN Reproduce the current MC_func ROE cycle and boundary circles.
% The stored library uses 5 kHz and 2*pi rad/s. Other f/Omega ratios are
% interpolated on the SAME geometric cycle. Generated paths are opt-in.
arguments
    cfg (1,1) struct
    width_rad (1,1) double {mustBePositive}
end
persistent cachedKey cachedScan
key = sprintf('%s|%s|%.15g|%.15g|%.15g',cfg.scanSource,cfg.scanLibrary, ...
    width_rad,cfg.system.prf_Hz,cfg.system.scanRate_rad_s);
if isequal(key,cachedKey)
    scan = cachedScan;
    return
end
baseF = 5000;
baseOmega = 2*pi;
if cfg.scanSource == "library"
    filename = fullfile(cfg.scanLibrary,sprintf( ...
        'beam_vector_fast_5000Hz_h3000m_D%gu.npy',round(width_rad*1e6)));
    assert(isfile(filename),'discussion:MissingScan', ...
        'Missing %s. Set cfg.scanSource="generated" explicitly for an analytic path.',filename);
    assert(~isempty(which('readNPY')),'discussion:MissingNpyReader', ...
        'readNPY is required for the existing scan library.');
    spiral = readNPY(filename);
    source = string(filename);
elseif cfg.scanSource == "generated"
    spiral = generateSpiral(width_rad,baseOmega/baseF);
    source = "Analytic spherical spiral; not the stored scan library";
else
    error('discussion:ScanSource','scanSource must be library or generated.');
end
spiral = spiral./vecnorm(spiral,2,2);
path1 = boundaryCircle(spiral(1,:),baseOmega/baseF);
path2 = boundaryCircle(spiral(end,:),baseOmega/baseF);
cycle = [spiral;path2;flipud(spiral);path1];
% Preserve MC_func's original initial-phase support. tmp_4 explicitly
% compares it with uniform sampling over the complete repeated cycle.
initialSupport = size(path1,1)+size(spiral,1)+size(path2,1);
ratio = (cfg.system.scanRate_rad_s/cfg.system.prf_Hz)/(baseOmega/baseF);
if abs(ratio-1) > 1e-12
    query = (0:ratio:size(cycle,1)-eps(size(cycle,1))).';
    old = (0:size(cycle,1)).';
    cycle = interp1(old,[cycle;cycle(1,:)],query,'linear');
    cycle = cycle./vecnorm(cycle,2,2);
    initialSupport = max(1,round(initialSupport/ratio));
end
horizontal = hypot(cycle(:,1),cycle(:,2));
short = [-cycle(:,2)./max(horizontal,eps), ...
    cycle(:,1)./max(horizontal,eps),zeros(size(horizontal))];
pole = horizontal < eps;
short(pole,:) = repmat([1,0,0],nnz(pole),1);
long = cross(cycle,short,2);
scan = struct('B',cycle,'long',long,'short',short, ...
    'length',size(cycle,1),'initialSupport',initialSupport, ...
    'source',source,'path1',path1,'spiral',spiral,'path2',path2);
cachedKey = key;
cachedScan = scan;
end

function path = boundaryCircle(point,step)
radius = hypot(point(1),point(2));
azimuth = atan2(point(2),point(1));
angle = (azimuth:step/max(radius,eps):azimuth+2*pi).';
path = [radius*cos(angle),radius*sin(angle),repmat(point(3),numel(angle),1)];
end

function path = generateSpiral(width,step)
% Same endpoint convention as lidarpath.m; no plotting side effects.
k = width/(2*pi);
phiHigh = pi/2-asin(width/2);
phiLow = asin(width/2);
arc = integral(@(phi)sqrt(k^2+sin(phi).^2),phiLow,phiHigh)/k;
solution = ode45(@(~,phi)-k./sqrt(k^2+sin(phi).^2),[0,arc], ...
    phiHigh,odeset('RelTol',1e-9,'AbsTol',1e-11));
phi = deval(solution,0:step:arc).';
theta = (phiHigh-phi)/k;
path = flipud([sin(phi).*cos(theta),sin(phi).*sin(theta),cos(phi)]);
end

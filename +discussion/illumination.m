function factor = illumination(cfg,beam,ranges,c3,cLong,cShort)
%ILLUMINATION Finite 7x7 target illumination with the current MC prefilters.
arguments
    cfg (1,1) struct
    beam (1,1) struct
    ranges (:,1) double
    c3 (:,1) double
    cLong (:,1) double
    cShort (:,1) double
end
s = cfg.system;
half = atan((s.targetWidth_m/2)./ranges);
factor = zeros(size(ranges));
width = beam.wd_rad/2;
nodes = ((1:s.targetSampleGrid)-0.5)/s.targetSampleGrid-0.5;
[xx,yy] = meshgrid(nodes);
xx = xx(:).';
yy = yy(:).';
if beam.name == "ring"
    theta = acos(max(-1,min(1,c3)));
    candidate = abs(theta-beam.wD_rad/2) <= ...
        atan(sqrt(2)*(s.targetWidth_m/2)./ranges)+3*width;
    angularSide = 2*half(candidate);
    sampleX = theta(candidate)+angularSide.*xx;
    sampleY = angularSide.*yy;
    sampleRadius = hypot(sampleX,sampleY);
    intensity = exp(-2*((sampleRadius-beam.wD_rad/2)/width).^2);
    if cfg.ringGap_deg > 0
        azimuth = atan2(cShort(candidate),cLong(candidate));
        sampleAzimuth = azimuth+atan2(sampleY,sampleX);
        wrapped = atan2(sin(sampleAzimuth),cos(sampleAzimuth));
        keep = abs(wrapped) >= deg2rad(cfg.ringGap_deg)/2;
        intensity = intensity.*keep;
        if cfg.gapRedistribute
            intensity = intensity/(1-cfg.ringGap_deg/360);
        end
    end
else
    rotation = deg2rad(cfg.lineRotation_deg);
    long = atan2(cLong*cos(rotation)+cShort*sin(rotation),c3);
    short = atan2(-cLong*sin(rotation)+cShort*cos(rotation),c3);
    candidate = c3 > 0 & abs(long) <= beam.wD_rad/2+half & ...
        abs(short) <= 3*beam.wd_rad+half;
    angularSide = 2*half(candidate);
    sampleLong = long(candidate)+angularSide.*xx;
    sampleShort = short(candidate)+angularSide.*yy;
    intensity = (abs(sampleLong) <= beam.wD_rad/2).* ...
        exp(-2*(sampleShort/width).^2);
end
factor(candidate) = mean(intensity,2);
end

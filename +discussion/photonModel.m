function model = photonModel(cfg,beam)
%PHOTONMODEL Parameterized copy of MC_ring_snr / MC_line_snr photon laws.
arguments
    cfg (1,1) struct
    beam (1,1) struct
end
s = cfg.system;
h = 6.62607015e-34;
c = 299792458;
energy = h*c/s.wavelength_m;
area = pi*(s.receiverDiameter_m/2)^2;
efficiency = s.opticalTransmission*s.quantumEfficiency;
nTransmit = s.pulseEnergy_J/energy;
gaussianWidth = beam.wd_rad/2;
if beam.name == "line"
    omegaBeam = beam.wD_rad*gaussianWidth*sqrt(pi/2);
    omegaMatched = beam.wD_rad*6*beam.wd_rad;
    envelope = hypot(beam.wD_rad/2,3*beam.wd_rad);
else
    radius = beam.wD_rad/2;
    u0 = -sqrt(2)*radius/gaussianWidth;
    radial = radius*gaussianWidth*sqrt(pi)/(2*sqrt(2))*erfc(u0) + ...
        gaussianWidth^2/4*exp(-2*(radius/gaussianWidth)^2);
    omegaBeam = 2*pi*radial;
    omegaMatched = 2*pi*(cos(max(radius-3*beam.wd_rad,0)) - ...
        cos(min(radius+3*beam.wd_rad,pi)));
    envelope = radius+3*beam.wd_rad;
end
switch cfg.receiverMode
    case "matched"
        omegaFov = omegaMatched;
    case "envelope"
        omegaFov = 2*pi*(1-cos(envelope));
    case "common"
        assert(cfg.commonFov_rad >= envelope,'discussion:FovTooSmall', ...
            'Common receiver half-angle must cover both beam envelopes.');
        omegaFov = 2*pi*(1-cos(cfg.commonFov_rad));
    otherwise
        error('discussion:ReceiverMode','Unknown receiver mode.');
end
background = s.skyRadiance_W_m2_sr_nm*s.filterBandwidth_nm*area* ...
    omegaFov*s.rangeGate_s*efficiency/energy;
ranges = (0:s.noiseLutStep_m:s.maxRange_m).';
offsets = linspace(-c*s.rangeGate_s/4,c*s.rangeGate_s/4,101);
r = max(ranges+offsets,1e-6);
survival = 1;
if beam.name == "ring" && ~cfg.gapRedistribute
    survival = 1-cfg.ringGap_deg/360;
end
backscatter = survival*nTransmit*efficiency*area*trapz(offsets, ...
    s.backscatter_m_inv_sr_inv*exp(-2*s.extinction_m_inv*r)./r.^2,2);
model = struct('omegaBeam',omegaBeam,'omegaFov',omegaFov, ...
    'omegaMatched',omegaMatched,'envelope_rad',envelope, ...
    'background',background,'dark',s.darkCountRate_Hz*s.rangeGate_s, ...
    'backscatter',backscatter,'lutStep',s.noiseLutStep_m, ...
    'signalConstant',nTransmit*s.targetWidth_m^2*s.targetReflectivity* ...
    area*efficiency/pi/omegaBeam,'extinction',s.extinction_m_inv);
end

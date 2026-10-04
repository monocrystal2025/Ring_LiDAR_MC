function cases = cases(cfg)
%CASES Predeclared one-factor variations covering the original ten factors.
arguments
    cfg (1,1) struct
end
names = ["Pulse energy","PRF, fixed power","Scan rate", ...
    "Entry range","Target speed","Target side","Reflectivity", ...
    "Extinction","Sky radiance","Receiver diameter"];
fields = ["pulseEnergy_J","prf_Hz","scanRate_rad_s","maxRange_m", ...
    "targetSpeed_m_s","targetWidth_m","targetReflectivity","extinction_m_inv", ...
    "skyRadiance_W_m2_sr_nm","receiverDiameter_m"];
scales = [0.5,2;0.5,2;0.5,1.5;0.75,1.25;0.5,2; ...
    2/3,4/3;0.6,1.4;0.5,10;0.2,1e4;0.75,1.25];
cases = repmat(struct('name',"Nominal",'factor',0,'level',0, ...
    'scale',1,'config',cfg),21,1);
for j = 1:10
    for level = 1:2
        index = 1+2*(j-1)+level;
        cases(index).name = names(j);
        cases(index).factor = j;
        cases(index).level = level;
        cases(index).scale = scales(j,level);
        cases(index).config.system.(fields(j)) = ...
            cfg.system.(fields(j))*scales(j,level);
        if j == 2
            cases(index).config.system.pulseEnergy_J = ...
                cfg.system.pulseEnergy_J/scales(j,level);
        end
    end
end
end

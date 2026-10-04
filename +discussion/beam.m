function b = beam(name,width_mrad,smallWidth_rad)
%BEAM Use the existing full 1/e^2 width convention.
arguments
    name (1,1) string {mustBeMember(name,["line","ring"])}
    width_mrad (1,1) double {mustBePositive,mustBeFinite}
    smallWidth_rad (1,1) double {mustBePositive,mustBeFinite} = 1e-3
end
b = struct('name',name,'wD_rad',width_mrad*1e-3,'wd_rad',smallWidth_rad);
end

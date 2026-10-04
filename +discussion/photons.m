function [signal,noise] = photons(model,ranges,illumination,returnScale)
%PHOTONS Expected signal and noise counts per pulse, including shot noise.
arguments
    model (1,1) struct
    ranges (:,1) double
    illumination (:,1) double
    returnScale (1,1) double {mustBeNonnegative} = 1
end
signal = zeros(size(ranges));
valid = ranges > 0 & illumination > 0;
signal(valid) = returnScale*model.signalConstant*illumination(valid).* ...
    exp(-2*model.extinction*ranges(valid))./ranges(valid).^4;
index = max(1,min(numel(model.backscatter),round(ranges/model.lutStep)+1));
noise = model.backscatter(index)+model.background+model.dark;
end

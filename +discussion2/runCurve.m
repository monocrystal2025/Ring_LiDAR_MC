function output = runCurve(cfg,bank,kind)
%RUNCURVE Save each complete curve with validated physical inputs.
arguments
    cfg (1,1) struct
    bank (1,1) struct
    kind (1,1) string {mustBeMember(kind,["dynamic","static"])} = "dynamic"
end
names = ["+discussion2/geometry.m","+discussion2/positions.m", ...
    "+discussion2/scoreEvents.m","+discussion2/evaluate.m","+discussion2/curve.m", ...
    "+discussion/illumination.m","+discussion/photonModel.m","+discussion/photons.m","+discussion/scan.m"];
source = strings(size(names));
for j = 1:numel(names)
    source(j) = string(fileread(fullfile(cfg.root,names(j))));
end
inputs = struct('config',cfg,'bank',bank,'kind',kind,'source',source);
key = discussion2.hash(inputs);
folder = fullfile(cfg.cacheRoot,"curves");
if ~isfolder(folder)
    mkdir(folder);
end
filename = fullfile(folder,kind+"_"+key+".mat");
if isfile(filename)
    loaded = load(filename,'inputs','output');
    if isequaln(loaded.inputs,inputs)
        output = loaded.output;
        return
    end
end
output = discussion2.curve(cfg,bank,kind);
save(filename,'inputs','output','-v7.3');
end

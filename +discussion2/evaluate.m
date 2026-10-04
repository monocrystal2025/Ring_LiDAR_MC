function output = evaluate(cfg,beam,bank,kind)
%EVALUATE Apply the original photon law to cached geometric event sequences.
arguments
    cfg (1,1) struct
    beam (1,1) struct
    bank (1,1) struct
    kind (1,1) string {mustBeMember(kind,["dynamic","static"])}
end
record = discussion2.geometry(cfg,beam,bank,kind);
model = discussion.photonModel(cfg,beam);
raw = cell(bank.N,1);
windows = cell(bank.N,1);
for i = 1:bank.N
    [raw{i},windows{i}] = discussion2.scoreEvents( ...
        record.events{i},model,cfg,record.windowPulses);
end
output = struct('metrics',struct2table(vertcat(raw{:})), ...
    'windows',{windows},'N',bank.N,'width_mrad',beam.wD_rad*1e3, ...
    'beam',beam.name,'kind',kind,'geometryKey',record.key, ...
    'cycleTime_s',record.cycleTime_s,'horizon_s',record.horizon_s, ...
    'model',rmfield(model,'backscatter'));
end

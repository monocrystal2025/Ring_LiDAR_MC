function result = runCase(cfg,beam,bank,tag)
%RUNCASE Checkpoint each parameter/beam case with exact input validation.
arguments
    cfg (1,1) struct
    beam (1,1) struct
    bank (1,1) struct
    tag (1,1) string
end
folder = fullfile(cfg.outputRoot,"checkpoints");
if ~isfolder(folder)
    mkdir(folder);
end
filename = fullfile(folder,sprintf('%s_%s_%gmrad.mat', ...
    matlab.lang.makeValidName(tag),beam.name,beam.wD_rad*1e3));
inputs = struct('config',cfg,'beam',beam,'bank',bank, ...
    'implementationSHA256',discussion.fingerprint());
if isfile(filename)
    saved = load(filename,'inputs','result');
    if isfield(saved,'inputs') && isequaln(saved.inputs,inputs)
        result = saved.result;
        return
    end
end
result = discussion.simulate(cfg,beam,bank);
save(filename,'inputs','result','-v7');
end

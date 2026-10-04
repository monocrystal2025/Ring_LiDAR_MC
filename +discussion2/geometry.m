function record = geometry(cfg,beam,bank,kind)
%GEOMETRY Cache nonzero illumination events over the entire observation.
% Photometric parameter sweeps reuse exactly the same event sequence.
% Zero-signal pulses have q=0 and may be omitted from optimal accumulation;
% integer pulse indices preserve every gap and the fixed window length.
arguments
    cfg (1,1) struct
    beam (1,1) struct
    bank (1,1) struct
    kind (1,1) string {mustBeMember(kind,["dynamic","static"])}
end
sourceFiles = ["+discussion2/geometry.m","+discussion2/positions.m", ...
    "+discussion/illumination.m","+discussion/scan.m"];
source = strings(size(sourceFiles));
for j = 1:numel(sourceFiles)
    source(j) = string(fileread(fullfile(cfg.root,sourceFiles(j))));
end
pathFile = fullfile(cfg.scanLibrary,sprintf( ...
    'beam_vector_fast_5000Hz_h3000m_D%gu.npy',round(beam.wD_rad*1e6)));
pathInfo = dir(pathFile);
librarySignature = string(pathFile);
if ~isempty(pathInfo)
    librarySignature = librarySignature+"|"+pathInfo.bytes+"|"+pathInfo.datenum;
end
g = struct('version',cfg.version,'kind',kind,'beam',beam,'bank',bank, ...
    'range',cfg.system.maxRange_m,'speed',cfg.system.targetSpeed_m_s, ...
    'prf',cfg.system.prf_Hz,'omega',cfg.system.scanRate_rad_s, ...
    'targetSide',cfg.system.targetWidth_m,'targetGrid',cfg.system.targetSampleGrid, ...
    'phase',cfg.phaseMode,'horizon',cfg.maxTime_s, ...
    'source',source,'library',librarySignature,'scanSource',cfg.scanSource);
if kind == "static"
    g.speed = 0;
    g.phase = "fixed-start-one-complete-cycle";
    g.horizon = Inf;
end
key = discussion2.hash(g);
folder = fullfile(cfg.cacheRoot,"geometry");
if ~isfolder(folder)
    mkdir(folder);
end
filename = fullfile(folder,kind+"_"+beam.name+"_"+key+".mat");
if isfile(filename)
    saved = load(filename,'record','inputs');
    if isequaln(saved.inputs,g)
        record = saved.record;
        return
    end
end
scan = discussion.scan(cfg,beam.wD_rad);
[position,velocity,horizon] = discussion2.positions(cfg,bank,kind);
if kind == "static"
    horizon(:) = (scan.length-1)/cfg.system.prf_Hz;
end
events = cell(bank.N,1);
timerValue = tic;
for i = 1:bank.N
    if kind == "static"
        initial = 1;
    else
        support = scan.initialSupport;
        if cfg.phaseMode == "cycle"
            support = scan.length;
        end
        initial = floor(bank.u(i,5)*support)+1;
    end
    events{i} = oneTrack(cfg,beam,scan,position(i,:),velocity(i,:),horizon(i),initial);
end
record = struct('events',{events},'horizon_s',horizon,'position_m',position, ...
    'velocity_m_s',velocity,'N',bank.N,'kind',kind,'key',key, ...
    'width_mrad',beam.wD_rad*1e3,'windowPulses', ...
    max(1,ceil(cfg.system.prf_Hz*beam.wD_rad/cfg.system.scanRate_rad_s)), ...
    'cycleTime_s',scan.length/cfg.system.prf_Hz,'elapsed_s',toc(timerValue));
inputs = g;
save(filename,'record','inputs','-v7');
fprintf('Geometry %s %s %g mrad N=%d: %.1f s\n', ...
    kind,beam.name,beam.wD_rad*1e3,bank.N,record.elapsed_s);
end

function events = oneTrack(cfg,beam,scan,p,v,horizon,initial)
lastStep = floor(horizon*cfg.system.prf_Hz);
starts = 0:cfg.blockSize:lastStep;
chunks = cell(numel(starts),1);
for j = 1:numel(starts)
    step = (starts(j):min(lastStep,starts(j)+cfg.blockSize-1)).';
    t = step/cfg.system.prf_Hz;
    x = p(1)+t*v(1);
    y = p(2)+t*v(2);
    z = p(3)+t*v(3);
    ranges = max(sqrt(x.^2+y.^2+z.^2),1e-9);
    index = mod(initial-1+step,scan.length)+1;
    c3 = project(scan.B,index,x,y,z)./ranges;
    cLong = project(scan.long,index,x,y,z)./ranges;
    cShort = project(scan.short,index,x,y,z)./ranges;
    if numel(step)==1
        f = discussion.illumination(cfg,beam,[ranges;ranges], ...
            [c3;c3],[cLong;cLong],[cShort;cShort]);
        factor = f(1);
    else
        factor = discussion.illumination(cfg,beam,ranges,c3,cLong,cShort);
    end
    active = factor>0;
    chunks{j} = [step(active),ranges(active),factor(active)];
end
events = vertcat(chunks{:});
end

function dotProduct = project(basis,index,x,y,z)
dotProduct = x.*basis(index,1)+y.*basis(index,2)+z.*basis(index,3);
end

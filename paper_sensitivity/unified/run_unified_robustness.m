function result = run_unified_robustness(options)
%RUN_UNIFIED_ROBUSTNESS Simulate one design for three parallel benefits.
arguments
    options.N (1,1) double {mustBeInteger,mustBePositive} = 512
    options.StaticN (1,1) double {mustBeInteger,mustBePositive} = 2048
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 20260930
    options.Factors (1,:) double {mustBePositive} = [.8 .9 1 1.1 1.2]
    options.Widths (1,:) double {mustBePositive} = [5:5:100 125 150]
    options.ReferenceWidths (1,:) double {mustBePositive} = [100 125 150]
    options.Deadline (1,1) double {mustBePositive} = 15
    options.Workers (1,1) double {mustBeInteger,mustBeNonnegative} = 4
    options.MaxTasks (1,1) double {mustBePositive} = Inf
    options.OutputDirectory (1,1) string = ""
end
folder=fileparts(mfilename('fullpath'));root=fileparts(fileparts(folder));
addpath(root,folder,fileparts(folder));
if options.OutputDirectory=="",options.OutputDirectory=fullfile(folder,'results');end
if ~isfolder(options.OutputDirectory),mkdir(options.OutputDirectory);end
assert(issorted(options.Widths) && all(ismember(options.ReferenceWidths,options.Widths)), ...
    'Unified:Grid','Sorted grid must contain reference widths.');
assert(any(options.Factors==1) && all(options.Factors<=2),'Unified:Factors','Include baseline; rho cannot exceed 1.');
source=load(fullfile(root,'para_sensative_N3000_range2','para_manifest.mat'),'manifest');m=source.manifest;
old=rng;cleanup=onCleanup(@()rng(old));
rng(options.Seed,'twister');
scale=ones(1,3);factorIndex=0;
for p=1:3
    for factor=options.Factors(options.Factors~=1)
        row=ones(1,3);row(p)=factor;scale(end+1,:)=row;factorIndex(end+1,1)=p; %#ok<AGROW>
    end
end
values=scale.*m.baseline([1 2 4]);
[position,velocity]=init_UAV(m.R,1,30,options.N);phase=rand(options.N,1);
az=2*pi*rand(options.StaticN,1);co=rand(options.StaticN,1);
radius=m.R*nthroot(max(rand(options.StaticN,1),realmin),3);
stationary=radius.*[sqrt(1-co.^2).*cos(az),sqrt(1-co.^2).*sin(az),co];
names=["MC_line_snr.m","MC_ring_snr.m","append_and_check_photon_window.m","init_UAV.m"];
sources=cell(numel(names)+1,1);sources{1}=fileread(mfilename('fullpath')+".m");
for k=1:numel(names),sources{k+1}=fileread(fullfile(root,names(k)));end
signature=struct('options',rmfield(options,{'Workers','MaxTasks'}),'sources',{sources});
result=struct('options',options,'signature',signature,'scale',scale,'values',values, ...
    'factorIndex',factorIndex,'position',position,'velocity',velocity,'phase',phase, ...
    'staticPosition',stationary,'complete',false(numel(options.Widths),size(scale,1)), ...
    'manifest',m);
manifestFile=fullfile(options.OutputDirectory,'design.mat');
if isfile(manifestFile)
    previous=load(manifestFile,'result');
    assert(isequaln(previous.result.signature,signature),'Unified:Cache','Configuration/source changed; use a new directory.');
    result=previous.result;
else
    save(manifestFile,'result','-v7.3');
end
if options.Workers>0 && isempty(gcp('nocreate')),parpool('threads',options.Workers);end
widthOrder=[find(ismember(options.Widths,options.ReferenceWidths)),find(~ismember(options.Widths,options.ReferenceWidths))];
started=tic;taskCount=0;
for d=widthOrder
    width=options.Widths(d);scanIndex=find(m.widthsMrad==width,1);
    assert(~isempty(scanIndex),'Unified:Scan','Missing trajectory width.');
    scan=readNPY(char(m.scanFiles(scanIndex)));
    path1=boundaryCircle(scan(1,:),m.f,m.omega);path2=boundaryCircle(scan(end,:),m.f,m.omega);
    phaseLength=size(path1,1)+size(scan,1)+size(path2,1);cycleLength=phaseLength+size(scan,1);
    initial=floor(phase*phaseLength)+1;isReference=ismember(width,options.ReferenceWidths);
    for j=1:size(scale,1)
        file=fullfile(options.OutputDirectory,sprintf('case_%02d_width_%03d.mat',j,width));
        if isfile(file)
            cached=load(file,'meta');
            assert(cached.meta.Case==j && cached.meta.Width==width && cached.meta.N==options.N ...
                && cached.meta.Seed==options.Seed,'Unified:RawCache','Inconsistent raw cache.');
            result.complete(d,j)=true;continue
        end
        if taskCount>=options.MaxTasks,return;end
        physics=[values(j,1:2),values(j,3),values(j,3)/m.alphaBetaRatio];
        dynamicLimit=floor(options.Deadline*m.f);if isReference,dynamicLimit=Inf;end
        hit=false(options.N,2);time=nan(options.N,2);pulse=nan(options.N,2);cluster=pulse;
        staticHit=false(options.StaticN,2);windows=cell(1,2);
        for b=1:2
            if b==1,kernel=@MC_line_snr;beam='line';else,kernel=@MC_ring_snr;beam='ring';end
            block=mc_ea_block_size(width*1e-3,beam);
            [hit(:,b),time(:,b),pulse(:,b),cluster(:,b),windows{b}]=simulateBank( ...
                kernel,m,scan,path1,path2,position,velocity,initial,width,block,dynamicLimit,'fixed',physics,options.Workers);
            if isReference
                staticHit(:,b)=simulateBank(kernel,m,scan,path1,path2,stationary, ...
                    zeros(options.StaticN,3),ones(options.StaticN,1),width,block,cycleLength-1,'variable',physics,options.Workers);
            end
        end
        meta=struct('Case',j,'Width',width,'N',options.N,'StaticN',options.StaticN,'Seed',options.Seed, ...
            'Physics',physics,'DynamicLimit',dynamicLimit,'StaticAvailable',isReference,'StaticLimit',cycleLength-1);
        temporary=file+".partial.mat";
        save(temporary,'meta','hit','time','pulse','cluster','staticHit','windows','-v7.3');movefile(temporary,file,'f');
        result.complete(d,j)=true;taskCount=taskCount+1;
        save(manifestFile,'result','-v7.3');
        fprintf('UNIFIED case %d/%d width %g: timely=%s pulse=%s, elapsed %.1fs\n', ...
            j,size(scale,1),width,mat2str(mean(hit & time<=options.Deadline),3),mat2str(mean(pulse,1,'omitnan'),3),toc(started));
    end
end
save(manifestFile,'result','-v7.3');
fprintf('UNIFIED_SIMULATION_COMPLETE\n');
end

function [hit,time,pulse,cluster,windows]=simulateBank(kernel,m,scan,path1,path2,P,V,initial,width,block,limit,mode,physics,workers)
n=size(P,1);hit=false(n,1);time=nan(n,1);pulse=nan(n,1);cluster=pulse;
windows=repmat(empty_detection_window_photons(),n,1);
f=m.f;omega=m.omega;small=m.smallWidthRad;R=m.R;
parfor (i=1:n,workers)
    [hit(i),step,~,win]=kernel(width*1e-3,small,path1,scan,path2,f,omega,P(i,:),V(i,:),initial(i),R,block,limit,mode,physics);
    time(i)=step/f;windows(i)=win;
    if hit(i)
        positive=isfinite(win.signal_photons(:)) & win.signal_photons(:)>0;
        pulse(i)=nnz(positive);cluster(i)=nnz(positive & [true;~positive(1:end-1)]);
    end
end
end

function circle=boundaryCircle(point,f,omega)
theta0=atan2(point(2),point(1));r0=hypot(point(1),point(2));theta=(theta0:omega/f/r0:theta0+2*pi).';
circle=[r0*cos(theta),r0*sin(theta),ones(numel(theta),1)*point(3)];
end

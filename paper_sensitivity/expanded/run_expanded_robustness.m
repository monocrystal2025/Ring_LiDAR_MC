function result = run_expanded_robustness(options)
%RUN_EXPANDED_ROBUSTNESS Extend symmetric ranges, reusing identical old trials.
arguments
    options.CaseIDs (1,:) double {mustBeInteger,mustBePositive} = []
    options.Widths (1,:) double {mustBePositive} = 125
    options.Workers (1,1) double {mustBeInteger,mustBeNonnegative} = 4
    options.MaxTasks (1,1) double {mustBePositive} = Inf
end
folder=fileparts(mfilename('fullpath'));root=fileparts(fileparts(folder));
addpath(root,folder,fullfile(root,'paper_sensitivity'));
out=fullfile(folder,'results');if ~isfolder(out),mkdir(out);end
old=load(fullfile(root,'paper_sensitivity','unified','results','design.mat'),'result');old=old.result;
m=old.manifest;grid=5:5:400;
scales=ones(1,4);parameter=0;
levels={.5:.1:1.5,.5:.1:1.5,.1:.1:1.9,.5:.1:1.5};
for p=1:4
    for level=levels{p}
        if abs(level-1)<1e-10,continue;end
        row=ones(1,4);row(p)=round(level,1);scales(end+1,:)=row;parameter(end+1,1)=p; %#ok<AGROW>
    end
end
values=scales.*[m.baseline([1 2 4]),m.omega];
sourceNames=[fullfile(folder,"run_expanded_robustness.m"),fullfile(folder,"expanded_spiral.m"), ...
    fullfile(root,"MC_line_snr.m"),fullfile(root,"MC_ring_snr.m"), ...
    fullfile(root,"append_and_check_photon_window.m"),fullfile(root,"init_UAV.m"), ...
    fullfile(root,"empty_detection_window_photons.m"),fullfile(root,"mc_ea_block_size.m")];
sources=cell(size(sourceNames));for k=1:numel(sources),sources{k}=fileread(sourceNames(k));end
% Confirm the old-bank reuse uses precisely the same numerical kernels.
for k=1:4,assert(strcmp(sources{k+2},old.signature.sources{k+1}),'Expanded:OldSource','Old simulation kernel changed.');end
signature=struct('sources',{sources},'oldSignature',old.signature,'scale',scales,'grid',grid);
result=struct('signature',signature,'sourceNames',sourceNames,'scale',scales,'values',values, ...
    'factorIndex',parameter,'manifest',m,'bank',old,'widths',grid,'files',strings(numel(grid),size(scales,1)), ...
    'rawCase',zeros(numel(grid),size(scales,1)),'complete',false(numel(grid),size(scales,1)));
file=fullfile(out,'design.mat');
if isfile(file)
    saved=load(file,'result');assert(isequaln(signature,saved.result.signature),'Expanded:Cache','Source or bank changed.');result=saved.result;
else
    for j=1:size(scales,1)
        if scales(j,4)~=1,continue;end
        reuse=find(all(abs(old.scale-scales(j,1:3))<1e-10,2),1);
        if isempty(reuse),continue;end
        for width=old.options.Widths
            d=find(grid==width,1);raw=fullfile(root,'paper_sensitivity','unified','results',sprintf('case_%02d_width_%03d.mat',reuse,width));
            assert(isfile(raw),'Expanded:OldMissing','Missing old raw file.');
            result.files(d,j)=raw;result.rawCase(d,j)=reuse;result.complete(d,j)=true;
        end
    end
    save(file,'result','-v7.3');
end
ids=options.CaseIDs;if isempty(ids),ids=1:size(scales,1);end
assert(all(ids<=size(scales,1)) && all(ismember(options.Widths,grid)),'Expanded:Tasks','Unsupported cases or widths.');
if options.Workers>0 && isempty(gcp('nocreate')),parpool('threads',options.Workers);end
P=old.position;V=old.velocity;stationary=old.staticPosition;phase=old.phase;
started=tic;done=0;
for width=options.Widths
    d=find(grid==width,1);base=readNPY(char(m.scanFiles(m.widthsMrad==width)));
    for j=ids
        if result.complete(d,j),continue;end
        if done>=options.MaxTasks,return;end
        omega=values(j,4);local=m;local.omega=omega;
        if scales(j,4)==1,scan=base;
        else,scan=expanded_spiral(m.f,2*1000*tan(width*1e-3/2),omega,1000);end
        path1=boundaryCircle(scan(1,:),m.f,omega);path2=boundaryCircle(scan(end,:),m.f,omega);
        phaseLength=size(path1,1)+size(scan,1)+size(path2,1);cycleLength=phaseLength+size(scan,1);
        initial=floor(phase*phaseLength)+1;
        physics=[values(j,1:2),values(j,3),values(j,3)/m.alphaBetaRatio];
        reference=ismember(width,[100 125 150]);limit=floor(old.options.Deadline*m.f);if reference,limit=Inf;end
        hit=false(old.options.N,2);time=nan(old.options.N,2);pulse=time;cluster=time;
        staticHit=false(old.options.StaticN,2);windows=cell(1,2);
        for b=1:2
            if b==1,kernel=@MC_line_snr;beam='line';else,kernel=@MC_ring_snr;beam='ring';end
            block=mc_ea_block_size(width*1e-3,beam);
            [hit(:,b),time(:,b),pulse(:,b),cluster(:,b),windows{b}]=simulateBank(kernel,local,scan,path1,path2,P,V,initial,width,block,limit,'fixed',physics,options.Workers);
            if reference
                staticHit(:,b)=simulateBank(kernel,local,scan,path1,path2,stationary,zeros(size(stationary)), ...
                    ones(size(stationary,1),1),width,block,cycleLength-1,'variable',physics,options.Workers);
            end
        end
        meta=struct('Case',j,'Width',width,'N',old.options.N,'StaticN',old.options.StaticN,'Seed',old.options.Seed, ...
            'Physics',physics,'Omega',omega,'DynamicLimit',limit,'StaticAvailable',reference, ...
            'StaticLimit',cycleLength-1,'PhaseLength',phaseLength,'CycleLength',cycleLength);
        raw=fullfile(out,sprintf('case_%02d_width_%03d.mat',j,width));temporary=raw+".partial.mat";
        save(temporary,'meta','hit','time','pulse','cluster','staticHit','windows','-v7.3');movefile(temporary,raw,'f');
        result.files(d,j)=raw;result.rawCase(d,j)=j;result.complete(d,j)=true;
        save(file,'result','-v7.3');done=done+1;
        fprintf('EXPANDED case %d scale=%s width=%g PD=%s pulse=%s blind=%s elapsed=%.1fs\n', ...
            j,mat2str(scales(j,:),2),width,mat2str(mean(hit & time<=15),3),mat2str(mean(pulse,1,'omitnan'),3),mat2str(mean(~staticHit),3),toc(started));
    end
end
fprintf('EXPANDED_REQUEST_COMPLETE\n');
end

function [hit,time,pulse,cluster,windows]=simulateBank(kernel,m,scan,path1,path2,P,V,initial,width,block,limit,mode,physics,workers)
n=size(P,1);hit=false(n,1);time=nan(n,1);pulse=time;cluster=time;
windows=repmat(empty_detection_window_photons(),n,1);f=m.f;omega=m.omega;small=m.smallWidthRad;R=m.R;
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

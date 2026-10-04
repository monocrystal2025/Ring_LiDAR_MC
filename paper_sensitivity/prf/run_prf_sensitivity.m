function result=run_prf_sensitivity(options)
%RUN_PRF_SENSITIVITY Paired PRF sweep at fixed pulse energy and scan speed.
arguments
    options.Frequencies (1,:) double {mustBePositive,mustBeInteger}=1000:1000:10000
    options.Widths (1,:) double {mustBePositive}=125
    options.Workers (1,1) double {mustBeInteger,mustBeNonnegative}=4
end
folder=fileparts(mfilename('fullpath'));root=fileparts(fileparts(folder));
addpath(root,fileparts(folder),fullfile(fileparts(folder),'expanded'));
out=fullfile(folder,'results');if ~isfolder(out),mkdir(out);end
s=load(fullfile(fileparts(folder),'expanded','results','design.mat'),'result');old=s.result;
bank=old.bank;m=old.manifest;frequencies=1000:1000:10000;widths=5:5:400;
names=[string(mfilename('fullpath'))+".m",old.sourceNames(2:end)];
sources=cell(size(names));for k=1:numel(names),sources{k}=fileread(names(k));end
for k=2:numel(names),assert(strcmp(sources{k},old.signature.sources{k}),'PRF:Source','Existing numerical source changed.');end
signature=struct('sources',{sources},'bankSignature',bank.signature,'frequencies',frequencies,'widths',widths,'energyMode','fixed pulse energy');
result=struct('signature',signature,'sourceNames',names,'frequencies',frequencies,'widths',widths,'bank',bank,'manifest',m, ...
    'files',strings(numel(widths),numel(frequencies)),'complete',false(numel(widths),numel(frequencies)));
designFile=fullfile(out,'design.mat');
if isfile(designFile)
    s=load(designFile,'result');assert(isequaln(signature,s.result.signature),'PRF:Cache','Source or bank mismatch.');result=s.result;
else
    base=find(frequencies==m.f);result.files(:,base)=old.files(:,1);result.complete(:,base)=old.complete(:,1);
    save(designFile,'result','-v7.3');
end
assert(all(ismember(options.Frequencies,frequencies)) && all(ismember(options.Widths,widths)));
if options.Workers>0 && isempty(gcp('nocreate')),parpool('threads',options.Workers);end
P=bank.position;V=bank.velocity;stationary=bank.staticPosition;
physics=[m.baseline(1:2),m.baseline(4),m.baseline(4)/m.alphaBetaRatio];started=tic;
for width=options.Widths
    d=find(widths==width);
    for f=options.Frequencies
        j=find(frequencies==f);if result.complete(d,j),continue;end
        scan=expanded_spiral(f,2000*tan(width*1e-3/2),m.omega,1000);
        path1=boundaryCircle(scan(1,:),f,m.omega);path2=boundaryCircle(scan(end,:),f,m.omega);
        phaseLength=size(path1,1)+size(scan,1)+size(path2,1);cycleLength=phaseLength+size(scan,1);
        initial=floor(bank.phase*phaseLength)+1;local=m;local.f=f;
        reference=width==125;limit=floor(bank.options.Deadline*f);if reference,limit=Inf;end
        hit=false(bank.options.N,2);time=nan(bank.options.N,2);pulse=time;
        staticHit=false(bank.options.StaticN,2);windows=cell(1,2);
        for b=1:2
            if b==1,kernel=@MC_line_snr;beam='line';else,kernel=@MC_ring_snr;beam='ring';end
            block=mc_ea_block_size(width*1e-3,beam);
            [hit(:,b),time(:,b),pulse(:,b),windows{b}]=simulateBank(kernel,local,scan,path1,path2,P,V,initial,width,block,limit,'fixed',physics,options.Workers);
            if reference
                staticHit(:,b)=simulateBank(kernel,local,scan,path1,path2,stationary,zeros(size(stationary)), ...
                    ones(size(stationary,1),1),width,block,cycleLength-1,'variable',physics,options.Workers);
            end
        end
        adjacent=sum(scan(1:end-1,:).*scan(2:end,:),2);angularStep=acos(max(-1,min(1,adjacent)));
        meta=struct('Frequency',f,'Width',width,'N',bank.options.N,'StaticN',bank.options.StaticN,'Seed',bank.options.Seed, ...
            'Physics',physics,'Omega',m.omega,'DynamicLimit',limit,'StaticAvailable',reference,'StaticLimit',cycleLength-1, ...
            'PhaseLength',phaseLength,'CycleLength',cycleLength,'CycleSeconds',cycleLength/f, ...
            'MaxRelativeScanSpeedError',max(abs(angularStep*f/m.omega-1)));
        raw=fullfile(out,sprintf('prf_%05d_width_%03d.mat',f,width));temporary=raw+".partial.mat";
        save(temporary,'meta','hit','time','pulse','staticHit','windows','-v7.3');movefile(temporary,raw,'f');
        result.files(d,j)=raw;result.complete(d,j)=true;save(designFile,'result','-v7.3');
        fprintf('PRF f=%g width=%g PD=%s pulse=%s',f,width,mat2str(mean(hit & time<=15),4),mat2str(mean(pulse,1,'omitnan'),4));
        if reference,fprintf(' blind=%s',mat2str(mean(~staticHit),4));end
        fprintf(' elapsed=%.1fs\n',toc(started));
    end
end
end

function [hit,time,pulse,windows]=simulateBank(kernel,m,scan,path1,path2,P,V,initial,width,block,limit,mode,physics,workers)
n=size(P,1);hit=false(n,1);time=nan(n,1);pulse=time;
windows=repmat(empty_detection_window_photons(),n,1);f=m.f;omega=m.omega;small=m.smallWidthRad;R=m.R;
parfor (i=1:n,workers)
    [hit(i),step,~,win]=kernel(width*1e-3,small,path1,scan,path2,f,omega,P(i,:),V(i,:),initial(i),R,block,limit,mode,physics);
    time(i)=step/f;windows(i)=win;
    if hit(i),pulse(i)=nnz(isfinite(win.signal_photons(:)) & win.signal_photons(:)>0);end
end
end

function circle=boundaryCircle(point,f,omega)
theta0=atan2(point(2),point(1));r0=hypot(point(1),point(2));theta=(theta0:omega/f/r0:theta0+2*pi).';
circle=[r0*cos(theta),r0*sin(theta),ones(numel(theta),1)*point(3)];
end

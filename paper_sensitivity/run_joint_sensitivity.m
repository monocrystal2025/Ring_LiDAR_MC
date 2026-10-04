function result = run_joint_sensitivity(options)
%RUN_JOINT_SENSITIVITY Paired joint perturbations using the manuscript kernels.
% All four physical factors vary simultaneously using a seeded Latin
% hypercube. The baseline is retained as an additional, separate scenario.
% This is an operating-box robustness experiment, not a Sobol analysis.
arguments
    options.N (1,1) double {mustBeInteger,mustBePositive} = 256
    options.ScenarioCount (1,1) double {mustBeInteger,mustBePositive} = 24
    options.RelativeHalfRange (1,1) double {mustBePositive} = .2
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 20260927
    options.WidthMrad (1,1) double {mustBePositive} = 100
    options.Deadline (1,1) double {mustBePositive} = 15
    options.UseParallel (1,1) logical = true
    options.Workers (1,1) double {mustBeInteger,mustBePositive} = 4
    options.OutputDirectory (1,1) string = ""
    options.DataDirectory (1,1) string = ""
    options.MaxScenarios (1,1) double {mustBePositive} = Inf
end
root=fileparts(fileparts(mfilename('fullpath')));
if options.DataDirectory=="",options.DataDirectory=fullfile(root,'para_sensative_N3000_range2');end
if options.OutputDirectory=="",options.OutputDirectory=fullfile(root,'paper_sensitivity','results',sprintf('joint_N%d',options.N));end
assert(options.RelativeHalfRange<1,'Sensitivity:Range','Half-range must be below 1.');
if ~isfolder(options.OutputDirectory),mkdir(options.OutputDirectory);end
saved=load(fullfile(options.DataDirectory,'para_manifest.mat'),'manifest');
m=saved.manifest;
assert(m.complete && ~m.truncated,'Sensitivity:Source','Complete original data required.');
d=find(m.widthsMrad==options.WidthMrad,1);
assert(~isempty(d),'Sensitivity:Width','Width must match a saved scan trajectory.');
oldRng=rng;cleanup=onCleanup(@()rng(oldRng)); %#ok<NASGU>
rng(options.Seed,'twister');
u=zeros(options.ScenarioCount,4);
for p=1:4,u(:,p)=(randperm(options.ScenarioCount).'-rand(options.ScenarioCount,1))/options.ScenarioCount;end
scale=[ones(1,4);1+options.RelativeHalfRange*(2*u-1)];
values=scale.*m.baseline;
[position,velocity]=init_UAV(m.R,1,30,options.N);
phase=rand(options.N,1);
az=2*pi*rand(options.N,1);co=rand(options.N,1);
radius=m.R*nthroot(max(rand(options.N,1),realmin),3);
stationary=radius.*[sqrt(1-co.^2).*cos(az),sqrt(1-co.^2).*sin(az),co];
scan=readNPY(char(m.scanFiles(d)));
path1=boundaryCircle(scan(1,:),m.f,m.omega);
path2=boundaryCircle(scan(end,:),m.f,m.omega);
phaseLength=size(path1,1)+size(scan,1)+size(path2,1);
cycleLength=phaseLength+size(scan,1);
indices=floor(phase*phaseLength)+1;
sourceNames=["run_joint_sensitivity.m","MC_line_snr.m","MC_ring_snr.m", ...
    "append_and_check_photon_window.m","init_UAV.m"];
sources=cell(size(sourceNames));
sources{1}=fileread(mfilename('fullpath')+".m");
for j=2:numel(sourceNames),sources{j}=fileread(fullfile(root,sourceNames(j)));end
signature=struct('options',rmfield(options,{'UseParallel','Workers','MaxScenarios'}), ...
    'sources',{sources},'scanFile',m.scanFiles(d),'scanSize',size(scan));
result=struct('options',options,'signature',signature,'values',values,'scale',scale, ...
    'position',position,'velocity30',velocity,'staticPosition',stationary, ...
    'phase',phase,'phaseLength',phaseLength,'cycleLength',cycleLength, ...
    'movingHit',false(options.N,2,size(values,1)), ...
    'movingTime',nan(options.N,2,size(values,1)), ...
    'staticHit',false(options.N,2,size(values,1)), ...
    'complete',false(size(values,1),1));
file=fullfile(options.OutputDirectory,'joint_results.mat');
if isfile(file)
    prev=load(file,'result');
    assert(isequaln(prev.result.signature,signature),'Sensitivity:Cache','Existing run differs; choose another OutputDirectory.');
    result=prev.result;
end
% Save the complete design before any responses are simulated.
save(file,'result','-v7.3');
workers=0;
if options.UseParallel
    if isempty(gcp('nocreate')),parpool('threads',options.Workers);end
    workers=options.Workers;
end
width=options.WidthMrad*1e-3;f=m.f;omega=m.omega;R=m.R;small=m.smallWidthRad;
movingLimit=floor(options.Deadline*f);staticLimit=cycleLength-1;
started=tic;
pending=find(~result.complete);
pending=pending(1:min(numel(pending),options.MaxScenarios));
for j=pending.'
    physics=[values(j,1:2),values(j,4),values(j,4)/m.alphaBetaRatio];
    movingV=velocity*(values(j,3)/30);
    for b=1:2
        if b==1,kernel=@MC_line_snr;beam='line';else,kernel=@MC_ring_snr;beam='ring';end
        for mode=1:2
            if mode==1
                P=position;V=movingV;initial=indices;limit=movingLimit;window='fixed';
            else
                P=stationary;V=zeros(options.N,3);initial=ones(options.N,1);limit=staticLimit;window='variable';
            end
            block=min(mc_ea_block_size(width,beam),limit);
            hit=false(options.N,1);time=nan(options.N,1);
            parfor (i=1:options.N,workers)
                [hit(i),time(i)]=kernel(width,small,path1,scan,path2,f,omega, ...
                    P(i,:),V(i,:),initial(i),R,block,limit,window,physics);
            end
            if mode==1,result.movingHit(:,b,j)=hit;result.movingTime(:,b,j)=time/f;
            else,result.staticHit(:,b,j)=hit;end
        end
    end
    result.complete(j)=true;
    temp=file+".partial.mat";save(temp,'result','-v7.3');movefile(temp,file,'f');
    fprintf('Joint %d/%d: E=%.1f rho=%.3f v=%.2f alpha=%.3g; PD=%s blind=%s; %.1f s\n', ...
        j,size(values,1),values(j,1)*1e6,values(j,2),values(j,3),values(j,4), ...
        mat2str(mean(result.movingHit(:,:,j)),3),mat2str(mean(~result.staticHit(:,:,j)),3),toc(started));
end
end

function circle=boundaryCircle(point,f,omega)
theta0=atan2(point(2),point(1));r0=hypot(point(1),point(2));
theta=(theta0:omega/f/r0:theta0+2*pi).';
circle=[r0*cos(theta),r0*sin(theta),ones(numel(theta),1)*point(3)];
end

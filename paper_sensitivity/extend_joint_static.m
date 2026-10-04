function result = extend_joint_static(options)
%EXTEND_JOINT_STATIC Uniformly extend all static scenarios to a fixed N.
% The original 512 static locations are retained; additional locations use
% one independent seeded sample bank shared across every beam and scenario.
% Dynamic trials and all parameter combinations remain unchanged.
arguments
    options.InputFile (1,1) string = ""
    options.OutputDirectory (1,1) string = ""
    options.StaticN (1,1) double {mustBeInteger,mustBePositive} = 2048
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 20260929
    options.Workers (1,1) double {mustBeInteger,mustBePositive} = 4
    options.UseParallel (1,1) logical = true
    options.MaxScenarios (1,1) double {mustBePositive} = Inf
end
folder=fileparts(mfilename('fullpath'));
if options.InputFile==""
    options.InputFile=fullfile(folder,'results','joint_N512','joint_results.mat');
end
if options.OutputDirectory==""
    options.OutputDirectory=fullfile(folder,'results',sprintf('joint_static_N%d',options.StaticN));
end
if ~isfolder(options.OutputDirectory)
    mkdir(options.OutputDirectory);
end
original=load(options.InputFile,'result');
original=original.result;
assert(all(original.complete),'Sensitivity:Incomplete','Finish the original joint experiment first.');
baseN=size(original.staticHit,1);
assert(options.StaticN>baseN,'Sensitivity:StaticN','StaticN must exceed the original static sample size.');
saved=load(fullfile(original.options.DataDirectory,'para_manifest.mat'),'manifest');
m=saved.manifest;
extension=struct('StaticN',options.StaticN,'Seed',options.Seed, ...
    'OriginalSignature',original.signature,'Source',fileread(mfilename('fullpath')+".m"));
file=fullfile(options.OutputDirectory,'joint_results.mat');
if isfile(file)
    saved=load(file,'result');
    result=saved.result;
    assert(isequaln(result.extensionSignature,extension),'Sensitivity:Cache','Extension differs; use a new output directory.');
else
    oldRng=rng;
    cleanup=onCleanup(@()rng(oldRng)); %#ok<NASGU>
    rng(options.Seed,'twister');
    extraN=options.StaticN-baseN;
    az=2*pi*rand(extraN,1);
    co=rand(extraN,1);
    radius=m.R*nthroot(max(rand(extraN,1),realmin),3);
    extra=radius.*[sqrt(1-co.^2).*cos(az),sqrt(1-co.^2).*sin(az),co];
    result=original;
    result.staticPosition=[original.staticPosition;extra];
    result.staticHit=false(options.StaticN,2,size(original.values,1));
    result.staticHit(1:baseN,:,:)=original.staticHit;
    result.options.StaticN=options.StaticN;
    result.extensionSignature=extension;
    result.staticComplete=false(size(result.complete));
    save(file,'result','-v7.3');
end
d=find(m.widthsMrad==original.options.WidthMrad,1);
scan=readNPY(char(m.scanFiles(d)));
path1=boundaryCircle(scan(1,:),m.f,m.omega);
path2=boundaryCircle(scan(end,:),m.f,m.omega);
limit=result.cycleLength-1;
workers=0;
if options.UseParallel
    if isempty(gcp('nocreate'))
        parpool('threads',options.Workers);
    end
    workers=options.Workers;
end
P=result.staticPosition(baseN+1:end,:);
extraN=size(P,1);
width=original.options.WidthMrad*1e-3;
f=m.f;omega=m.omega;R=m.R;small=m.smallWidthRad;
pending=find(~result.staticComplete);
pending=pending(1:min(numel(pending),options.MaxScenarios));
started=tic;
for j=pending.'
    value=result.values(j,:);
    physics=[value(1:2),value(4),value(4)/m.alphaBetaRatio];
    for b=1:2
        if b==1
            kernel=@MC_line_snr;beam='line';
        else
            kernel=@MC_ring_snr;beam='ring';
        end
        block=min(mc_ea_block_size(width,beam),limit);
        hit=false(extraN,1);
        parfor (i=1:extraN,workers)
            hit(i)=kernel(width,small,path1,scan,path2,f,omega,P(i,:), ...
                [0 0 0],1,R,block,limit,'variable',physics);
        end
        result.staticHit(baseN+1:end,b,j)=hit;
    end
    result.staticComplete(j)=true;
    temporary=file+".partial.mat";
    save(temporary,'result','-v7.3');
    movefile(temporary,file,'f');
    fprintf('Static extension %d/%d, N=%d, blind=%s, %.1f s\n', ...
        j,numel(result.staticComplete),options.StaticN,mat2str(mean(~result.staticHit(:,:,j)),4),toc(started));
end
end

function circle=boundaryCircle(point,f,omega)
theta0=atan2(point(2),point(1));r0=hypot(point(1),point(2));
theta=(theta0:omega/f/r0:theta0+2*pi).';
circle=[r0*cos(theta),r0*sin(theta),ones(numel(theta),1)*point(3)];
end

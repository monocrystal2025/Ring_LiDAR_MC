function result = receiver_noise_data(varargin)
%RECEIVER_NOISE_DATA Independent receiver-noise robustness study.
%
%   eta_R = sqrt(B0/B) * sqrt(tGate0/tGate) * sqrt(DCR0/DCR).
%   A single-mechanism curve changes one factor. The coupled curve changes
%   all three factors, with each contributing eta_R^(1/3).
%
%   This file independently constructs cases, runs the Monte Carlo model,
%   evaluates confidence intervals and saves its MAT output. It does not
%   call the other normalized-study data files. wd is fixed at 1.0 mrad.

scriptDir=fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir=pwd; end
moduleDir=fullfile(scriptDir,'ea_sensitivity_analysis');
defaultOutput=fullfile(scriptDir,'LUBANG_data','receiver_noise_data.mat');
p=inputParser; p.FunctionName=mfilename;
addParameter(p,'Mode','preview',@isTextScalar);
addParameter(p,'N',[],@optionalInteger);
addParameter(p,'WDGrid_mrad',[],@optionalVector);
addParameter(p,'TimeLimit_s',15,@nonnegativeScalar);
addParameter(p,'TargetProbability',0.80,@probabilityScalar);
addParameter(p,'UseParallel',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'Seed',20260817,@nonnegativeInteger);
addParameter(p,'OutputFile',defaultOutput,@isTextScalar);
addParameter(p,'SaveData',true,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
cfg=makeRunConfig(p.Results,scriptDir,moduleDir);
verifySimulator(moduleDir);
addedPath=addSimulatorPath(moduleDir);
pathCleanup=onCleanup(@()restorePath(moduleDir,addedPath));

definitions=buildReceiverCases();
baseConfig=makeSimulationConfig(cfg);
bank=ea_sens_build_scenario_bank(cfg.N,cfg.Seed);
cases=repmat(emptyCase(),numel(definitions),1);
fprintf(['Generating independent receiver-noise study: N=%d, ', ...
    'wd=1.0 mrad, %d wD values, %d cases.\n'], ...
    cfg.N,numel(cfg.WDGrid_mrad),numel(definitions));
t=tic;
for k=1:numel(definitions)
    fprintf('[%d/%d] %s, eta_R=%.4g\n',k,numel(definitions), ...
        definitions(k).caseLabel,definitions(k).xValue);
    simConfig=applyOverrides(baseConfig,definitions(k).overrides);
    line=simulateCurve(simConfig,bank,'line',cfg);
    annular=simulateCurve(simConfig,bank,'ring',cfg);
    cases(k)=packCase(definitions(k),line,annular, ...
        equalPerformance(line,annular,cfg.TargetProbability));
end
plotData=makePlotData(cases);
summary=makeSummary(plotData); disp(summary);
result=packResult(cfg,definitions,cases,plotData,summary,toc(t));
if cfg.SaveData
    folder=fileparts(cfg.OutputFile);
    if ~isfolder(folder), mkdir(folder); end
    save(cfg.OutputFile,'result','-v7.3');
    fprintf('Saved MAT data: %s\n',cfg.OutputFile);
end
fprintf('Receiver-noise study finished in %.1f s.\n',result.elapsed_s);
end

function definitions=buildReceiverCases()
baseline=baselineSystem();
% Two-decade logarithmic-style grid shared by all four mechanisms.
qualityGrid=[0.10,0.20,0.50,1.00,2.00,5.00,10.0];
definitions=repmat(emptyDefinition(),4*numel(qualityGrid),1); k=0;
for quality=qualityGrid
    bandwidth=baseline.filterBandwidth_nm/quality^2;
    overrides=struct('filterBandwidth_nm',bandwidth); k=k+1;
    definitions(k)=receiverDefinition(k,sprintf('Bandwidth %.3g nm',bandwidth), ...
        'Filter bandwidth',quality,overrides);
end
for quality=qualityGrid
    gate=baseline.rangeGate_s/quality^2;
    overrides=struct('rangeGate_s',gate); k=k+1;
    definitions(k)=receiverDefinition(k,sprintf('Gate %.3g ns',gate*1e9), ...
        'Time gate',quality,overrides);
end
for quality=qualityGrid
    darkRate=baseline.darkCountRate_Hz/quality^2;
    overrides=struct('darkCountRate_Hz',darkRate); k=k+1;
    definitions(k)=receiverDefinition(k,sprintf('DCR %.3g cps',darkRate), ...
        'Dark count',quality,overrides);
end
for quality=qualityGrid
    componentQuality=quality^(1/3);
    overrides=struct('filterBandwidth_nm', ...
        baseline.filterBandwidth_nm/componentQuality^2, ...
        'rangeGate_s',baseline.rangeGate_s/componentQuality^2, ...
        'darkCountRate_Hz',baseline.darkCountRate_Hz/componentQuality^2);
    k=k+1; definitions(k)=receiverDefinition(k, ...
        sprintf('Receiver quality %.3g',quality),'Coupled receiver noise', ...
        quality,overrides);
end
end

function d=receiverDefinition(index,label,series,quality,overrides)
d=struct('caseId',"receiver_"+index,'caseLabel',string(label), ...
    'series',string(series),'xValue',double(quality), ...
    'overrides',overrides);
end
function s=baselineSystem()
s=struct('wavelength_m',1550e-9,'receiverDiameter_m',0.0508, ...
    'opticalTransmission',0.80,'quantumEfficiency',0.80, ...
    'rangeGate_s',50e-9,'filterBandwidth_nm',10, ...
    'darkCountRate_Hz',400,'skyRadiance_W_m2_sr_nm',5e-9);
end

function cfg=makeRunConfig(o,scriptDir,moduleDir)
switch lower(string(o.Mode))
    case "smoke", n0=4; wd0=[80,140,200];
    case "preview", n0=200; wd0=[50,75,100,125,150,190,230];
    case "paper", n0=3000; wd0=40:10:260;
    otherwise, error('receiver_noise_data:UnknownMode', ...
            'Mode must be smoke, preview or paper.');
end
if isempty(o.N), n=n0; else, n=o.N; end
if isempty(o.WDGrid_mrad), wd=wd0; else, wd=o.WDGrid_mrad; end
cfg=struct('mode',lower(string(o.Mode)),'N',double(n), ...
    'WDGrid_mrad',unique(double(wd(:).'),'sorted'), ...
    'TimeLimit_s',double(o.TimeLimit_s), ...
    'TargetProbability',double(o.TargetProbability), ...
    'UseParallel',o.UseParallel,'Seed',double(o.Seed), ...
    'OutputFile',char(o.OutputFile),'SaveData',o.SaveData, ...
    'smallWidth_mrad',1.0,'scriptDir',scriptDir,'moduleDir',moduleDir, ...
    'wilsonZ',1.95996398454005);
end
function config=makeSimulationConfig(cfg)
config=ea_sens_default_config();
config.seed=cfg.Seed; config.execution.useParallel=cfg.UseParallel;
config.execution.showProgress=false; config.execution.blockSize=25000;
s=config.system; s.maxRange_m=2000; s.prf_Hz=5000;
s.scanRate_rad_s=2*pi; s.targetSpeed_m_s=30; s.beamWidth_rad=1e-3;
s.pulseEnergy_J=150e-6; s.wavelength_m=1550e-9;
s.receiverDiameter_m=0.0508; s.opticalTransmission=0.80;
s.quantumEfficiency=0.80; s.rangeGate_s=50e-9;
s.filterBandwidth_nm=10; s.darkCountRate_Hz=400;
s.extinction_m_inv=1.5e-5; s.backscatter_m_inv_sr_inv=0.3e-6;
s.skyRadiance_W_m2_sr_nm=5e-9; s.targetWidth_m=0.30;
s.targetReflectivity=0.50; s.snrThreshold=2;
s.gaussianCutoffWidths=3; s.targetSampleGrid=7; s.noiseLutStep_m=0.5;
config.system=s;
end
function config=applyOverrides(config,overrides)
names=fieldnames(overrides);
for k=1:numel(names), config.system.(names{k})=overrides.(names{k}); end
config.system.beamWidth_rad=1e-3;
end

function curve=simulateCurve(config,bank,beamName,cfg)
nWD=numel(cfg.WDGrid_mrad); counts=zeros(1,nWD);
for j=1:nWD
    beam=struct('name',string(beamName),'label',string(beamName), ...
        'wD_rad',cfg.WDGrid_mrad(j)*1e-3,'wd_rad',1e-3);
    raw=ea_sens_simulate_case(config,beam,bank);
    firstTime=raw.firstTime_s(:);
    timely=logical(raw.detect(:))&isfinite(firstTime)&firstTime>=0& ...
        firstTime<=cfg.TimeLimit_s;
    counts(j)=nnz(timely);
end
prob=counts./bank.N; [low,high]=wilson(counts,bank.N,cfg.wilsonZ);
curve=struct('beam',string(beamName),'wD_mrad',cfg.WDGrid_mrad, ...
    'wd_mrad',1.0,'N',bank.N,'timelyCounts',counts, ...
    'probability',prob,'probabilityLow95',low,'probabilityHigh95',high);
end
function equal=equalPerformance(line,ring,target)
lm=crossing(line.wD_mrad,line.probability,target);
rm=crossing(ring.wD_mrad,ring.probability,target);
ll=crossing(line.wD_mrad,line.probabilityHigh95,target);
lh=crossing(line.wD_mrad,line.probabilityLow95,target);
rl=crossing(ring.wD_mrad,ring.probabilityHigh95,target);
rh=crossing(ring.wD_mrad,ring.probabilityLow95,target);
equal=struct('lineMinWD_mrad',lm,'ringMinWD_mrad',rm, ...
    'kappa',safeRatio(rm,lm),'kappaLow95',safeRatio(rl,lh), ...
    'kappaHigh95',safeRatio(rh,ll));
equal.robustAdvantage=isfinite(equal.kappaHigh95)&&equal.kappaHigh95<1;
equal.advantageDisappears=isfinite(equal.kappa)&&equal.kappa>=1;
end
function value=crossing(wD,p,target)
wD=wD(:); p=p(:); [wD,order]=sort(wD); p=p(order);
valid=isfinite(wD)&isfinite(p);
wD=wD(valid); p=cummax(p(valid)); index=find(p>=target,1);
if isempty(index), value=nan; return; end
if index==1, value=wD(1); return; end
p1=p(index-1); p2=p(index);
if p2<=p1+eps(max(1,abs(p2))), value=wD(index); else
    value=wD(index-1)+(target-p1)*(wD(index)-wD(index-1))/(p2-p1);
end
end
function [low,high]=wilson(successes,n,z)
p=successes./n; d=1+z^2/n; center=(p+z^2/(2*n))./d;
half=z.*sqrt(p.*(1-p)./n+z^2/(4*n^2))./d;
low=max(0,center-half); high=min(1,center+half);
end
function value=safeRatio(a,b)
if isfinite(a)&&isfinite(b)&&b>0, value=a/b; else, value=nan; end
end

function item=packCase(def,line,annular,equal)
item=emptyCase(); item.caseId=def.caseId; item.caseLabel=def.caseLabel;
item.series=def.series; item.xValue=def.xValue; item.overrides=def.overrides;
item.line=line; item.annular=annular; item.equalPerformance=equal;
end
function data=makePlotData(cases)
n=numel(cases); data=struct('caseId',strings(n,1),'caseLabel',strings(n,1), ...
    'series',strings(n,1),'x',nan(n,1),'kappa',nan(n,1), ...
    'kappaLow95',nan(n,1),'kappaHigh95',nan(n,1), ...
    'lineMinWD_mrad',nan(n,1),'ringMinWD_mrad',nan(n,1), ...
    'robustAdvantage',false(n,1),'advantageDisappears',false(n,1));
for k=1:n
    e=cases(k).equalPerformance; data.caseId(k)=cases(k).caseId;
    data.caseLabel(k)=cases(k).caseLabel; data.series(k)=cases(k).series;
    data.x(k)=cases(k).xValue; data.kappa(k)=e.kappa;
    data.kappaLow95(k)=e.kappaLow95; data.kappaHigh95(k)=e.kappaHigh95;
    data.lineMinWD_mrad(k)=e.lineMinWD_mrad;
    data.ringMinWD_mrad(k)=e.ringMinWD_mrad;
    data.robustAdvantage(k)=e.robustAdvantage;
    data.advantageDisappears(k)=e.advantageDisappears;
end
end
function summary=makeSummary(d)
summary=table(d.caseId,d.series,d.x,d.lineMinWD_mrad,d.ringMinWD_mrad, ...
    d.kappa,d.kappaLow95,d.kappaHigh95,d.robustAdvantage, ...
    d.advantageDisappears,'VariableNames',{'CaseId','Series', ...
    'NormalizedParameter','LineMinWD_mrad','AnnularMinWD_mrad','Kappa', ...
    'KappaLow95','KappaHigh95','RobustAdvantage','AdvantageDisappears'});
summary=sortrows(summary,{'Series','NormalizedParameter'});
end
function result=packResult(cfg,definitions,cases,plotData,summary,elapsed)
metadata=struct('description',"Independent receiver-noise paired ROE Monte Carlo", ...
    'studyName',"receiver_noise_quality", ...
    'panelTitle',"Receiver-noise suppression", ...
    'xLabel',"\eta_R (normalized receiver quality)", ...
    'xDescription',"Product of bandwidth, gate and dark-count quality factors", ...
    'generatedAt',string(datetime('now','Format','yyyy-MM-dd HH:mm:ss Z')), ...
    'kappaDefinition',"min(wD_annular)/min(wD_line)", ...
    'historicalDataLoaded',false,'independentEntryPoint',true);
result=struct('metadata',metadata,'config',cfg,'caseDefinitions',definitions, ...
    'cases',cases,'plotData',plotData,'summary',summary,'elapsed_s',elapsed);
end

function d=emptyDefinition()
d=struct('caseId',"",'caseLabel',"",'series',"",'xValue',nan,'overrides',struct());
end
function c=emptyCase()
c=struct('caseId',"",'caseLabel',"",'series',"",'xValue',nan, ...
    'overrides',struct(),'line',struct(),'annular',struct(),'equalPerformance',struct());
end
function verifySimulator(moduleDir)
names={'ea_sens_default_config.m','ea_sens_build_scenario_bank.m','ea_sens_simulate_case.m'};
for k=1:numel(names)
    assert(isfile(fullfile(moduleDir,names{k})),'receiver_noise_data:MissingSimulator', ...
        'Required simulator file is missing: %s',fullfile(moduleDir,names{k}));
end
end
function added=addSimulatorPath(moduleDir)
entries=string(strsplit(path,pathsep)); added=~any(strcmpi(entries,moduleDir));
if added, addpath(moduleDir); end
end
function restorePath(moduleDir,added), if added, rmpath(moduleDir); end, end
function tf=isTextScalar(x), tf=ischar(x)||(isstring(x)&&isscalar(x)); end
function tf=optionalInteger(x), tf=isempty(x)||(isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&fix(x)==x); end
function tf=optionalVector(x), tf=isempty(x)||(isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>0)); end
function tf=nonnegativeScalar(x), tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0; end
function tf=probabilityScalar(x), tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<=1; end
function tf=nonnegativeInteger(x), tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&fix(x)==x; end

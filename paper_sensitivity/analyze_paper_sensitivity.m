function results = analyze_paper_sensitivity(options)
%ANALYZE_PAPER_SENSITIVITY Paired OAT width uncertainty and joint contrasts.
arguments
    options.DataDirectory (1,1) string = ""
    options.OutputDirectory (1,1) string = ""
    options.JointFile (1,1) string = ""
    options.Deadline (1,1) double {mustBePositive} = 15
    options.RequiredProbability (1,1) double {mustBePositive} = .6
    options.BootstrapReplicates (1,1) double {mustBeInteger,mustBePositive} = 1000
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 20260928
end
root=fileparts(fileparts(mfilename('fullpath')));
if options.DataDirectory=="",options.DataDirectory=fullfile(root,'para_sensative_N3000_range2');end
if options.OutputDirectory=="",options.OutputDirectory=fullfile(root,'paper_sensitivity','results');end
if options.JointFile=="",options.JointFile=fullfile(options.OutputDirectory,'joint_static_N2048','joint_results.mat');end
if ~isfolder(options.OutputDirectory),mkdir(options.OutputDirectory);end
saved=load(fullfile(options.DataDirectory,'para_summary.mat'));
m=saved.manifest;w=m.widthsMrad;n=m.N;
assert(m.complete && ~m.truncated,'Sensitivity:Source','Incomplete source data.');
for sourceIndex=1:numel(m.signature.sourceNames)
    assert(strcmp(fileread(fullfile(root,m.signature.sourceNames(sourceIndex))),m.signature.sourceText{sourceIndex}), ...
        'Sensitivity:SourceChanged','The original simulation source has changed since the OAT run.');
end
oldRng=rng;cleanup=onCleanup(@()rng(oldRng)); %#ok<NASGU>
rng(options.Seed,'twister');
weights=bootstrapWeights(n,options.BootstrapReplicates);
rep=options.BootstrapReplicates;
records=cell(numel(m.cases),1);
responses=cell(numel(m.cases),1);
bootKappa=nan(rep,numel(m.cases));
baseIndex=find([m.cases.parameterIndex]==1 & [m.cases.value]==m.baseline(1),1);
order=[baseIndex,setdiff(1:numel(m.cases),baseIndex,'stable')];
codes='LR';
for j=order
    if isequal(m.cases(j).dynamicFiles,m.cases(baseIndex).dynamicFiles) && j~=baseIndex
        r=records{baseIndex};r.Case=m.cases(j).id;r.ParameterIndex=m.cases(j).parameterIndex;r.Value=m.cases(j).value;
        records{j}=r;responses{j}=responses{baseIndex};bootKappa(:,j)=bootKappa(:,baseIndex);continue
    end
    event=false(n,numel(w),2);
    for d=1:numel(w)
        for b=1:2
            code=codes(b);
            raw=load(m.cases(j).dynamicFiles(d,b),['detect_' code],['first_time_' code],'meta');
            hit=raw.(['detect_' code]);time=raw.(['first_time_' code]);
            assert(isequal(size(hit),[n 1]) && raw.meta.N==n && raw.meta.widthMrad==w(d) ...
                && isequal(raw.meta.physics,m.cases(j).physics) && raw.meta.speed==m.cases(j).speed, ...
                'Sensitivity:Metadata','Raw data do not match manifest.');
            event(:,d,b)=hit & time<=options.Deadline;
        end
    end
    p=squeeze(mean(event,1));
    [minimum,gridMinimum]=sensitivity_min_width(w,p,options.RequiredProbability);
    bootstrapP=reshape(double(reshape(event,n,[])).'*weights/n,numel(w),[]);
    minimumBoot=reshape(sensitivity_min_width(w,bootstrapP,options.RequiredProbability),2,rep).';
    bootKappa(:,j)=minimumBoot(:,2)./minimumBoot(:,1);
    ci=nan(2,2);
    for b=1:2
        finite=isfinite(minimumBoot(:,b));
        if isfinite(minimum(b)) && any(finite),ci(b,:)=quantile(minimumBoot(finite,b),[.025 .975]);end
    end
    finite=isfinite(bootKappa(:,j));kci=[nan nan];
    if all(isfinite(minimum)) && any(finite),kci=quantile(bootKappa(finite,j),[.025 .975]);end
    r=struct('Case',m.cases(j).id,'ParameterIndex',m.cases(j).parameterIndex, ...
        'Value',m.cases(j).value,'LineWidth',minimum(1),'RingWidth',minimum(2), ...
        'LineGridWidth',gridMinimum(1),'RingGridWidth',gridMinimum(2), ...
        'LineLow',ci(1,1),'LineHigh',ci(1,2),'RingLow',ci(2,1),'RingHigh',ci(2,2), ...
        'LineFeasibleBootstrap',mean(isfinite(minimumBoot(:,1))), ...
        'RingFeasibleBootstrap',mean(isfinite(minimumBoot(:,2))), ...
        'Kappa',minimum(2)/minimum(1),'KappaLow',kci(1),'KappaHigh',kci(2), ...
        'BothFeasibleBootstrap',mean(finite),'LineMaxP',max(p(:,1)),'RingMaxP',max(p(:,2)));
    records{j}=r;responses{j}=p;
    fprintf('OAT %s: minimum L/R=%s, feasible bootstrap=%s\n',r.Case,mat2str(minimum,4),mat2str(mean(isfinite(minimumBoot)),3));
end
oat=struct2table(vertcat(records{:}));
% Cross-check all imported empirical curves with the earlier full summary.
requirementIndex=find(saved.options.Requirements(:,1)==options.Deadline,1);
if ~isempty(requirementIndex)
    for j=1:numel(m.cases)
        assert(max(abs(responses{j}-saved.curves{j}.pd(:,:,requirementIndex)),[],'all')<1e-12, ...
            'Sensitivity:Reconstruction','Imported probabilities differ from saved summary.');
    end
end
% Secondary requirements are fixed before plotting, so conclusions need not
% depend on the single selected 60% criterion. Export all, including failures.
targets=[.5 .6 .7 .8];secondary=cell(numel(m.cases)*numel(targets),1);k=0;
for j=1:numel(m.cases)
    for target=targets
        k=k+1;minimum=sensitivity_min_width(w,responses{j},target);
        secondary{k}=struct('Case',m.cases(j).id,'TargetP',target, ...
            'LineWidth',minimum(1),'RingWidth',minimum(2),'Kappa',minimum(2)/minimum(1));
    end
end
secondary=struct2table(vertcat(secondary{:}));
theory=checkPhysicalDirections(m,responses,saved.curves);
theory.MinimumWidthOrderViolations=zeros(1,4);
for p=1:4
    selected=oat(oat.ParameterIndex==p,:);
    width=[selected.LineWidth,selected.RingWidth];
    width(isnan(width))=inf;
    direction=1;
    if p<=2,direction=-1;end
    theory.MinimumWidthOrderViolations(p)=nnz(direction*diff(width,1,1)<-1e-12);
end
joint=table();jointMetadata=struct();
if isfile(options.JointFile)
    raw=load(options.JointFile,'result');a=raw.result;
    assert(all(a.complete),'Sensitivity:JointIncomplete','Joint experiment is still running.');
    if isfield(a,'staticComplete')
        assert(all(a.staticComplete),'Sensitivity:StaticIncomplete','The static extension is still running.');
    end
    staticN=size(a.staticHit,1);
    jointWeights=bootstrapWeights(a.options.N,rep);
    staticWeights=bootstrapWeights(staticN,rep);
    jr=cell(size(a.values,1),1);
    jointGain=zeros(rep,size(a.values,1));
    jointReduction=jointGain;
    for j=1:size(a.values,1)
        hit=double(a.movingHit(:,:,j));blind=double(~a.staticHit(:,:,j));
        pd=mean(hit);b=mean(blind);
        gain=(hit(:,2)-hit(:,1)).'*jointWeights/a.options.N;
        reduction=(blind(:,1)-blind(:,2)).'*staticWeights/staticN;
        jointGain(:,j)=gain.';
        jointReduction(:,j)=reduction.';
        pci=quantile(gain,[.025 .975]);bci=quantile(reduction,[.025 .975]);
        jr{j}=struct('Scenario',j-1,'Energy_uJ',a.values(j,1)*1e6,'Reflectivity',a.values(j,2), ...
            'Speed_mps',a.values(j,3),'Alpha',a.values(j,4),'LinePd',pd(1),'RingPd',pd(2), ...
            'PdGain',pd(2)-pd(1),'PdGainLow',pci(1),'PdGainHigh',pci(2), ...
            'LineBlind',b(1),'RingBlind',b(2),'BlindReduction',b(1)-b(2), ...
            'BlindReductionLow',bci(1),'BlindReductionHigh',bci(2));
    end
    joint=struct2table(vertcat(jr{:}));
    % A single simultaneous band covers the 24 perturbations, the separate
    % baseline, and both response contrasts. Baseline is not a LHS draw.
    boot=[jointGain,jointReduction];
    estimates=[joint.PdGain;joint.BlindReduction].';
    standardError=std(boot,0,1);
    standardized=abs((boot-estimates)./max(standardError,eps));
    critical=quantile(max(standardized,[],2),.95);
    deltaPd=critical*std(jointGain,0,1).';
    deltaBlind=critical*std(jointReduction,0,1).';
    joint.PdSimLow=joint.PdGain-deltaPd;
    joint.PdSimHigh=joint.PdGain+deltaPd;
    joint.BlindSimLow=joint.BlindReduction-deltaBlind;
    joint.BlindSimHigh=joint.BlindReduction+deltaBlind;
    jointMetadata=struct('N',a.options.N,'StaticN',staticN,'ScenarioCount',a.options.ScenarioCount, ...
        'HalfRange',a.options.RelativeHalfRange,'WidthMrad',a.options.WidthMrad, ...
        'Deadline',a.options.Deadline,'Seed',a.options.Seed, ...
        'BaselineSeparate',true,'AllComplete',all(a.complete), ...
        'SimultaneousCriticalValue',critical,'SimultaneousResponses',2);
    writetable(joint,fullfile(options.OutputDirectory,'joint_contrasts.csv'));
end
results=struct('options',options,'manifest',m,'oat',oat,'secondary',secondary, ...
    'joint',joint,'jointMetadata',jointMetadata,'responses',{responses}, ...
    'bootstrapKappa',bootKappa,'theory',theory,'rawCurvesMatchSaved',true);
writetable(oat,fullfile(options.OutputDirectory,'oat_minimum_width.csv'));
writetable(secondary,fullfile(options.OutputDirectory,'secondary_requirements.csv'));
save(fullfile(options.OutputDirectory,'paper_sensitivity.mat'),'results','-v7.3');
fid=fopen(fullfile(options.OutputDirectory,'physical_checks.json'),'w');
fprintf(fid,'%s',jsonencode(theory,PrettyPrint=true));fclose(fid);
end

function weights=bootstrapWeights(n,reps)
weights=zeros(n,reps);
for k=1:reps,weights(:,k)=histcounts(randi(n,n,1),.5:1:n+.5).';end
end

function report=checkPhysicalDirections(m,p,curves)
report=struct('Method','Raw paired curves; no monotonicity correction', ...
    'energyPdViolations',0,'reflectivityPdViolations',0,'atmospherePdViolations',0, ...
    'energyBlindViolations',0,'reflectivityBlindViolations',0,'atmosphereBlindViolations',0);
names={'energy','reflectivity','atmosphere'};
parameters=[1 2 4];
for k=1:3
    indices=find([m.cases.parameterIndex]==parameters(k));
    signExpected=1;if parameters(k)==4,signExpected=-1;end
    for z=2:numel(indices)
        diffP=signExpected*(p{indices(z)}-p{indices(z-1)});
        diffB=signExpected*(curves{indices(z)}.blind-curves{indices(z-1)}.blind);
        report.([names{k} 'PdViolations'])=report.([names{k} 'PdViolations'])+nnz(diffP<-1e-12);
        report.([names{k} 'BlindViolations'])=report.([names{k} 'BlindViolations'])+nnz(diffB>1e-12);
    end
end
report.SpeedNote='No universal monotonicity is imposed on moving-target performance or relative gain.';
end

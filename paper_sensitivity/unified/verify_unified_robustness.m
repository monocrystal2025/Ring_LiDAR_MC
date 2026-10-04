function report = verify_unified_robustness(options)
%VERIFY_UNIFIED_ROBUSTNESS Verify coverage of the user's three equal metrics.
arguments
    options.Directory (1,1) string = ""
end
folder=fileparts(mfilename('fullpath'));
if options.Directory=="",options.Directory=fullfile(folder,'results');end
saved=load(fullfile(options.Directory,'analysis.mat'),'analysis');a=saved.analysis;r=a.design;
root=fileparts(fileparts(folder));
sourceFiles=[fullfile(folder,"run_unified_robustness.m"),fullfile(root,"MC_line_snr.m"), ...
    fullfile(root,"MC_ring_snr.m"),fullfile(root,"append_and_check_photon_window.m"),fullfile(root,"init_UAV.m")];
for sourceIndex=1:numel(sourceFiles)
    assert(strcmp(fileread(sourceFiles(sourceIndex)),r.signature.sources{sourceIndex}), ...
        'Unified:Sources','Simulation source differs from its saved provenance.');
end
k=find(r.options.ReferenceWidths==125,1);main=a.point(:,:,k);
assert(all(r.complete,'all'),'Unified:Complete','An intended simulation configuration is incomplete.');
assert(isequal(size(main),[13 3]),'Unified:Design','The three metrics must use the same 13 physical configurations.');
assert(all(isfinite(main),'all'),'Unified:Finite','The primary figure still has unattained/undefined values.');
assert(all(isfinite(a.ci(:,:,:,k)),'all'),'Unified:Intervals','The primary intervals are incomplete.');
assert(all(a.feasibility(:,3)>=.975),'Unified:Feasibility','Threshold reachability uncertainty must be resolved.');
assert(max(a.brackets(:,:,2)-a.brackets(:,:,1),[],'all')<=5, ...
    'Unified:Resolution','A threshold crossing uses a grid gap larger than 5 mrad.');
assert(all(a.theory.TimelyViolations==0) && all(a.theory.StaticViolations==0), ...
    'Unified:Physics','Expected model direction is violated; investigate before delivery.');
assert(isequal(r.scale(r.factorIndex==0,:),ones(1,3)), 'Unified:Baseline','Missing common baseline.');
for p=1:3
    ids=find(r.factorIndex==p | r.factorIndex==0);
    assert(isequal(sort(r.scale(ids,p)).',r.options.Factors),'Unified:Levels','A perturbation level was omitted.');
end
files=["three_benefits_w125","three_benefits_w100","three_benefits_w150", ...
    "absolute_responses","wide_range_boundaries"];
for base=files
    for ext=[".png",".pdf",".fig"]
        info=dir(fullfile(folder,'figures',base+ext));
        assert(~isempty(info) && info.bytes>1000,'Unified:Artifacts','Missing figure artifact.');
    end
end
wide=readtable(fullfile(options.Directory,'wide_range_boundaries.csv'));
assert(height(wide)==20 && any(wide.BlindReduction<0) && any(isnan(wide.AngularReduction)), ...
    'Unified:Boundaries','Wide-range failures and reversals must remain disclosed.');
report=struct('ChecksPassed',true,'Configurations',size(main,1),'EqualMetrics',3, ...
    'MainFinitePoints',numel(main),'MainPositiveEstimates',sum(main>0,1), ...
    'MainPointwiseLowerAboveZero',sum(a.ci(:,:,1,k)>0,1), ...
    'MainSimultaneousLowerAboveZero',sum(a.sim(:,:,1,k)>0,1), ...
    'MinimumJointBootstrapReachability',min(a.feasibility(:,3)), ...
    'NMoving',r.options.N,'NStatic',r.options.StaticN,'PhysicalChecks',a.theory, ...
    'FullEarlierRangeRetained',true,'PrimaryReferenceWidth',125,'SupplementReferenceWidths',[100 150]);
fid=fopen(fullfile(options.Directory,'verification.json'),'w');fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));fclose(fid);
disp(report);
end

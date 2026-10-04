function report=verify_expanded_robustness()
%VERIFY_EXPANDED_ROBUSTNESS Verify complete symmetric parameter domains.
folder=fileparts(mfilename('fullpath'));s=load(fullfile(folder,'results','analysis.mat'),'analysis');a=s.analysis;r=a.design;
for k=1:numel(r.sourceNames)
assert(strcmp(fileread(r.sourceNames(k)),r.signature.sources{k}),'Expanded:Source','Simulation source changed.');
end
assert(all(a.halfRange(1:3)>=30),'Expanded:Range','A required range did not exceed twenty percent.');
for p=find(a.halfRange>0)
ids=find(a.chosen & (r.factorIndex==p | r.factorIndex==0));
assert(isequal(sort(round(100*(r.scale(ids,p)-1))).',-a.halfRange(p):10:a.halfRange(p)),'Expanded:Symmetry','Incomplete symmetric interval.');
end
assert(all(isfinite(a.point(a.chosen,:)) & a.point(a.chosen,:)>0,'all'),'Expanded:Main','Missing/nonpositive main estimate.');
assert(all(isfinite(a.ci(a.chosen,:,:)),'all') && all(a.feasibility(a.chosen)>=.975),'Expanded:CI','Incomplete intervals.');
assert(max(a.brackets(a.chosen,:,2)-a.brackets(a.chosen,:,1),[],'all')<=5,'Expanded:Grid','Threshold grid is too coarse.');
assert(all(a.theory.TimelyViolations==0) && all(a.theory.StaticViolations==0),'Expanded:Physics','Unexpected physical direction.');
assert(all(r.complete(r.widths==125,:)),'Expanded:Screen','Missing candidate reference response.');
for name=["expanded_three_benefits","expanded_absolute_responses","candidate_boundaries"]
for ext=[".png",".pdf",".fig"]
info=dir(fullfile(folder,'figures',name+ext));assert(~isempty(info)&&info.bytes>1000,'Expanded:Export','Missing figure.');
end
end
errors=zeros(1,3);widths=[5 60 125];
for k=1:3
d=find(r.manifest.widthsMrad==widths(k));original=readNPY(char(r.manifest.scanFiles(d)));
generated=expanded_spiral(r.manifest.f,2000*tan(widths(k)*1e-3/2),r.manifest.omega,1000);
assert(isequal(size(original),size(generated)),'Expanded:Trajectory','Baseline size differs.');
errors(k)=max(abs(original-generated),[],'all');assert(errors(k)<1e-12,'Expanded:Trajectory','Baseline scan differs.');
end
report=struct('ChecksPassed',true,'HalfRangePercent',a.halfRange,'CandidateConfigurations',size(r.scale,1),'SelectedConfigurations',nnz(a.chosen),'MainFinitePoints',nnz(isfinite(a.point(a.chosen,:))),'PointwiseLowerAboveZero',sum(a.ci(a.chosen,:,1)>0,1),'SimultaneousLowerAboveZero',sum(a.sim(a.chosen,:,1)>0,1),'MinimumBootstrapReachability',min(a.feasibility(a.chosen)),'BaselineTrajectoryMaxErrors',errors,'PhysicalChecks',a.theory,'AllReferenceCandidatesRetained',true);
fid=fopen(fullfile(folder,'results','verification.json'),'w');fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));fclose(fid);disp(report);
end

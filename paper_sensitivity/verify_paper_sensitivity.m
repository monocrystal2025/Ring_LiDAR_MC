function report = verify_paper_sensitivity()
%VERIFY_PAPER_SENSITIVITY Check data provenance and delivered analysis outputs.
folder=fileparts(mfilename('fullpath'));
data=load(fullfile(folder,'results','paper_sensitivity.mat'),'results');
r=data.results;
assert(r.rawCurvesMatchSaved && height(r.oat)==20,'Sensitivity:Verify','Incomplete OAT reconstruction.');
joint=load(r.options.JointFile,'result');
joint=joint.result;
assert(all(joint.complete),'Sensitivity:Verify','Incomplete dynamic configurations.');
if isfield(joint,'staticComplete')
    assert(all(joint.staticComplete),'Sensitivity:Verify','Incomplete static extension.');
    originalFile=fullfile(joint.options.OutputDirectory,'joint_results.mat');
    original=load(originalFile,'result');
    original=original.result;
    assert(isequaln(original.movingHit,joint.movingHit) && isequaln(original.movingTime,joint.movingTime), ...
        'Sensitivity:Verify','Static extension changed dynamic observations.');
    assert(isequal(original.staticHit,joint.staticHit(1:size(original.staticHit,1),:,:)), ...
        'Sensitivity:Verify','Original static observations were not retained.');
end
scale=joint.scale(2:end,:);
half=joint.options.RelativeHalfRange;
u=(scale-(1-half))/(2*half);
count=size(scale,1);
assert(all(u>0 & u<1,'all'),'Sensitivity:Verify','Perturbation outside the specified box.');
strata=floor(count*u)+1;
assert(isequal(sort(strata,1),repmat((1:count).',1,4)), ...
    'Sensitivity:Verify','The parameter design is not a Latin hypercube.');
assert(size(joint.movingHit,1)==r.jointMetadata.N ...
    && size(joint.staticHit,1)==r.jointMetadata.StaticN, ...
    'Sensitivity:Verify','Reported sample sizes disagree with observations.');
times=joint.movingTime(joint.movingHit);
assert(all(isfinite(times) & times>=0 & times<=joint.options.Deadline), ...
    'Sensitivity:Verify','A detected moving target exceeds the prescribed deadline.');
pd=squeeze(mean(joint.movingHit,1)).';
blind=squeeze(mean(~joint.staticHit,1)).';
assert(max(abs(pd-[r.joint.LinePd r.joint.RingPd]),[],'all')<1e-12 ...
    && max(abs(blind-[r.joint.LineBlind r.joint.RingBlind]),[],'all')<1e-12, ...
    'Sensitivity:Verify','Exported response values disagree with raw observations.');
assert(all(isnan(r.oat.LineLow(isnan(r.oat.LineWidth)))) ...
    && all(isnan(r.oat.RingLow(isnan(r.oat.RingWidth)))), ...
    'Sensitivity:Verify','An unattained point is assigned a misleading finite interval.');
files=["figures_final/paper_parameter_sensitivity.png", ...
    "figures_final/paper_parameter_sensitivity.pdf","figures_final/paper_parameter_sensitivity.fig", ...
    "results/oat_minimum_width.csv","results/secondary_requirements.csv", ...
    "results/joint_contrasts.csv","manuscript_final.md","README.md"];
for j=1:numel(files)
    info=dir(fullfile(folder,files(j)));
    assert(~isempty(info) && info.bytes>100,'Sensitivity:Verify','Missing or empty deliverable: %s',files(j));
end
imageInfo=imfinfo(fullfile(folder,files(1)));
assert(imageInfo.Width>2000 && imageInfo.Height>1500,'Sensitivity:Verify','Raster export is too small.');
report=struct('AllChecksPassed',true,'OATConfigurations',height(r.oat), ...
    'OATTargetsPerConfiguration',r.manifest.N,'JointPerturbations',count, ...
    'JointMovingN',r.jointMetadata.N,'JointStaticN',r.jointMetadata.StaticN, ...
    'OriginalObservationsRetained',true,'LatinHypercubeValid',true, ...
    'RawResponseReconstructionPassed',true,'FigureSize',[imageInfo.Width imageInfo.Height], ...
    'PhysicalDirectionChecks',r.theory,'Deliverables',files);
fid=fopen(fullfile(folder,'results','verification.json'),'w');
fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));
fclose(fid);
disp(report);
end

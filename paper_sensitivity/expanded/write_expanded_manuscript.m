function write_expanded_manuscript()
%WRITE_EXPANDED_MANUSCRIPT Export measured operating ranges and limitations.
folder=fileparts(mfilename('fullpath'));
saved=load(fullfile(folder,'results','analysis.mat'),'analysis');a=saved.analysis;
fid=fopen(fullfile(folder,'manuscript.md'),'w');cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# Expanded symmetric operating domains\n\n');
fprintf(fid,'Half ranges for energy, reflectivity, extinction and scan speed: %s percent. All sampled levels are multiples of 10 percent.\n\n',mat2str(a.halfRange));
ids=find(a.chosen);limits=[min(a.point(ids,:),[],1);max(a.point(ids,:),[],1)];
names={'Required-divergence reduction','Relative blind-zone reduction','Effective-pulse increase'};
fprintf(fid,'|Metric|Benefit range|Pointwise lower > 0|Simultaneous lower > 0|\n|---|---|---|---|\n');
for p=1:3
fprintf(fid,'|%s|%.2f to %.2f percent|%d/%d|%d/%d|\n',names{p},limits(1,p),limits(2,p),sum(a.ci(ids,p,1)>0),numel(ids),sum(a.sim(ids,p,1)>0),numel(ids));
end
fprintf(fid,'\nThese are exploratory domains selected from the disclosed candidate sweep, not preregistered tolerance bounds. Intervals describe target-sampling uncertainty conditional on selected configurations; they exclude domain selection, angular-grid and model uncertainty. Every candidate remains in all_candidates.csv, with decisions in range_decisions.csv. No claim is made about every continuous parameter value or arbitrary joint perturbations.\n\n');
fprintf(fid,'Energy/reflectivity/speed candidates extend to +/-50 percent; extinction to +/-90 percent. Levels are separated by 10 percent. The largest tested contiguous symmetric domain with positive finite estimates and at least 97.5 percent bootstrap joint reachability is displayed. Statistical significance is not an inclusion rule. All interior levels remain.\n\n');
fprintf(fid,'The same 512 moving and 2048 stationary targets are used throughout. Required divergence satisfies detection fraction >=60 percent within 15 s. Blind-zone and successful-window pulse statistics use 125 mrad. Existing identical trials are reused after source checks. Scan-speed perturbations regenerate every pulse direction, boundary circles, phase support, cycle and accumulation window at fixed PRF and pulse energy. Static coverage is measured over one complete cycle at each speed, whose duration varies.\n\n');
fprintf(fid,'Figure: Three parallel annular-beam benefits over expanded symmetric parameter perturbations. Ranges are shown above columns. Positive values favor annular illumination. Whiskers are pointwise 95 percent paired-bootstrap intervals; diamonds mark the shared baseline. Supplemental figures retain all candidates and absolute responses.\n');
end

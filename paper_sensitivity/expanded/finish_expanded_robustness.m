function finish_expanded_robustness()
%FINISH_EXPANDED_ROBUSTNESS Extend unattained searches, then export evidence.
folder=fileparts(mfilename('fullpath'));addpath(folder,fileparts(folder));
saved=load(fullfile(folder,'results','design.mat'),'result');r=saved.result;
plan=load(fullfile(folder,'results','angular_search_plan.mat'),'selected');
assert(all(r.complete(r.widths<=150,plan.selected),'all'),'Expanded:Incomplete','Wait for the full initial campaign.');
unattained=[];
for j=plan.selected
    probabilities=zeros(30,2);
    for d=1:30
        raw=load(r.files(d,j),'hit','time');probabilities(d,:)=mean(raw.hit & raw.time<=15);
    end
    if any(max(probabilities,[],1)<.65),unattained(end+1)=j;end %#ok<AGROW>
end
if ~isempty(unattained)
    fprintf('Extending unattained cases to 400 mrad: %s\n',mat2str(unattained));
    run_expanded_robustness(CaseIDs=unattained,Widths=155:5:400);
end
analyze_expanded_robustness;
plot_expanded_robustness;
write_expanded_manuscript;
verify_expanded_robustness;
fprintf('EXPANDED_DELIVERABLES_COMPLETE\n');
end


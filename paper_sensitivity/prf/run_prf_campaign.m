function run_prf_campaign()
%RUN_PRF_CAMPAIGN Reference screen then complete angular searches.
folder=fileparts(mfilename('fullpath'));addpath(folder);
run_prf_sensitivity(Widths=125);
analyze_prf_sensitivity;
% No candidates are filtered on their advantage sign.
run_prf_sensitivity(Widths=5:5:150);
s=load(fullfile(folder,'results','design.mat'),'result');r=s.result;extend=[];
for j=1:numel(r.frequencies)
    pd=zeros(30,2);
    for d=1:30
        raw=load(r.files(d,j),'hit','time');pd(d,:)=mean(raw.hit & raw.time<=15);
    end
    if any(max(pd,[],1)<.65),extend(end+1)=r.frequencies(j);end %#ok<AGROW>
end
if ~isempty(extend),run_prf_sensitivity(Frequencies=extend,Widths=155:5:400);end
analyze_prf_sensitivity;
fprintf('PRF_CAMPAIGN_COMPLETE\n');
end

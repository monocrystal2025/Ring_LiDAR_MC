function run_expanded_campaign()
%RUN_EXPANDED_CAMPAIGN Screen the fixed reference, then complete angular grids.
folder=fileparts(mfilename('fullpath'));addpath(folder);
r=run_expanded_robustness(Widths=125);
reference=find(r.widths==125);nc=size(r.scale,1);summary=zeros(nc,5);
for j=1:nc
    raw=load(r.files(reference,j));blind=mean(~raw.staticHit);pulse=mean(raw.pulse,1,'omitnan');
    summary(j,:)=[j,100*(1-blind(2)/blind(1)),100*(pulse(2)/pulse(1)-1),mean(raw.hit(:,1)&raw.time(:,1)<=15),mean(raw.hit(:,2)&raw.time(:,2)<=15)];
end
t=array2table(summary,'VariableNames',{'Case','BlindGain','PulseGain','LinePd125','RingPd125'});
t.ParameterIndex=r.factorIndex;t.EnergyScale=r.scale(:,1);t.ReflectivityScale=r.scale(:,2);
t.AlphaScale=r.scale(:,3);t.OmegaScale=r.scale(:,4);
writetable(t,fullfile(folder,'results','reference_screen.csv'));
% A non-positive reference benefit already violates the requested operating
% domain criterion. Retain it in the table instead of simulating redundant
% angular searches for that rejected case. This is not a significance gate.
selected=find(all(isfinite(summary(:,2:3)) & summary(:,2:3)>0,2)).';
save(fullfile(folder,'results','angular_search_plan.mat'),'selected','t');
fprintf('REFERENCE_SCREEN_COMPLETE: angular searches for %d/%d candidates\n',numel(selected),nc);
run_expanded_robustness(CaseIDs=selected,Widths=[5:5:100 105:5:150]);
fprintf('EXPANDED_CAMPAIGN_COMPLETE\n');
end

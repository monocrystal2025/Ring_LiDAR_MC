function study=compare_prf_criteria()
%COMPARE_PRF_CRITERIA Compare declared deadlines/thresholds without imputing failures.
folder=fileparts(mfilename('fullpath'));addpath(fileparts(folder));
s=load(fullfile(folder,'results','analysis.mat'),'a');original=s.a;r=original.design;
out=fullfile(folder,'criteria_study');if ~isfolder(out),mkdir(out);end
times=[10 15];targets=[.2 .3 .4 .5 .6];reps=1000;n=r.bank.options.N;
previous=rng;cleanup=onCleanup(@()rng(previous));rng(20260928,'twister');
wm=zeros(n,reps);for b=1:reps,wm(:,b)=histcounts(randi(n,n,1),.5:1:n+.5).';end
minimum=nan(10,2,numel(times),numel(targets));gain=nan(10,numel(times),numel(targets));
ci=nan(10,2,numel(times),numel(targets));reach=zeros(10,numel(times),numel(targets));
boot=nan(reps,10,numel(times),numel(targets));reference=zeros(10,2,6);referenceTimes=[15 20 25 30 40 60];
for j=1:10
    first=find(~r.complete(:,j),1);last=numel(r.widths);if ~isempty(first),last=first-1;end
    assert(last>=30);hit=false(n,last,2);when=nan(n,last,2);
    for d=1:last
        raw=load(r.files(d,j),'meta','hit','time');
        assert(raw.meta.DynamicLimit/r.frequencies(j)>=max(times));
        hit(:,d,:)=reshape(raw.hit,n,1,2);when(:,d,:)=reshape(raw.time,n,1,2);
    end
    for ti=1:numel(times)
        event=hit & when<=times(ti);pd=squeeze(mean(event,1));
        bp=reshape(double(reshape(event,n,[])).'*wm/n,last,[]);
        for qi=1:numel(targets)
            w=sensitivity_min_width(r.widths(1:last),pd,targets(qi));
            bw=reshape(sensitivity_min_width(r.widths(1:last),bp,targets(qi)),2,reps);
            values=100*(1-bw(2,:)./bw(1,:));good=isfinite(values);
            minimum(j,:,ti,qi)=w;gain(j,ti,qi)=100*(1-w(2)/w(1));
            reach(j,ti,qi)=mean(good);boot(:,j,ti,qi)=values.';
            if mean(good)>=.975,ci(j,:,ti,qi)=quantile(values(good),[.025 .975]);end
        end
    end
    raw=load(r.files(r.widths==125,j),'meta','hit','time');assert(isinf(raw.meta.DynamicLimit));
    for k=1:numel(referenceTimes),reference(j,:,k)=mean(raw.hit & raw.time<=referenceTimes(k));end
end
rows=cell(0,1);details=cell(0,1);
for ti=1:numel(times)
    for qi=1:numel(targets)
        y=gain(:,ti,qi);finite=isfinite(y);
        rows{end+1}=struct('DeadlineSeconds',times(ti),'TargetProbability',targets(qi),'BothFeasible',nnz(finite), ...
            'Configurations',10,'CompleteCI',nnz(all(isfinite(ci(:,:,ti,qi)),2)), ...
            'MinimumBootstrapReachability',min(reach(:,ti,qi)),'PositivePointEstimates',nnz(y>0), ...
            'MinimumFiniteGain',min(y,[],'omitnan'),'MaximumFiniteGain',max(y,[],'omitnan')); %#ok<AGROW>
        for j=1:10
            details{end+1}=struct('FrequencyKHz',r.frequencies(j)/1000,'DeadlineSeconds',times(ti),'TargetProbability',targets(qi), ...
                'LineMinMrad',minimum(j,1,ti,qi),'RingMinMrad',minimum(j,2,ti,qi),'GainPercent',y(j), ...
                'Low',ci(j,1,ti,qi),'High',ci(j,2,ti,qi),'BootstrapReachability',reach(j,ti,qi)); %#ok<AGROW>
        end
    end
end
summary=struct2table([rows{:}]);detail=struct2table([details{:}]);
writetable(summary,fullfile(out,'criterion_comparison.csv'));writetable(detail,fullfile(out,'angular_results.csv'));
study=struct('times',times,'targets',targets,'minimum',minimum,'gain',gain,'ci',ci,'reach',reach,'bootstrap',boot, ...
    'referenceTimes',referenceTimes,'referencePd',reference,'summary',summary,'detail',detail);
save(fullfile(out,'study.mat'),'study','-v7.3');disp(summary);
% Check the original 15 s / 60% result is reproduced by the common code.
assert(isequaln(minimum(:,:,2,5),original.absolute(:,:,1)));
assert(max(abs(gain(:,2,5)-original.point(:,1)),[],'omitnan')<1e-12);
assert(isequaln(squeeze(ci(:,:,2,5)),squeeze(original.ci(:,1,:))));
% Export explicit alternatives. The original analysis and data stay intact.
for q=[.3 .5]
    qi=find(abs(targets-q)<1e-10);ti=2;a=original;
    a.criterion=struct('DeadlineSeconds',15,'TargetProbability',q,'Exploratory',true);
    a.absolute(:,:,1)=minimum(:,:,ti,qi);a.point(:,1)=gain(:,ti,qi);
    a.ci(:,1,:)=reshape(ci(:,:,ti,qi),10,1,2);a.bootstrap(:,:,1)=boot(:,:,ti,qi);a.reachability=reach(:,ti,qi);
    a.table.LineMinMrad=a.absolute(:,1,1);a.table.RingMinMrad=a.absolute(:,2,1);
    a.table.AngularGain=a.point(:,1);a.table.AngularLow=a.ci(:,1,1);a.table.AngularHigh=a.ci(:,1,2);
    a.table.BothBootstrapReachable=a.reachability;
    a.table.AngularStatus(:)="Both beams attainable";
    a.table.AngularStatus(~isfinite(a.point(:,1)))="Requirement unattained in searched grid";
    a.table.DeadlineSeconds=repmat(15,10,1);a.table.TargetProbability=repmat(q,10,1);
    targetFolder=fullfile(out,sprintf('T15_P%02d',round(100*q)));if ~isfolder(targetFolder),mkdir(targetFolder);end
    save(fullfile(targetFolder,'analysis.mat'),'a','-v7.3');writetable(a.table,fullfile(targetFolder,'summary.csv'));
    plot_prf_sensitivity(false,targetFolder,fullfile(targetFolder,'figures'));
end
fprintf('CRITERION_COMPARISON_COMPLETE\n');
end

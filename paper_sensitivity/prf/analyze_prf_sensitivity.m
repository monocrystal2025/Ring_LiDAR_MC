function a=analyze_prf_sensitivity()
%ANALYZE_PRF_SENSITIVITY Retain failures and quantify paired uncertainty.
folder=fileparts(mfilename('fullpath'));addpath(fileparts(folder));out=fullfile(folder,'results');
s=load(fullfile(out,'design.mat'),'result');r=s.result;n=r.bank.options.N;ns=r.bank.options.StaticN;
for k=1:numel(r.sourceNames),assert(strcmp(fileread(r.sourceNames(k)),r.signature.sources{k}),'PRF:Source','Source changed.');end
reps=1000;previous=rng;cleanup=onCleanup(@()rng(previous));rng(20260928,'twister');
wm=weights(n,reps);ws=weights(ns,reps);nf=numel(r.frequencies);
point=nan(nf,4);absolute=nan(nf,2,4);boot=nan(reps,nf,4);ci=nan(nf,4,2);
blindDiff=nan(nf,3);commonPulse=nan(nf,1);reach=nan(nf,1);searchMax=zeros(nf,1);
counts=zeros(nf,3);statuses=strings(nf,1);pd=nan(numel(r.widths),2,nf);cycles=zeros(nf,1);scanErrors=zeros(nf,1);
ref=find(r.widths==125);assert(all(r.complete(ref,:)),'PRF:Reference','Reference incomplete.');
for j=1:nf
    raw=load(r.files(ref,j));validate(raw,r,j,125);
    b=mean(~raw.staticHit);bb=double(~raw.staticHit).'*ws/ns;
    pm=mean(raw.pulse,1,'omitnan');pv=raw.pulse;pv(~raw.hit)=0;pb=(pv.'*wm)./(double(raw.hit).'*wm);
    assert(all(isfinite(raw.pulse(raw.hit))) && all(isnan(raw.pulse(~raw.hit))));
    for beam=1:2
        windows=raw.windows{beam};cnt=arrayfun(@(x)nnz(isfinite(x.signal_photons)&x.signal_photons>0),windows);
        assert(isequal(cnt(raw.hit(:,beam)),raw.pulse(raw.hit(:,beam),beam)));
    end
    common=all(raw.hit,2);cm=mean(raw.pulse(common,:),1);commonPulse(j)=100*(cm(2)/cm(1)-1);
    counts(j,:)=[sum(raw.hit),sum(common)];
    event=raw.hit & raw.time<=15;det=mean(event);db=double(event).'*wm/n;
    absolute(j,:,2)=b;absolute(j,:,3)=pm;absolute(j,:,4)=det;
    if b(1)>0,point(j,2)=100*(1-b(2)/b(1));boot(:,j,2)=100*(1-bb(2,:)./bb(1,:)).';end
    point(j,3)=100*(pm(2)/pm(1)-1);boot(:,j,3)=100*(pb(2,:)./pb(1,:)-1).';
    point(j,4)=100*(det(2)-det(1));boot(:,j,4)=100*(db(2,:)-db(1,:)).';
    blindDiff(j,:)=[100*(b(1)-b(2)),quantile(100*(bb(1,:)-bb(2,:)),[.025 .975])];
    if isfield(raw.meta,'CycleLength'),cycles(j)=raw.meta.CycleLength/r.frequencies(j);
    else,cycles(j)=(raw.meta.StaticLimit+1)/r.frequencies(j);end
    if isfield(raw.meta,'MaxRelativeScanSpeedError'),scanErrors(j)=raw.meta.MaxRelativeScanSpeedError;end
    first=find(~r.complete(:,j),1);last=numel(r.widths);if ~isempty(first),last=first-1;end
    statuses(j)="Angular search not yet complete";
    if last>=30
        searchMax(j)=r.widths(last);events=false(n,last,2);
        for d=1:last
            raw=load(r.files(d,j),'meta','hit','time');validate(raw,r,j,r.widths(d));
            events(:,d,:)=reshape(raw.hit & raw.time<=15,n,1,2);
        end
        probabilities=squeeze(mean(events,1));pd(1:last,:,j)=probabilities;
        minimum=sensitivity_min_width(r.widths(1:last),probabilities,.6);
        bp=reshape(double(reshape(events,n,[])).'*wm/n,last,[]);
        bw=reshape(sensitivity_min_width(r.widths(1:last),bp,.6),2,reps);
        reach(j)=mean(all(isfinite(bw),1));absolute(j,:,1)=minimum;
        point(j,1)=100*(1-minimum(2)/minimum(1));boot(:,j,1)=100*(1-bw(2,:)./bw(1,:)).';
        statuses(j)="Both beams attainable";
        if any(~isfinite(minimum)),statuses(j)="Requirement unattained in searched grid";end
    end
    for metric=1:4
        values=boot(:,j,metric);good=isfinite(values);
        if mean(good)>=.975 && isfinite(point(j,metric)),ci(j,metric,:)=quantile(values(good),[.025 .975]);end
    end
end
T=table(r.frequencies.'/1000,r.frequencies.'*r.manifest.baseline(1),cycles, ...
    absolute(:,1,1),absolute(:,2,1),point(:,1),ci(:,1,1),ci(:,1,2), ...
    100*absolute(:,1,2),100*absolute(:,2,2),point(:,2),ci(:,2,1),ci(:,2,2), ...
    blindDiff(:,1),blindDiff(:,2),blindDiff(:,3),absolute(:,1,3),absolute(:,2,3),point(:,3),ci(:,3,1),ci(:,3,2), ...
    100*absolute(:,1,4),100*absolute(:,2,4),point(:,4),ci(:,4,1),ci(:,4,2),commonPulse,reach,searchMax,statuses, ...
    'VariableNames',{'FrequencyKHz','AveragePowerW','CycleSeconds','LineMinMrad','RingMinMrad','AngularGain','AngularLow','AngularHigh', ...
    'LineBlindPercent','RingBlindPercent','BlindGain','BlindLow','BlindHigh','BlindDifferencePP','BlindDifferenceLow','BlindDifferenceHigh', ...
    'LinePulses','RingPulses','PulseGain','PulseLow','PulseHigh','LinePdPercent','RingPdPercent','PdDifferencePP','PdDifferenceLow','PdDifferenceHigh', ...
    'CommonPulseGain','BothBootstrapReachable','AngularSearchMax','AngularStatus'});
writetable(T,fullfile(out,'summary.csv'));
a=struct('design',r,'point',point,'absolute',absolute,'ci',ci,'bootstrap',boot,'blindDifference',blindDiff, ...
    'commonPulseGain',commonPulse,'reachability',reach,'counts',counts,'pd',pd,'table',T,'scanErrors',scanErrors);
save(fullfile(out,'analysis.mat'),'a','-v7.3');disp(T(:,{'FrequencyKHz','AngularGain','BlindGain','BlindDifferencePP','PulseGain','PdDifferencePP'}));
end

function validate(raw,r,j,width)
assert(raw.meta.Width==width && raw.meta.N==r.bank.options.N && raw.meta.Seed==r.bank.options.Seed);
if isfield(raw.meta,'Frequency'),assert(raw.meta.Frequency==r.frequencies(j));else,assert(r.frequencies(j)==r.manifest.f);end
expected=[r.manifest.baseline(1:2),r.manifest.baseline(4),r.manifest.baseline(4)/r.manifest.alphaBetaRatio];
assert(max(abs(raw.meta.Physics-expected))<1e-15);
if isfield(raw.meta,'Omega'),assert(raw.meta.Omega==r.manifest.omega);end
end

function w=weights(n,reps)
w=zeros(n,reps);for k=1:reps,w(:,k)=histcounts(randi(n,n,1),.5:1:n+.5).';end
end

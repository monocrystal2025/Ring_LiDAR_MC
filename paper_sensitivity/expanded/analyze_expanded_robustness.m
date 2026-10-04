function analysis = analyze_expanded_robustness()
%ANALYZE_EXPANDED_ROBUSTNESS Retain candidates and select contiguous symmetry.
folder=fileparts(mfilename('fullpath'));addpath(fileparts(folder));
out=fullfile(folder,'results');saved=load(fullfile(out,'design.mat'),'result');r=saved.result;

nr=size(r.scale,1);n=r.bank.options.N;ns=r.bank.options.StaticN;reps=1000;
old=rng;cleanup=onCleanup(@()rng(old));rng(20261002,'twister');
wm=bootstrapWeights(n,reps);ws=bootstrapWeights(ns,reps);
point=nan(nr,3);absolute=nan(nr,2,3);boot=nan(reps,nr,3);ci=nan(nr,3,2);
feasibility=zeros(nr,1);brackets=nan(nr,2,2);commonGain=nan(nr,1);
pd=nan(numel(r.widths),2,nr);static=false(ns,2,nr);counts=zeros(nr,3);
timely=cell(nr,1);status=strings(nr,1);searchMax=zeros(nr,1);
for j=1:nr
    ref=find(r.widths==125);assert(r.complete(ref,j),'Expanded:Reference','Incomplete reference screen.');
    raw=load(r.files(ref,j));validateRaw(raw,r,ref,j);
    static(:,:,j)=raw.staticHit;
    blind=double(~raw.staticHit);b=mean(blind);bb=blind.'*ws/ns;
    pulse=raw.pulse;valid=raw.hit;pv=pulse;pv(~valid)=0;
    assert(all(isnan(pulse(~valid))) && all(isfinite(pulse(valid))),'Expanded:Pulses','Invalid conditional pulse data.');
    for beam=1:2
        windows=raw.windows{beam};actual=arrayfun(@(x)nnz(isfinite(x.signal_photons) & x.signal_photons>0),windows);
        assert(isequal(actual(valid(:,beam)),pulse(valid(:,beam),beam)),'Expanded:Count','Photon-window count mismatch.');
        scores=arrayfun(@windowScore,windows(valid(:,beam)));
        assert(all(scores>=4-1e-8),'Expanded:SNR','Reported window below SNR threshold.');
    end
    pm=mean(pulse,1,'omitnan');pb=(pv.'*wm)./(double(valid).'*wm);
    common=all(valid,2);cm=mean(pulse(common,:),1);commonGain(j)=100*(cm(2)/cm(1)-1);counts(j,:)=[sum(valid),nnz(common)];
    absolute(j,:,2)=b;absolute(j,:,3)=pm;
    point(j,2:3)=[100*(1-b(2)/b(1)),100*(pm(2)/pm(1)-1)];
    boot(:,j,2)=100*(1-bb(2,:)./bb(1,:)).';boot(:,j,3)=100*(pb(2,:)./pb(1,:)-1).';

    if b(1)==0,point(j,2)=nan;boot(:,j,2)=nan;end
    firstMissing=find(~r.complete(:,j),1);last=numel(r.widths);if ~isempty(firstMissing),last=firstMissing-1;end
    if last<30
        assert(any(~isfinite(point(j,2:3)) | point(j,2:3)<=0),'Expanded:Incomplete','Eligible angular search is incomplete.');
        status(j)="Reference benefit nonpositive or undefined; angular search not performed";
    else
        available=1:last;searchMax(j)=r.widths(last);event=false(n,last,2);
        for d=available
            raw=load(r.files(d,j),'meta','hit','time');validateRaw(raw,r,d,j);
            event(:,d,:)=reshape(raw.hit & raw.time<=15,n,1,2);
        end
        timely{j}=event;prob=squeeze(mean(event,1));pd(available,:,j)=prob;
        widths=r.widths(available);minimum=sensitivity_min_width(widths,prob,.6);
        bp=reshape(double(reshape(event,n,[])).'*wm/n,last,[]);
        bw=reshape(sensitivity_min_width(widths,bp,.6),2,reps);
        feasibility(j)=mean(all(isfinite(bw),1));absolute(j,:,1)=minimum;
        point(j,1)=100*(1-minimum(2)/minimum(1));boot(:,j,1)=100*(1-bw(2,:)./bw(1,:)).';
        for beam=1:2
            d=find(prob(:,beam)>=.6,1);if ~isempty(d),brackets(j,beam,:)=[widths(max(1,d-1)),widths(d)];end
        end
        status(j)="Complete angular search";
        if any(~isfinite(minimum)),status(j)="Requirement unattained within searched grid";end
    end
    for metric=1:3
        values=boot(:,j,metric);good=isfinite(values);
        if mean(good)>=.975 && isfinite(point(j,metric)),ci(j,metric,:)=quantile(values(good),[.025 .975]);end
    end
end
% Select the largest contiguous tested symmetric interval. No cherry-picked
% interior points; significance is NOT an inclusion rule.
eligible=all(isfinite(point) & point>0,2) & all(isfinite(reshape(ci,nr,[])),2) & feasibility>=.975;
halfRange=zeros(1,4);chosen=false(nr,1);chosen(1)=true;decision=cell(0,1);
for p=1:4
    maximum=round(100*max(abs(r.scale(r.factorIndex==p,p)-1)));
    for h=10:10:maximum
        ids=find(r.factorIndex==0 | (r.factorIndex==p & abs(r.scale(:,p)-1)<=h/100+1e-10));
        expected=2*(h/10)+1;present=numel(ids)==expected;
        pass=present && all(eligible(ids));
        decision{end+1,1}=struct('ParameterIndex',p,'HalfRangePercent',h,'ExpectedCases',expected, ...
            'AllPointsPresent',present,'Pass',pass,'FailedCases',join(string(ids(~eligible(ids))).',',')); %#ok<AGROW>
        if pass,halfRange(p)=h;end
    end
    ids=r.factorIndex==p & abs(r.scale(:,p)-1)<=halfRange(p)/100+1e-10;chosen=chosen|ids;
end
qualifiedHalfRange=halfRange;
% The optional kinematic column must cover at least +/-30 percent to be
% useful in this expanded-range figure; retain smaller domains in the audit.
if halfRange(4)<30
    chosen(r.factorIndex==4)=false;halfRange(4)=0;
end
included=find(chosen);sim=nan(nr,3,2);
if ~isempty(included)
    x=reshape(boot(:,included,:),reps,[]);est=reshape(point(included,:),1,[]);good=all(isfinite(x),2);
    if mean(good)>=.975
        se=std(x(good,:),0,1);critical=quantile(max(abs((x(good,:)-est)./max(se,eps)),[],2),.95);
        sim(included,:,1)=reshape(est-critical*se,[],3);sim(included,:,2)=reshape(est+critical*se,[],3);
    end
end
record=cell(nr*3,1);names=["AngularReduction","BlindReduction","PulseIncrease"];
for j=1:nr
    for metric=1:3
        record{(j-1)*3+metric}=struct('Case',j,'ParameterIndex',r.factorIndex(j), ...
            'EnergyScale',r.scale(j,1),'ReflectivityScale',r.scale(j,2),'AlphaScale',r.scale(j,3),'OmegaScale',r.scale(j,4), ...
            'Metric',names(metric),'LineValue',absolute(j,1,metric),'RingValue',absolute(j,2,metric), ...
            'GainPercent',point(j,metric),'Low',ci(j,metric,1),'High',ci(j,metric,2), ...
            'SimLow',sim(j,metric,1),'SimHigh',sim(j,metric,2),'MainIncluded',chosen(j), ...
            'BootstrapBothFeasible',feasibility(j),'AngularSearchMax',searchMax(j),'Status',status(j));
    end
end
tableOut=struct2table(vertcat(record{:}));writetable(tableOut,fullfile(out,'all_candidates.csv'));
writetable(struct2table(vertcat(decision{:})),fullfile(out,'range_decisions.csv'));
writetable(table((1:nr).',commonGain,counts(:,1),counts(:,2),counts(:,3), ...
    'VariableNames',{'Case','CommonPulseGain','LineSuccessN','RingSuccessN','CommonSuccessN'}),fullfile(out,'pulse_controls.csv'));
theory=struct('TimelyViolations',zeros(1,3),'StaticViolations',zeros(1,3));
for p=1:3
    ids=find(r.factorIndex==p | r.factorIndex==0);[~,order]=sort(r.scale(ids,p));ids=ids(order);
    direction=1;if p==3,direction=-1;end
    theory.StaticViolations(p)=nnz(direction*diff(double(static(:,:,ids)),1,3)<0);
    for q=2:numel(ids)
        a=ids(q-1);b=ids(q);if isempty(timely{a}) || isempty(timely{b}),continue;end
        count=min(size(timely{a},2),size(timely{b},2));
        theory.TimelyViolations(p)=theory.TimelyViolations(p)+nnz(direction*(double(timely{b}(:,1:count,:))-double(timely{a}(:,1:count,:)))<0);
    end
end
analysis=struct('design',r,'point',point,'absolute',absolute,'ci',ci,'sim',sim,'bootstrap',boot, ...
    'feasibility',feasibility,'brackets',brackets,'pd',pd,'commonPulseGain',commonGain, ...
    'counts',counts,'status',status,'halfRange',halfRange,'qualifiedHalfRange',qualifiedHalfRange,'chosen',chosen,'table',tableOut,'theory',theory);
save(fullfile(out,'analysis.mat'),'analysis','-v7.3');
fprintf('SELECTED_HALF_RANGES_PERCENT=%s\n',mat2str(halfRange));disp(theory);
end

function validateRaw(raw,r,d,j)
expected=[r.values(j,1:2),r.values(j,3),r.values(j,3)/r.manifest.alphaBetaRatio];
assert(raw.meta.Case==r.rawCase(d,j) && raw.meta.Width==r.widths(d) ...
    && raw.meta.N==r.bank.options.N && raw.meta.Seed==r.bank.options.Seed ...
    && max(abs(raw.meta.Physics-expected))<1e-15,'Expanded:Metadata','Raw metadata mismatch.');
if isfield(raw.meta,'Omega'),assert(raw.meta.Omega==r.values(j,4),'Expanded:Omega','Incorrect scan speed.');
else,assert(r.scale(j,4)==1,'Expanded:Omega','Nonbaseline speed requires explicit metadata.');end
end

function value=windowScore(w)
variance=w.signal_photons+w.backscatter_photons+w.background_photons+w.dark_count_photons;
valid=w.signal_photons>0 & variance>0 & isfinite(variance);
value=sum(w.signal_photons(valid).^2./variance(valid));
end

function weights=bootstrapWeights(n,reps)
weights=zeros(n,reps);for b=1:reps,weights(:,b)=histcounts(randi(n,n,1),.5:1:n+.5).';end
end




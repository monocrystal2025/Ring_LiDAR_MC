function analysis = analyze_unified_robustness(options)
%ANALYZE_UNIFIED_ROBUSTNESS Paired uncertainty for three parallel advantages.
arguments
    options.Directory (1,1) string = ""
    options.Replicates (1,1) double {mustBeInteger,mustBePositive} = 1000
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 20261001
    options.RequiredProbability (1,1) double {mustBePositive} = .6
end
folder=fileparts(mfilename('fullpath'));
addpath(fileparts(folder));
if options.Directory=="",options.Directory=fullfile(folder,'results');end
saved=load(fullfile(options.Directory,'design.mat'),'result');r=saved.result;o=r.options;
assert(all(r.complete,'all'),'Unified:Incomplete','All configurations must finish before analysis.');
old=rng;cleanup=onCleanup(@()rng(old));rng(options.Seed,'twister');
wm=bootstrapWeights(o.N,options.Replicates);ws=bootstrapWeights(o.StaticN,options.Replicates);
nc=size(r.scale,1);nw=numel(o.Widths);nr=numel(o.ReferenceWidths);reps=options.Replicates;
point=zeros(nc,3,nr);boots=nan(reps,nc,3,nr);absolute=nan(nc,2,3,nr);
counts=zeros(nc,3,nr);commonGain=nan(nc,nr);feasibility=zeros(nc,3);
brackets=nan(nc,2,2);timely=false(o.N,nw,2,nc);static=false(o.StaticN,2,nc,nr);
pulseCounts=nan(o.N,2,nc,nr);allHits=false(o.N,2,nc,nr);
for j=1:nc
    for d=1:nw
        raw=load(fullfile(options.Directory,sprintf('case_%02d_width_%03d.mat',j,o.Widths(d))));
        assert(raw.meta.Case==j && raw.meta.Width==o.Widths(d) && raw.meta.N==o.N ...
            && raw.meta.Seed==o.Seed && isequal(raw.meta.Physics,[r.values(j,1:2),r.values(j,3),r.values(j,3)/r.manifest.alphaBetaRatio]), ...
            'Unified:Raw','Raw metadata inconsistent with saved design.');
        timely(:,d,:,j)=reshape(raw.hit & raw.time<=o.Deadline,o.N,1,2);
        k=find(o.ReferenceWidths==o.Widths(d),1);
        if isempty(k),continue;end
        assert(raw.meta.StaticAvailable && isinf(raw.meta.DynamicLimit),'Unified:Horizon','Reference observations must be complete.');
        static(:,:,j,k)=raw.staticHit;pulseCounts(:,:,j,k)=raw.pulse;allHits(:,:,j,k)=raw.hit;
        assert(all(isfinite(raw.pulse(raw.hit))) && all(isnan(raw.pulse(~raw.hit))), ...
            'Unified:Pulse','Pulse statistics must be conditioned on successful reports.');
        % Recount the saved photon windows; no display scaling is permitted.
        for b=1:2
            wins=raw.windows{b};actual=arrayfun(@(x)nnz(isfinite(x.signal_photons) & x.signal_photons>0),wins);
            assert(isequal(actual(raw.hit(:,b)),raw.pulse(raw.hit(:,b),b)), ...
                'Unified:PulseCount','Stored pulse count differs from the actual successful window.');
            score=arrayfun(@windowScore,wins(raw.hit(:,b)));
            assert(all(score>=4-1e-8),'Unified:SNR','A successful window does not satisfy SNR >= 2.');
        end
        blind=double(~raw.staticHit);bmean=mean(blind);bboot=blind.'*ws/o.StaticN;
        pmean=mean(raw.pulse,1,'omitnan');pv=raw.pulse;pv(~raw.hit)=0;
        pboot=(pv.'*wm)./(double(raw.hit).'*wm);
        common=all(raw.hit,2);commonMeans=mean(raw.pulse(common,:),1);
        commonGain(j,k)=100*(commonMeans(2)/commonMeans(1)-1);
        counts(j,:,k)=[sum(raw.hit),nnz(common)];
        absolute(j,:,2,k)=bmean;absolute(j,:,3,k)=pmean;
        point(j,2,k)=100*(1-bmean(2)/bmean(1));
        point(j,3,k)=100*(pmean(2)/pmean(1)-1);
        boots(:,j,2,k)=100*(1-bboot(2,:)./bboot(1,:)).';
        boots(:,j,3,k)=100*(pboot(2,:)./pboot(1,:)-1).';
    end
    event=timely(:,:,:,j);prob=squeeze(mean(event,1));
    width=sensitivity_min_width(o.Widths,prob,options.RequiredProbability);
    bootProb=reshape(double(reshape(event,o.N,[])).'*wm/o.N,nw,[]);
    bootWidth=reshape(sensitivity_min_width(o.Widths,bootProb,options.RequiredProbability),2,reps);
    for b=1:2
        d=find(prob(:,b)>=options.RequiredProbability,1);
        if ~isempty(d),brackets(j,b,:)=[o.Widths(max(1,d-1)),o.Widths(d)];end
    end
    feasibility(j,:)=[mean(isfinite(bootWidth),2).',mean(all(isfinite(bootWidth),1))];
    for k=1:nr
        absolute(j,:,1,k)=width;
        point(j,1,k)=100*(1-width(2)/width(1));
        boots(:,j,1,k)=100*(1-bootWidth(2,:)./bootWidth(1,:)).';
    end
end
ci=nan(nc,3,2,nr);sim=ci;
for k=1:nr
    valid=all(isfinite(reshape(boots(:,:,:,k),reps,[])),2);
    for j=1:nc
        for metric=1:3
            x=boots(:,j,metric,k);finite=isfinite(x);
            % An interval based on too few attainable replicates is not reported.
            if mean(finite)>=.975 && isfinite(point(j,metric,k))
                ci(j,metric,:,k)=quantile(x(finite),[.025 .975]);
            end
        end
    end
    if mean(valid)>=.975
        x=reshape(boots(valid,:,:,k),nnz(valid),[]);est=reshape(point(:,:,k),1,[]);
        se=std(x,0,1);critical=quantile(max(abs((x-est)./max(se,eps)),[],2),.95);
        sim(:,:,1,k)=reshape(est-critical*se,nc,3);sim(:,:,2,k)=reshape(est+critical*se,nc,3);
    end
end
records=cell(nc*nr*3,1);index=0;
names=["AngularReduction","BlindReduction","PulseIncrease"];
for k=1:nr
    for j=1:nc
        for metric=1:3
            index=index+1;records{index}=struct('Case',j,'FactorIndex',r.factorIndex(j), ...
                'Energy_uJ',r.values(j,1)*1e6,'Reflectivity',r.values(j,2),'Alpha',r.values(j,3), ...
                'ReferenceWidth',o.ReferenceWidths(k),'Metric',names(metric), ...
                'LineValue',absolute(j,1,metric,k),'RingValue',absolute(j,2,metric,k), ...
                'GainPercent',point(j,metric,k),'Low',ci(j,metric,1,k),'High',ci(j,metric,2,k), ...
                'SimLow',sim(j,metric,1,k),'SimHigh',sim(j,metric,2,k), ...
                'WidthBootstrapBothFeasible',feasibility(j,3), ...
                'LineWidthBracketLow',brackets(j,1,1),'LineWidthBracketHigh',brackets(j,1,2), ...
                'RingWidthBracketLow',brackets(j,2,1),'RingWidthBracketHigh',brackets(j,2,2));
        end
    end
end
tab=struct2table(vertcat(records{:}));writetable(tab,fullfile(options.Directory,'three_benefits.csv'));
commonTable=table();
for k=1:nr
    one=table((1:nc).',repmat(o.ReferenceWidths(k),nc,1),counts(:,1,k),counts(:,2,k), ...
        counts(:,3,k),point(:,3,k),commonGain(:,k),'VariableNames', ...
        {'Case','ReferenceWidth','LineSuccessfulN','RingSuccessfulN','CommonSuccessfulN', ...
        'SeparateSuccessPulseGain','CommonSuccessPulseGain'});
    commonTable=[commonTable;one]; %#ok<AGROW>
end
writetable(commonTable,fullfile(options.Directory,'common_pulse_comparison.csv'));
% Model-implied directions are checked on absolute detection/coverage data.
theory=struct('TimelyViolations',zeros(1,3),'StaticViolations',zeros(1,3));
for p=1:3
    ids=find(r.factorIndex==p | r.factorIndex==0);[~,ord]=sort(r.scale(ids,p));ids=ids(ord);
    direction=1;if p==3,direction=-1;end
    theory.TimelyViolations(p)=nnz(direction*diff(double(timely(:,:,:,ids)),1,4)<0);
    theory.StaticViolations(p)=nnz(direction*diff(double(static(:,:,ids,:)),1,3)<0);
end
analysis=struct('options',options,'design',r,'point',point,'ci',ci,'sim',sim, ...
    'absolute',absolute,'counts',counts,'commonPulseGain',commonGain,'feasibility',feasibility, ...
    'brackets',brackets,'table',tab,'theory',theory,'timely',timely,'staticHit',static, ...
    'pulseCounts',pulseCounts,'allHits',allHits);
save(fullfile(options.Directory,'analysis.mat'),'analysis','-v7.3');
fid=fopen(fullfile(options.Directory,'physical_checks.json'),'w');fprintf(fid,'%s',jsonencode(theory,PrettyPrint=true));fclose(fid);
disp(tab(tab.ReferenceWidth==125,:));
end

function score=windowScore(window)
signal=window.signal_photons;
variance=signal+window.backscatter_photons+window.background_photons+window.dark_count_photons;
valid=signal>0 & variance>0 & isfinite(variance);
score=sum(signal(valid).^2./variance(valid));
end

function weights=bootstrapWeights(n,reps)
weights=zeros(n,reps);
for b=1:reps,weights(:,b)=histcounts(randi(n,n,1),.5:1:n+.5).';end
end

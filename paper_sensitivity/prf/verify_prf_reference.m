function report=verify_prf_reference()
%VERIFY_PRF_REFERENCE Audit reference data and simultaneous uncertainty.
folder=fileparts(mfilename('fullpath'));s=load(fullfile(folder,'results','analysis.mat'),'a');a=s.a;r=a.design;
previous=rng;cleanup=onCleanup(@()rng(previous));rng(20260928,'twister');
n=r.bank.options.N;ns=r.bank.options.StaticN;reps=size(a.bootstrap,1);
% Recreate the same static weights, after dynamic-weight RNG consumption.
for k=1:reps,randi(n,n,1);end
ws=zeros(ns,reps);for k=1:reps,ws(:,k)=histcounts(randi(ns,ns,1),.5:1:ns+.5).';end
bd=zeros(reps,10);ref=find(r.widths==125);raws=cell(10,1);
for j=1:10
    raw=load(r.files(ref,j));raws{j}=raw;
    dv=100*(double(~raw.staticHit(:,1))-double(~raw.staticHit(:,2)));
    bd(:,j)=(dv.'*ws/ns).';
    assert(isequal(size(raw.staticHit),[ns 2]));
    for b=1:2
        windows=raw.windows{b};scores=zeros(n,1);
        for i=find(raw.hit(:,b)).'
            w=windows(i);v=w.signal_photons+w.backscatter_photons+w.background_photons+w.dark_count_photons;
            good=w.signal_photons>0 & isfinite(v) & v>0;scores(i)=sum(w.signal_photons(good).^2./v(good));
        end
        assert(all(scores(raw.hit(:,b))>=4-1e-8));
    end
end
x=[bd,a.bootstrap(:,:,3)];est=[a.blindDifference(:,1).',a.point(:,3).'];
assert(all(isfinite(x),'all'));se=std(x,0,1);
critical=quantile(max(abs((x-est)./max(se,eps)),[],2),.95);
lo=est-critical*se;hi=est+critical*se;
T=table(r.frequencies.'/1000,lo(1:10).',hi(1:10).',lo(11:20).',hi(11:20).', ...
    'VariableNames',{'FrequencyKHz','BlindDifferenceLow','BlindDifferenceHigh','PulseLow','PulseHigh'});
writetable(T,fullfile(folder,'results','reference_simultaneous.csv'));disp(T);
old=load(fullfile(fileparts(folder),'expanded','results','case_40_width_125.mat'),'staticHit');
assert(isequal(raws{10}.staticHit,old.staticHit),'PRF:Ratio','Static f/omega identity failed.');
assert(max(a.scanErrors)<.001,'PRF:Scan','Incorrect scan angular speed.');
assert(max(a.table.CycleSeconds)-min(a.table.CycleSeconds)<.005,'PRF:Period','Period changed materially.');
report=struct('ReferenceChecksPassed',true,'FrequenciesKHz',r.frequencies/1000,'MovingN',n,'StaticN',ns, ...
    'EnergyMicrojoules',r.manifest.baseline(1)*1e6,'Omega',r.manifest.omega, ...
    'MaxRelativeScanSpeedError',max(a.scanErrors),'StaticRatioIdentity',true, ...
    'SimultaneousScope','20 reference endpoints: 10 absolute blind differences and 10 pulse gains', ...
    'SimultaneousCritical',critical,'PulseReversalKHz',r.frequencies(T.PulseHigh<0)/1000, ...
    'BlindReversalKHz',r.frequencies(T.BlindDifferenceHigh<0)/1000);
fid=fopen(fullfile(folder,'results','reference_verification.json'),'w');fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));fclose(fid);
disp(report);
end

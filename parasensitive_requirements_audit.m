function audit = parasensitive_requirements_audit()
%PARASENSITIVE_REQUIREMENTS_AUDIT Compare requirements without new simulation.
root = fileparts(mfilename('fullpath'));
s = load(fullfile(root,'parasensitive_data.mat'),'data');
d = s.data;
times = 5:5:60;
probabilities = 0.30:0.05:0.95;
requirements = [15 .6;30 .7;60 .9;19 .88;20 .88];
[tt,pp] = ndgrid(times,probabilities);
requirements = unique([requirements;tt(:) pp(:)],'rows');
reductions = nan(size(requirements,1),20);
minima = nan(size(requirements,1),20,2);
for r = 1:size(requirements,1)
    for c = 1:20
        j = d.caseMap(c);
        for b = 1:2
            probability = zeros(numel(d.config.WidthsMrad),1);
            for k = 1:numel(probability)
                raw = d.dynamic{k,j,b};
                t = requirements(r,1);
                probability(k) = mean(raw.detect~=0 & isfinite(raw.firstTime) ...
                    & raw.firstTime>=0 & raw.firstTime<t+eps(t));
            end
            minima(r,c,b) = minimumWidth(d.config.WidthsMrad,probability,requirements(r,2));
        end
    end
    reductions(r,:) = 100*(1-minima(r,:,2)./minima(r,:,1));
end
summary = table(requirements(:,1),requirements(:,2),sum(isfinite(reductions),2), ...
    reductions(:,3),median(reductions,2,'omitnan'),sum(reductions>0,2), ...
    min(reductions,[],2,'omitnan'),max(reductions,[],2,'omitnan'), ...
    'VariableNames',{'Time_s','PD','ValidOf20','BaselinePct','MedianPct', ...
    'PositiveOf20','MinPct','MaxPct'});
audit = struct('requirements',requirements,'reductionsPct',reductions, ...
    'minimumWidthsMrad',minima,'summary',summary,'N',d.config.N);
writetable(summary,fullfile(root,'parasensitive_requirements_audit.csv'));
save(fullfile(root,'parasensitive_requirements_audit.mat'),'audit');
fprintf('Original requirements:\n');
disp(summary(ismember(requirements,[15 .6;30 .7;60 .9],'rows'),:));
full = summary(summary.ValidOf20==20,:);
full = sortrows(full,{'MedianPct','BaselinePct'},{'descend','descend'});
fprintf('Largest median reductions among fully covered requirements (exploratory):\n');
disp(full(1:min(12,height(full)),:));
fprintf('Selected alternatives:\n');
selected = [19 .88;20 .88;20 .85;30 .85;40 .9;50 .9;60 .8;60 .85;60 .95];
for k = 1:size(selected,1)
    index = find(all(abs(requirements-selected(k,:))<1e-10,2),1);
    if ~isempty(index)
        disp(summary(index,:));
    end
end
fprintf('Baseline raw cumulative detection probabilities:\n');
j = d.caseMap(3,1);
out = zeros(numel(d.config.WidthsMrad),7);
out(:,1) = d.config.WidthsMrad(:);
originalTimes = [15 30 60];
for k = 1:size(out,1)
    for b = 1:2
        raw = d.dynamic{k,j,b};
        for r = 1:3
            out(k,1+(r-1)*2+b) = mean(raw.detect & raw.firstTime<=originalTimes(r));
        end
    end
end
disp(out);
paper = load(fullfile(root,'minimum_thetaD.mat'));
comparisonReq = [15 .6;30 .7;60 .9;19 .88];
comparison = nan(4,7);
for r = 1:4
    t = find(paper.timeThresholds_s==comparisonReq(r,1));
    [~,q] = min(abs(paper.targetProbabilities-comparisonReq(r,2)));
    idx = find(all(abs(requirements-comparisonReq(r,:))<1e-10,2),1);
    coarseMin = nan(1,2);
    for b = 1:2
        series = paper.results(1).Series(b+1);
        [exists,indices] = ismember(d.config.WidthsMrad,series.WD_mrad);
        assert(all(exists));
        coarseMin(b) = minimumWidth(d.config.WidthsMrad,series.PD(t,indices),comparisonReq(r,2));
    end
    comparison(r,:) = [comparisonReq(r,:),paper.results(1).LineMinWD_mrad(q,t), ...
        paper.results(1).RingMinWD_mrad(q,t),100*(1-paper.results(1).Kappa(q,t)), ...
        100*(1-coarseMin(2)/coarseMin(1)),reductions(idx,3)];
end
comparison = array2table(comparison,'VariableNames',{'Time_s','PD','PaperLineMrad', ...
    'PaperRingMrad','PaperReductionPct','PaperCoarseGridPct','CurrentBaselinePct'});
disp(comparison);
audit.paperComparison = comparison;
save(fullfile(root,'parasensitive_requirements_audit.mat'),'audit');
writetable(comparison,fullfile(root,'parasensitive_paper_comparison.csv'));
end

function width = minimumWidth(w,p,target)
p = cummax(p(:));
k = find(p>=target,1);
width = nan;
if ~isempty(k)
    if k==1
        width = w(1);
    elseif abs(diff(p(k-1:k)))<=eps(max(abs(p(k-1:k))))
        width = w(k);
    else
        width = interp1(p(k-1:k),w(k-1:k),target);
    end
end
end

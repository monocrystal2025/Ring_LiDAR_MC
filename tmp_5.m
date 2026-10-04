function result = tmp_5(mode,options)
%TMP_5 Operating-region frontier with independent design-selection validation.
% Designs are selected using training trajectories only, then evaluated
% using an independent validation ensemble. No forced monotonic response.
arguments
    mode (1,1) string = "preview"
    options.N (1,1) double {mustBeInteger,mustBeNonnegative} = 0
end
cfg = discussion.config(mode);
if options.N>0
    cfg.N = options.N;
end
cfg.outputRoot = fullfile(cfg.outputRoot,"tmp_5");
trainingBank = discussion.bank(cfg.N,cfg.seed+100);
validationBank = discussion.bank(2*cfg.N,cfg.seed+200);
training = discussion.runCurve(cfg,trainingBank,"training");
validation = discussion.runCurve(cfg,validationBank,"validation");
timeGrid = 3:2:15;
probabilityGrid = 0.5:0.05:0.95;
kappa = nan(numel(probabilityGrid),numel(timeGrid));
leftBoundary = false(size(kappa));
certified = false(size(kappa));
frontiers = cell(size(kappa));
for a = 1:numel(probabilityGrid)
    for b = 1:numel(timeGrid)
        f = discussion.frontier(training,timeGrid(b),probabilityGrid(a));
        frontiers{a,b} = f;
        kappa(a,b) = f.kappa;
        leftBoundary(a,b) = f.lineStatus=="left-boundary" || f.ringStatus=="left-boundary";
        certified(a,b) = f.certifiedGridAdvantage;
    end
end
fieldBudgets = [45,60,75,90,105,125]; % Conservative half field in mrad.
trainHit = false(cfg.N,numel(fieldBudgets),2);
validationHit = false(validationBank.N,numel(fieldBudgets),2);
selectedWidths = nan(numel(fieldBudgets),2);
validBudget = false(numel(fieldBudgets),1);
rows = cell(numel(fieldBudgets),1);
for a = 1:numel(fieldBudgets)
    feasible = find(cfg.widthGrid_mrad/2+3*cfg.system.beamWidth_rad*1e3 <= fieldBudgets(a));
    if isempty(feasible)
        rows{a} = struct('HalfFieldBudget_mrad',fieldBudgets(a),'LineSelectedWidth_mrad',NaN, ...
            'RingSelectedWidth_mrad',NaN,'LineValidationP',NaN,'RingValidationP',NaN, ...
            'DeltaP',NaN,'Low95',NaN,'High95',NaN,'ValidationN',validationBank.N);
        continue
    end
    validBudget(a) = true;
    for k = 1:2
        p = mean(training.times(:,feasible,k)<=cfg.maxTime_s,1);
        % Among tied maxima choose the smallest width, using training only.
        [~,best] = max(p);
        index = feasible(best);
        selectedWidths(a,k) = cfg.widthGrid_mrad(index);
        trainHit(:,a,k) = training.times(:,index,k)<=cfg.maxTime_s;
        validationHit(:,a,k) = validation.times(:,index,k)<=cfg.maxTime_s;
    end
    s = discussion.paired(double(validationHit(:,a,1)),double(validationHit(:,a,2)),cfg);
    rows{a} = struct('HalfFieldBudget_mrad',fieldBudgets(a), ...
        'LineSelectedWidth_mrad',selectedWidths(a,1),'RingSelectedWidth_mrad',selectedWidths(a,2), ...
        'LineValidationP',s.lineMean,'RingValidationP',s.ringMean, ...
        'DeltaP',s.delta,'Low95',s.low,'High95',s.high,'ValidationN',validationBank.N);
end
summary = struct2table(vertcat(rows{:}));
[fig,layout,colors] = discussion.figure("Validated angular-field tradeoff",2,2);
ax = nexttile(layout);
h = imagesc(ax,timeGrid,probabilityGrid,kappa);
set(h,'AlphaData',isfinite(kappa));
set(ax,'Color',[0.85,0.85,0.85]); axis(ax,'xy');
colorbar(ax);
hold(ax,'on');
[aa,bb] = find(leftBoundary & isfinite(kappa));
plot(ax,timeGrid(bb),probabilityGrid(aa),'kx','MarkerSize',5);
[aa,bb] = find(certified & ~leftBoundary);
plot(ax,timeGrid(bb),probabilityGrid(aa),'k.','MarkerSize',8);
discussion.styleAxes(ax,"Warning deadline T (s)","Required P_T", ...
    "(a) Training \kappa (x: grid boundary)");
ax = nexttile(layout);
handles = discussion.probabilityPlot(ax,cfg.widthGrid_mrad, ...
    training.times<=cfg.maxTime_s,colors);
hold(ax,'on');
for k = 1:2
    plot(ax,cfg.widthGrid_mrad,100*mean(validation.times(:,:,k)<=cfg.maxTime_s,1), ...
        ':','Color',colors(k,:),'LineWidth',1.4,'HandleVisibility','off');
end
legend(ax,handles,["Line","Annular"],'Box','off','Location','southeast');
discussion.styleAxes(ax,"w_D (mrad)","P_T (%)","(b) Training and validation");
ax = nexttile(layout);
discussion.probabilityPlot(ax,fieldBudgets(validBudget),validationHit(:,validBudget,:),colors);
discussion.styleAxes(ax,"Half-field budget (mrad)","Validation P_T (%)", ...
    "(c) Selected designs: validation");
ax = nexttile(layout);
hold(ax,'on');
idx = find(validBudget);
errorbar(ax,fieldBudgets(idx),100*summary.DeltaP(idx), ...
    100*(summary.DeltaP(idx)-summary.Low95(idx)), ...
    100*(summary.High95(idx)-summary.DeltaP(idx)), ...
    '-o','Color',colors(2,:),'LineWidth',1.2,'MarkerSize',4);
yline(ax,0,':','HandleVisibility','off');
discussion.styleAxes(ax,"Half-field budget (mrad)","\Delta P_T (pp)", ...
    "(d) Paired validation difference");
result = struct('training',training,'validation',validation, ...
    'timeGrid_s',timeGrid,'probabilityGrid',probabilityGrid, ...
    'frontiers',{frontiers},'kappa',kappa,'leftBoundary',leftBoundary, ...
    'certifiedGridAdvantage',certified,'selectedWidths_mrad',selectedWidths, ...
    'fieldBudgets_mrad',fieldBudgets,'trainHit',trainHit, ...
    'validationHit',validationHit,'summary',summary);
discussion.finish(cfg,"tmp_5",result,fig,summary);
end

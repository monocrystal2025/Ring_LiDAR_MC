function result = tmp_3(mode,options)
%TMP_3 Mechanism decomposition and matched-trajectory accumulation ablations.
% All indicators include misses. Removing accumulation changes only the
% decision rule on identical photon sequences; it does not change scanning.
arguments
    mode (1,1) string = "preview"
    options.N (1,1) double {mustBeInteger,mustBeNonnegative} = 0
end
cfg = discussion.config(mode);
if options.N>0
    cfg.N = options.N;
end
cfg.outputRoot = fullfile(cfg.outputRoot,"tmp_3");
bank = discussion.bank(cfg.N,cfg.seed);
widths = [90,120,150,195];
raw = cell(numel(widths),2);
curves = struct('widths_mrad',widths,'times',inf(cfg.N,numel(widths),2), ...
    'singleTimes',inf(cfg.N,numel(widths),2),'uniformTimes',inf(cfg.N,numel(widths),2), ...
    'N',cfg.N);
rows = cell(numel(widths),1);
decomposition = zeros(numel(widths),2,2);
for j = 1:numel(widths)
    for k = 1:2
        names = ["line","ring"];
        raw{j,k} = discussion.runCase(cfg,discussion.beam(names(k),widths(j)),bank,"mechanism");
        m = raw{j,k}.metrics;
        curves.times(:,j,k) = m.firstTime;
        curves.singleTimes(:,j,k) = m.singleTime;
        curves.uniformTimes(:,j,k) = m.uniformTime;
        encountered = isfinite(m.encounterTime);
        detected = isfinite(m.firstTime);
        decomposition(j,k,1) = mean(encountered);
        decomposition(j,k,2) = nnz(detected)/max(1,nnz(encountered));
        assert(all(~detected | encountered),'discussion:Decomposition','Detection must imply encounter.');
    end
    lo = raw{j,1}.metrics; ri = raw{j,2}.metrics;
    weighted = discussion.paired(double(isfinite(lo.firstTime)),double(isfinite(ri.firstTime)),cfg);
    single = discussion.paired(double(isfinite(lo.singleTime)),double(isfinite(ri.singleTime)),cfg);
    % Difference in differences: how much more the annular comparison
    % benefits from accumulation; compute it on each paired trial.
    gainL = double(isfinite(lo.firstTime))-double(isfinite(lo.singleTime));
    gainR = double(isfinite(ri.firstTime))-double(isfinite(ri.singleTime));
    interaction = discussion.paired(gainL,gainR,cfg);
    rows{j} = struct('Width_mrad',widths(j),'WeightedDeltaP',weighted.delta, ...
        'WeightedLow95',weighted.low,'WeightedHigh95',weighted.high, ...
        'SingleDeltaP',single.delta,'AccumulationInteraction',interaction.delta, ...
        'InteractionLow95',interaction.low,'InteractionHigh95',interaction.high, ...
        'MeanParticipationLine',mean(lo.participation), ...
        'MeanParticipationRing',mean(ri.participation), ...
        'MeanSignalLine',mean(lo.signalPhotons),'MeanSignalRing',mean(ri.signalPhotons), ...
        'N',cfg.N);
end
summary = struct2table(vertcat(rows{:}));
% Closure ablation: remove half the annulus. Loss and equal-total-energy
% versions distinguish closure from simply throwing away half the photons.
ablationLabels = ["Full annulus","Half annulus: loss","Half annulus: equal energy"];
ablation = cell(1,3);
for j = 1:3
    caseCfg = cfg;
    if j>1
        caseCfg.ringGap_deg = 180;
        caseCfg.gapRedistribute = j==3;
    end
    ablation{j} = discussion.runCase(caseCfg,discussion.beam("ring",120), ...
        bank,"closure_"+j);
end
[fig,layout,colors] = discussion.figure("Mechanism and accumulation",2,2);
ax = nexttile(layout);
handles = discussion.probabilityPlot(ax,widths,isfinite(curves.times),colors);
hold(ax,'on');
for k = 1:2
    plot(ax,widths,100*mean(isfinite(curves.singleTimes(:,:,k)),1),':', ...
        'Color',colors(k,:),'LineWidth',1.5,'HandleVisibility','off');
end
legend(ax,handles,["Line","Annular"],'Box','off','Location','southeast');
discussion.styleAxes(ax,"w_D (mrad)","P_T (%)","(a) Weighted vs. single-pulse");
ax = nexttile(layout);
errorbar(ax,widths,100*summary.AccumulationInteraction, ...
    100*(summary.AccumulationInteraction-summary.InteractionLow95), ...
    100*(summary.InteractionHigh95-summary.AccumulationInteraction), ...
    '-o','Color',colors(2,:),'LineWidth',1.3,'MarkerSize',4);
yline(ax,0,':');
discussion.styleAxes(ax,"w_D (mrad)","Accumulation interaction (pp)", ...
    "(b) Accumulation interaction");
ax = nexttile(layout);
hold(ax,'on');
for k = 1:2
    plot(ax,widths,100*decomposition(:,k,1),'-o','Color',colors(k,:), ...
        'MarkerSize',3,'LineWidth',1.2);
    plot(ax,widths,100*decomposition(:,k,2),'--s','Color',colors(k,:), ...
        'MarkerSize',3,'LineWidth',1.2);
end
ylim(ax,[0,100]);
discussion.styleAxes(ax,"w_D (mrad)","Probability (%)", ...
    "(c) Encounter (solid), conversion (dash)");
ax = nexttile(layout);
hold(ax,'on');
index = 2; % Predeclared 120 mrad example, not chosen for favorable outcomes.
for k = 1:2
    q = raw{index,k}.metrics.participation;
    q = sort(q);
    plot(ax,q,(1:numel(q))/numel(q),'Color',colors(k,:), ...
        'LineWidth',1.3);
end
discussion.styleAxes(ax,"Q participation ratio","Empirical CDF", ...
    "(d) Participation, all trials");
[fig2,layout2,~] = discussion.figure("Closure ablation and window dependence",1,2);
set(fig2,'Position',[2,2,18,8]);
ax = nexttile(layout2);
p = cellfun(@(x)mean(isfinite(x.metrics.firstTime)),ablation);
counts = p*cfg.N;
[lower,upper] = discussion.wilson(counts,cfg.N);
bar(ax,1:3,p*100,'FaceColor',colors(2,:),'FaceAlpha',0.7);
hold(ax,'on');
errorbar(ax,1:3,p*100,(p-lower)*100,(upper-p)*100,'k.','LineWidth',0.8);
ylim(ax,[0,100]); xticks(ax,1:3);
xticklabels(ax,["Full ring","Half, lossy","Half, fixed E"]);
xtickangle(ax,20);
discussion.styleAxes(ax,"Annular geometry","P_T (%)","(a) Closure and energy controls");
ax = nexttile(layout2);
rules = ["Single pulse","Uniform window","Optimal weights"];
pRules = zeros(3,2);
for k = 1:2
    m = raw{2,k}.metrics;
    pRules(:,k) = [mean(isfinite(m.singleTime)); ...
        mean(isfinite(m.uniformTime));mean(isfinite(m.firstTime))];
end
bar(ax,pRules*100);
xticks(ax,1:3); xticklabels(ax,rules); ylim(ax,[0,100]);
lg = legend(ax,["Line","Annular"],'Box','off','Orientation','horizontal');
lg.Layout.Tile = 'south';
lg.ItemTokenSize = [12,8];
discussion.styleAxes(ax,"Decision rule","P_T (%)","(b) Decision-rule ablation");
result = struct('summary',summary,'raw',{raw},'curves',curves, ...
    'decomposition',decomposition,'closureLabels',ablationLabels, ...
    'closureRaw',{ablation},'rules',rules,'ruleProbability',pRules);
discussion.finish(cfg,"tmp_3",result,[fig,fig2],summary);
end

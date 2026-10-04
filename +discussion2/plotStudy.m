function figures = plotStudy(result)
%PLOTSTUDY Show all three claims and physical sensitivity at equal ability.
arguments
    result (1,1) struct
end
cfg = result.config;
w = cfg.widthGrid_mrad;
colors = [0,0.447,0.698;0.835,0.369,0];
styles = [":","-","--"];
beamMarkers = ["s","o"];
energyLabels = ["0.5 E_0","E_0","2 E_0"];
[fig1,layout] = discussion.figure("Three manuscript advantages",2,2);
fig1.Position = [2,2,22,15];
legendHandles = gobjects(6,1);
legendLabels = strings(6,1);
for panel = 1:4
    ax = nexttile(layout);
    hold(ax,'on');
    for j = 1:3
        curve = result.curves{result.exampleIndices(j)};
        for k = 1:2
            if panel==2 && j~=2, continue; end
            if panel == 1
                values = 100*curve.probability(:,k);
            elseif panel == 2
                values = 100*result.geometricCoverage.curves{1}.blind(:,k);
            elseif panel == 3
                values = curve.pulseMean(:,k);
            else
                values = curve.clusterMean(:,k);
            end
            h = plot(ax,w,values,'LineStyle',styles(j),'Color',colors(k,:), ...
                'LineWidth',1.3,'Marker',beamMarkers(k),'MarkerSize',3);
            if panel == 1
                legendHandles(k+2*(j-1)) = h;
                beamLabels = ["Line","Annular"];
                legendLabels(k+2*(j-1)) = beamLabels(k)+", "+energyLabels(j);
            end
        end
    end
    if panel == 1
        yline(ax,100*cfg.targetProbability,':','80%','HandleVisibility','off');
        iso = result.minima{1};
        for k = 1:2
            if isfinite(iso.widths_mrad(k))
                xline(ax,iso.widths_mrad(k),':','Color',colors(k,:),'HandleVisibility','off');
            end
        end
        ylabelText = "Overall P_D (%)";
        titleText = "(a) Same detection requirement";
        ylim(ax,[0,100]);
    elseif panel == 2
        ylabelText = "Geometric blind fraction (%)";
        titleText = "(b) Geometric blind zone; one full scan";
        ymax = max(5,5*ceil(max(result.geometricCoverage.curves{1}.blind,[],'all')*20)+5);
        ylim(ax,[0,ymax]);
    elseif panel == 3
        ylabelText = "Mean N_{eff} | detection";
        titleText = "(c) Pulses in first detection window";
    else
        ylabelText = "Mean C_{eff} | detection";
        titleText = "(d) Clusters in first detection window";
        yline(ax,1,':','HandleVisibility','off');
    end
    discussion.styleAxes(ax,"w_D (mrad)",ylabelText,titleText);
    xlim(ax,[min(w),max(w)]);
end
lg = legend(legendHandles,legendLabels,'NumColumns',3,'Box','off','FontSize',8);
lg.Layout.Tile = 'south';
title(layout,sprintf('ROE: N=%d moving targets; N=%d static probes',cfg.N,cfg.staticN), ...
    'FontWeight','normal','FontSize',10);

% Every sensitivity metric uses each beam's FIRST TESTED width reaching the
% SAME 80% overall Pd, rather than comparing at a fixed common width.
[fig2,layout2] = discussion.figure("Sensitivity of the three advantages",2,2);
fig2.Position = [2,2,25,19];
s = result.summary;
xValues = [s.Kappa,100*(s.LineGeometricBlind-s.RingGeometricBlind), ...
    s.RingPulses./s.LinePulses,s.RingClusters-s.LineClusters];
labels = ["\kappa = w_{D,R}^* / w_{D,L}^*", ...
    "Blind reduction: B_L - B_R (pp)", ...
    "Pulse ratio: N_{eff,R} / N_{eff,L}", ...
    "Cluster gain: C_{eff,R} - C_{eff,L}"];
titles = ["(a) Smaller required angular width","(b) Geometric blind zone at selected widths", ...
    "(c) Pulse use at equal required ability","(d) Cluster count at equal required ability"];
names = string({result.cases(2:2:20).name});
markers = ["o","s"];
for panel = 1:4
    ax = nexttile(layout2);
    hold(ax,'on');
    available = xValues(:,panel);
    finiteValues = available(isfinite(available));
    reference = 0;
    if ismember(panel,[1,3])
        reference = 1;
    end
    bounds = [min([finiteValues;reference]),max([finiteValues;reference])];
    padding = max(diff(bounds)*0.12,0.12);
    missingX = bounds(2)+padding;
    for j = 2:height(s)
        y = s.Factor(j)+(s.Level(j)-1.5)*0.23;
        k = s.Level(j);
        if isfinite(available(j))
            plot(ax,available(j),y,markers(k),'Color',colors(k,:), ...
                'MarkerFaceColor',colors(k,:),'MarkerSize',4);
        else
            plot(ax,missingX,y,'x','Color',colors(k,:),'MarkerSize',4);
        end
    end
    xline(ax,reference,':','HandleVisibility','off');
    if isfinite(available(1))
        xline(ax,available(1),'--','Color',[0.5,0.5,0.5],'HandleVisibility','off');
    end
    xlim(ax,[bounds(1)-padding,bounds(2)+2*padding]);
    ylim(ax,[0.4,10.6]);
    yticks(ax,1:10);
    yticklabels(ax,names);
    ax.YDir = 'reverse';
    discussion.styleAxes(ax,labels(panel),"",titles(panel));
    text(ax,0.98,0.02,'o low; square high; x unavailable', ...
        'Units','normalized','HorizontalAlignment','right','FontSize',7);
end
title(layout2,'P_D^* = 80%; discrete-width point estimates; grey dash = nominal', ...
    'FontWeight','normal','FontSize',10);

fig3 = clusterFigure(result,colors);
fig4 = requirementFigure(result,colors);
fig5 = blindFigure(result,colors);
figures = [fig1,fig2,fig3,fig4,fig5];
end

function fig = clusterFigure(result,colors)
[fig,layout] = discussion.figure("Pulse-cluster evidence and timing baseline",2,2);
fig.Position = [2,2,22,15];
iso = result.minima{1};
iso.widths_mrad = [195,120];
iso.indices = nan(1,2);
for k = 1:2
    found = find(result.config.widthGrid_mrad==iso.widths_mrad(k),1);
    if ~isempty(found), iso.indices(k) = found; end
end
ax1 = nexttile(layout,1);
ax2 = nexttile(layout,2);
ax3 = nexttile(layout,3);
ax4 = nexttile(layout,4);
hold(ax1,'on');
hold(ax3,'on');
if any(~isfinite(iso.indices))
    text(ax1,0.5,0.5,'Nominal iso-performance pair unavailable','HorizontalAlignment','center');
    return
end
curves = result.curves{1};
bars = zeros(3,2);
barN = zeros(1,2);
raw = cell(1,2);
for k = 1:2
    raw{k} = curves.raw{iso.indices(k),k};
    m = raw{k}.metrics;
    hit = isfinite(m.firstTime);
    barN(k) = nnz(hit);
    c = m.clusters(hit);
    thresholds = 1:max(4,max(c));
    plot(ax1,thresholds,arrayfun(@(x)mean(c>=x),thresholds)*100, ...
        '-o','Color',colors(k,:),'LineWidth',1.3,'MarkerSize',4);
    bars(:,k) = [mean(c>=2);mean(m.qualifiedClusters(hit)>=2); ...
        mean(hit & m.qualifiedClusters>=2)]*100;
    span = sort(m.clusterSpan_ms(isfinite(m.clusterSpan_ms)));
    if ~isempty(span)
        stairs(ax3,span,(1:numel(span))/numel(span),'Color',colors(k,:),'LineWidth',1.3);
    end
end
lg = legend(ax1,["Line, "+iso.widths_mrad(1)+" mrad", ...
    "Annular, "+iso.widths_mrad(2)+" mrad"],'Box','off','FontSize',8);
lg.Location = 'northeast';
ylim(ax1,[0,100]);
discussion.styleAxes(ax1,"Cluster threshold c","P(C_{eff} >= c | detection) (%)", ...
    "(a) First-window cluster distribution");
bar(ax2,bars);
xticks(ax2,1:3);
xticklabels(ax2,["C >= 2 | det.","C_Q >= 2 | det.","Det. and C_Q >= 2"]);
xtickangle(ax2,15);
ylim(ax2,[0,100]);
discussion.styleAxes(ax2,"","Probability (%)","(b) Geometric and Q-qualified clusters");
discussion.styleAxes(ax3,"Qualified-cluster timing span (ms)","Conditional CDF", ...
    "(c) Temporal baseline; C_Q >= 2");
ylim(ax3,[0,1]);
% Illustrative annular trial: nearest the median qualified timing span.
% It is explicitly not selected for maximum count or maximum signal.
m = raw{2}.metrics;
candidates = find(isfinite(m.clusterSpan_ms));
if ~isempty(candidates)
    [~,localIndex] = min(abs(m.clusterSpan_ms(candidates)-median(m.clusterSpan_ms(candidates))));
    example = candidates(localIndex);
    window = raw{2}.windows{example};
    x = (window.steps-window.steps(end))/result.config.system.prf_Hz*1000;
    stem(ax4,x,window.q,'filled','Color',colors(2,:),'MarkerSize',3,'LineWidth',1);
    titleText = sprintf('(d) Annular example #%d: C=%d, C_Q=%d', ...
        example,m.clusters(example),m.qualifiedClusters(example));
else
    text(ax4,0.5,0.5,'No qualified multi-cluster detections', ...
        'HorizontalAlignment','center','Units','normalized');
    titleText = "(d) No qualified timing example";
end
discussion.styleAxes(ax4,"Time relative to first detection (ms)","q_i = S_i^2/(S_i+B_i)",titleText);
title(layout,sprintf('Manuscript pair L=195, R=120 mrad; detected counts L=%d, R=%d',barN), ...
    'FontWeight','normal','FontSize',10);
end

function fig = requirementFigure(result,colors)
% All six requirements are shown, including the original 80% null result.
[fig,layout] = discussion.figure("Dependence on the common detection requirement",1,3);
fig.Position = [2,2,26,9];
energy = [0.5,1,2];
maximumWidth = max(result.config.widthGrid_mrad);
for a = 1:3
    ax = nexttile(layout);
    hold(ax,'on');
    index = result.exampleIndices(a);
    widths = nan(numel(result.requirementLevels),2);
    for b = 1:numel(result.requirementLevels)
        widths(b,:) = result.requirements{index,b}.widths_mrad;
    end
    h = gobjects(1,2);
    for k = 1:2
        h(k) = plot(ax,100*result.requirementLevels,widths(:,k),'-o', ...
            'Color',colors(k,:),'LineWidth',1.3,'MarkerSize',4);
        missing = ~isfinite(widths(:,k));
        plot(ax,100*result.requirementLevels(missing), ...
            repmat(maximumWidth*(1+0.04*k),nnz(missing),1),'x', ...
            'Color',colors(k,:),'MarkerSize',5,'HandleVisibility','off');
        boundary = widths(:,k)==min(result.config.widthGrid_mrad);
        plot(ax,100*result.requirementLevels(boundary),widths(boundary,k),'v', ...
            'Color',colors(k,:),'MarkerFaceColor','w','HandleVisibility','off');
    end
    legend(ax,h,["Line","Annular"],'Box','off','Location','northwest','FontSize',8);
    ylim(ax,[0,maximumWidth*1.18]);
    xlim(ax,[68,97]);
    xticks(ax,70:5:95);
    discussion.styleAxes(ax,'Required overall P_D (%)','Minimum tested w_D (mrad)', ...
        sprintf('(%c) Pulse energy = %g E_0','a'+a-1,energy(a)));
end
title(layout,sprintf('x above %g: unattained; open triangle: lower grid boundary',maximumWidth), ...
    'FontWeight','normal','FontSize',10);
end

function fig = blindFigure(result,colors)
% Preserve both meanings of blind zone instead of conflating them.
[fig,layout] = discussion.figure("Geometric blindness versus SNR-limited blindness",1,2);
fig.Position = [2,2,25,10];
styles = [":","-","--"];
w = result.config.widthGrid_mrad;
ax = nexttile(layout);
hold(ax,'on');
indices = [2,1,3]; % Coverage cases: PRF half, nominal, double.
scales = [0.5,1,2];
handles = gobjects(6,1);
labels = strings(6,1);
beamNames = ["Line","Annular"];
for j = 1:3
    curve = result.geometricCoverage.curves{indices(j)};
    for k = 1:2
        handles(k+2*(j-1)) = plot(ax,w,100*curve.blind(:,k), ...
            'LineStyle',styles(j),'Color',colors(k,:),'LineWidth',1.3);
        labels(k+2*(j-1)) = beamNames(k)+", f/f_0="+scales(j);
    end
end
legend(ax,handles,labels,'Box','off','FontSize',7,'NumColumns',2,'Location','northwest');
ylim(ax,[0,100]);
discussion.styleAxes(ax,"w_D (mrad)","Geometric blind fraction (%)", ...
    "(a) Point-centre 1/e^2 footprint; no SNR test");
ax = nexttile(layout);
hold(ax,'on');
for j = 1:3
    curve = result.staticCurves{j};
    for k = 1:2
        handles(k+2*(j-1)) = plot(ax,w,100*(1-curve.probability(:,k)), ...
            'LineStyle',styles(j),'Color',colors(k,:),'LineWidth',1.3);
        labels(k+2*(j-1)) = beamNames(k)+", E/E_0="+scales(j);
    end
end
legend(ax,handles,labels,'Box','off','FontSize',7,'NumColumns',2,'Location','northwest');
ylim(ax,[0,100]);
discussion.styleAxes(ax,"w_D (mrad)","SNR-blind volume fraction (%)", ...
    "(b) Finite target; one scan; SNR threshold 2");
title(layout,'Geometric and SNR blind sets use different criteria; neither is the 15 s moving-target miss rate', ...
    'FontWeight','normal','FontSize',10);
end

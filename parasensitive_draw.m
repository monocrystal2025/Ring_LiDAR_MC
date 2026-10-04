
function fig = parasensitive_draw(options)
%PARASENSITIVE_DRAW Separate line/annular metrics with shared row scales.
arguments
    options.Visible (1,1) string {mustBeMember(options.Visible,["on","off"])} = "on"
    options.CompactLegends (1,1) logical = true
end
root = fileparts(mfilename('fullpath'));
s = load(fullfile(root,'parasensitive_predraw.mat'),'plotData');
p = s.plotData;
assert(p.schemaVersion==5,'Run parasensitive_predraw again.');
assert(isfield(p,'meanPulseCount'),'Run parasensitive_predraw again.');
fig = gobjects(3,1);
colors = [0 .4470 .7410;.8500 .3250 .0980]; % advantage_draw.m: line / annular
styles = {'-','--',':','-.'};
markers = {'o','s','^','d'};
xLimits = {[100 300],[2 10],[20 40],[.1e-4 .9e-4]};
xTicks = {100:50:300,2:2:10,20:5:40,.1e-4:.2e-4:.9e-4};
xLabels = {'Pulse energy (\muJ)','Repetition rate (kHz)', ...
    'Target speed (m/s)','\alpha (m^{-1})'};
yLabels = {'\kappa','Blind fraction (%)','Pulse count'};
clusterTop = 2;
pulseTop = 12;
kappaValues = p.kappa(:,[1 2 4],:);
yMax = max(kappaValues,[],'all','omitnan');
if isempty(yMax) || ~isfinite(yMax), yMax=1; end
kappaTop=1.2;%max(1.2,ceil(2.3*yMax*2)/2);
yLimits = {[0.3 kappaTop],[0 80],[0 pulseTop]};
yTicks = {0.3:.3:kappaTop,0:20:80,0:3:12};
displayColumns = [1 2 4];
outputNames = {'parasensitive_kappa.png','parasensitive_blind.png','parasensitive_pulses.png'};
figureSizes = [
    780 250;   % κ 图
    780 250;   % 盲区率图
    800 250    % 脉冲数／脉冲簇数图
];
for row=1:3
    fig(row)=figure('Color','w','Position',[40+30*(row-1), 40+30*(row-1), figureSizes(row,:)], ...
        'Visible',options.Visible,'Name',sprintf('ROE sensitivity: %s',yLabels{row}));
    layout=tiledlayout(fig(row),1,3,'TileSpacing','tight','Padding','compact');
    for col=1:3
        sourceCol=displayColumns(col);
        ax = nexttile(layout,col);
        hold(ax,'on');
        x = p.values{sourceCol}*p.displayScale(sourceCol);
        if sourceCol==4, x=p.values{sourceCol}; end
        if row==1
            handles=gobjects(size(p.requirements,1),1);
            labels=cell(size(p.requirements,1),1);
            conditionColors=lines(size(p.requirements,1));
            for r=1:size(p.requirements,1)
                k=mod(r-1,numel(styles))+1;
                handles(r)=plot(ax,x,p.kappa(:,sourceCol,r), ...
                    'LineStyle',styles{k},'Marker',markers{k}, ...
                    'Color',conditionColors(r,:),'LineWidth',3,'MarkerSize',5);
                labels{r}=sprintf('Detection probability>%g%%\nfirst-warning time<%gs', ...
                    100*p.requirements(r,2),p.requirements(r,1));
                if options.CompactLegends
                    labels{r}=sprintf('P_D>%g%%, t<%gs', ...
                        100*p.requirements(r,2),p.requirements(r,1));
                end
            end
            if col==1, sharedHandles=handles; sharedLabels=labels; end
        elseif row==2
            handles=gobjects(2,1);
            for b=1:2
                if row==2, y=p.meanBlindPct(:,sourceCol,b); else, y=p.meanClusterCount(:,sourceCol,b); end
                handles(b)=plot(ax,x,y,'-o','Color',colors(b,:),'LineWidth',3, ...
                    'MarkerSize',5,'MarkerFaceColor',colors(b,:));
            end
            if col==1, sharedHandles=handles; sharedLabels={'Line','Annular'}; end
        else
            handles=gobjects(4,1);
            yyaxis(ax,'left');
            for b=1:2
                handles(b)=plot(ax,x,p.meanPulseCount(:,sourceCol,b),'-o', ...
                    'Color',colors(b,:),'LineWidth',3,'MarkerSize',5);
            end
            yyaxis(ax,'right');
            for b=1:2
                handles(2+b)=plot(ax,x,p.meanClusterCount(:,sourceCol,b),'--s', ...
                    'Color',colors(b,:),'LineWidth',3,'MarkerSize',5);
            end
            ylim(ax,[0 clusterTop]); yticks(ax,0:0.5:2);
            if col==3, ylabel(ax,'Cluster count'); else, yticklabels(ax,{}); end
            ax.YColor='k'; ax.YTickLabelRotation=0;
            yyaxis(ax,'left'); ax.YColor='k';
            if options.CompactLegends
                % Separate color and line-style keys avoid four verbose labels.
                keys=gobjects(4,1);
                for b=1:2
                    keys(b)=plot(ax,nan,nan,'-','Color',colors(b,:),'LineWidth',3);
                end
                keys(3)=plot(ax,nan,nan,'-o','Color',[.2 .2 .2],'LineWidth',3);
                keys(4)=plot(ax,nan,nan,'--s','Color',[.2 .2 .2],'LineWidth',3);
                if col==1, sharedHandles=keys; sharedLabels={'Line','Annular','Pulses','Clusters'}; end
            else
                if col==1, sharedHandles=handles; sharedLabels={'Line: pulses','Annular: pulses','Line: clusters','Annular: clusters'}; end
            end
        end

        xlim(ax,xLimits{sourceCol}); xticks(ax,xTicks{sourceCol});
        ylim(ax,yLimits{row});
        yticks(ax,yTicks{row}); % Identical manual ticks/grid in every column.
        if col==1, ylabel(ax,yLabels{row}); else, yticklabels(ax,{}); end
        xlabel(ax,xLabels{sourceCol});
        if sourceCol==4, ax.XAxis.Exponent=-4; end
        ax.FontName='Helvetica'; ax.FontSize=14;
        ax.LabelFontSizeMultiplier=1; ax.TitleFontSizeMultiplier=1;
        ax.XTickLabelRotation=0; ax.YTickLabelRotation=0;
        ax.Box='on'; ax.LineWidth=2; grid(ax,'on'); ax.GridAlpha=.10;
        ax.XMinorGrid='off'; ax.YMinorGrid='off';
        ax.XMinorTick='off'; ax.YMinorTick='off';
        ax.Tag=sprintf('parasensitive_r%d_c%d',row,col);
        ax.Toolbar.Visible='off';
        if sourceCol==4
            % Avoid R2025a export clipping of small SI-valued X coordinates.
            set(findall(ax,'Type','line'),'Clipping','off');
        end
    end
    lgd=legend(sharedHandles,sharedLabels,'Orientation','horizontal', ...
        'Box','off','AutoUpdate','off','FontName','Helvetica','FontSize',14);
    lgd.Layout.Tile='north';
    lgd.ItemTokenSize=[25 12];
    lgd.Tag=sprintf('shared_legend_r%d',row);
    set(findall(fig(row),'-property','FontName'),'FontName','Helvetica');
    set(findall(fig(row),'-property','FontSize'),'FontSize',14);
    drawnow;
    % Freeze the laid-out axes in figure coordinates before positioning the
    % shared legend, so tiledlayout cannot restore its default vertical gap.
    axesList=findall(fig(row),'Type','axes');
    axesPositions=zeros(numel(axesList),4);
    for k=1:numel(axesList)
        axesPositions(k,:)=getpixelposition(axesList(k),true);
    end
    legendPosition=getpixelposition(lgd,true);
    for k=1:numel(axesList)
        axesList(k).Parent=fig(row);
        axesList(k).Units='pixels';
        axesList(k).Position=axesPositions(k,:);
    end
    lgd.Parent=fig(row);
    delete(layout);
    firstAxes=findobj(fig(row),'Tag',sprintf('parasensitive_r%d_c1',row));
    lgd.Units='pixels';
    legendPosition(2)=firstAxes.Position(2)+firstAxes.Position(4)+12;
    lgd.Position=legendPosition;
    drawnow;
    exportgraphics(fig(row),fullfile(root,outputNames{row}),'Resolution',300);
    fprintf('Saved %s (300 dpi, N=%d).\n',outputNames{row},p.N);
end
end

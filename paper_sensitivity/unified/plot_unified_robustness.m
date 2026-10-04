function figures = plot_unified_robustness(options)
%PLOT_UNIFIED_ROBUSTNESS Equal visual treatment for three manuscript benefits.
arguments
    options.Directory (1,1) string = ""
    options.Visible (1,1) string {mustBeMember(options.Visible,["on","off"])} = "off"
end
folder=fileparts(mfilename('fullpath'));root=fileparts(fileparts(folder));
if options.Directory=="",options.Directory=fullfile(folder,'results');end
saved=load(fullfile(options.Directory,'analysis.mat'),'analysis');a=saved.analysis;
out=fullfile(folder,'figures');if ~isfolder(out),mkdir(out);end
titles={'Pulse energy','Target reflectivity','Atmospheric extinction'};
labels={'Angular-width reduction (%)','Blind-zone reduction (%)','Effective-pulse increase (%)'};
colors=[0 .45 .36;0 .447 .698;.835 .369 0];
figures=gobjects(0);
for k=1:numel(a.design.options.ReferenceWidths)
    width=a.design.options.ReferenceWidths(k);
    fig=figure('Color','w','Units','centimeters','Position',[2 2 19 16], ...
        'Visible',options.Visible,'Renderer','painters');
    layout=tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
    for row=1:3
        lo=min([0;reshape(a.ci(:,row,1,k),[],1)],[],'omitnan');
        hi=max([5;reshape(a.ci(:,row,2,k),[],1)],[],'omitnan');
        span=max(hi-lo,10);limits=[lo-.08*span,hi+.12*span];
        for p=1:3
            ax=nexttile(layout);hold(ax,'on');
            patch(ax,[-23 23 23 -23],[0 0 limits(2) limits(2)],[.94 .98 .96], ...
                'EdgeColor','none','HandleVisibility','off');
            ids=find(a.design.factorIndex==p | a.design.factorIndex==0);
            [scale,order]=sort(a.design.scale(ids,p));ids=ids(order);x=100*(scale-1);
            y=a.point(ids,row,k);low=a.ci(ids,row,1,k);high=a.ci(ids,row,2,k);
            errorbar(ax,x,y,y-low,high-y,'-o','Color',colors(1,:), ...
                'LineWidth',1.2,'MarkerSize',4,'MarkerFaceColor','w','CapSize',4);
            yline(ax,0,'--','Color',[.45 .45 .45],'LineWidth',.9);
            xline(ax,0,':','Color',[.65 .65 .65]);
            baseline=find(abs(x)<1e-10);
            plot(ax,x(baseline),y(baseline),'d','MarkerSize',5,'Color',colors(1,:),'MarkerFaceColor',colors(1,:));
            xlim(ax,[-23 23]);ylim(ax,limits);xticks(ax,[-20 -10 0 10 20]);
            grid(ax,'on');ax.Layer='top';ax.FontSize=8;ax.FontName='Arial';ax.Box='off';
            letter=char('a'+(row-1)*3+p-1);
            if row==1,title(ax,sprintf('(%s) %s',letter,titles{p}),'FontWeight','normal');
            else,title(ax,sprintf('(%s)',letter),'FontWeight','normal');end
            if p==1,ylabel(ax,labels{row});end
            if row==3,xlabel(ax,'Parameter change (%)');end
            if any(~isfinite(y)),text(ax,0,limits(2)-.05*span,'NA present','HorizontalAlignment','center');end
        end
    end
    title(layout,sprintf('Three parallel benefits | Local perturbations around the baseline'), ...
        'FontSize',10,'FontWeight','normal');
    subtitle(layout,sprintf('Width: P_D(15 s) \\geq 60%%; blind zone and pulses: w_D = %g mrad',width),'FontSize',8);
    xlabel(layout,{sprintf('Positive = annular advantage. Paired 95%% pointwise intervals. N_{moving} = %d; N_{static} = %d.', ...
        a.design.options.N,a.design.options.StaticN), ...
        'Baseline: E_p = 150 \muJ; \rho = 0.5; \alpha = 1.5 \times 10^{-5} m^{-1}.'},'FontSize',8);
    name=sprintf('three_benefits_w%d',width);exportFigure(fig,fullfile(out,name));figures(end+1)=fig; %#ok<AGROW>
end
% Absolute quantities disclose what the normalized benefits represent.
k=find(a.design.options.ReferenceWidths==125,1);
fig=figure('Color','w','Units','centimeters','Position',[2 2 19 16],'Visible',options.Visible,'Renderer','painters');
layout=tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
absoluteLabels={'Required divergence (mrad)','Blind fraction (%)','Pulses / successful window'};
for row=1:3
    factor=1;if row==2,factor=100;end
    allValues=reshape(a.absolute(:,:,row,k),[],1)*factor;
    pad=.15*max(max(allValues)-min(allValues),1);limits=[max(0,min(allValues)-pad),max(allValues)+pad];
    for p=1:3
        ax=nexttile(layout);hold(ax,'on');ids=find(a.design.factorIndex==p | a.design.factorIndex==0);
        [scale,ord]=sort(a.design.scale(ids,p));ids=ids(ord);
        plot(ax,100*(scale-1),a.absolute(ids,1,row,k)*factor,'-o','Color',colors(2,:),'MarkerSize',4,'LineWidth',1.2);
        plot(ax,100*(scale-1),a.absolute(ids,2,row,k)*factor,'-^','Color',colors(3,:),'MarkerSize',4,'LineWidth',1.2);
        xlim(ax,[-23 23]);ylim(ax,limits);xticks(ax,[-20 -10 0 10 20]);grid(ax,'on');ax.FontSize=8;ax.FontName='Arial';
        if row==1,title(ax,titles{p},'FontWeight','normal');end
        if row==1 && p==1,legend(ax,{'Line','Annular'},'Location','best','Box','off');end
        if p==1,ylabel(ax,absoluteLabels{row});end
        if row==3,xlabel(ax,'Parameter change (%)');end
    end
end
title(layout,'Absolute responses underlying the three benefits','FontSize',10,'FontWeight','normal');
exportFigure(fig,fullfile(out,'absolute_responses'));figures(end+1)=fig;
% Preserve all earlier wide-range results, including failures and reversals.
old=load(fullfile(root,'para_sensative_N3000_range2','para_summary.mat'));
previous=load(fullfile(root,'paper_sensitivity','results','paper_sensitivity.mat'),'results');
wide=nan(numel(old.manifest.cases),3);d=find(old.manifest.widthsMrad==125,1);
for j=1:size(wide,1)
    c=old.curves{j};id=old.manifest.cases(j).id;
    match=find(previous.results.oat.Case==id,1);
    wide(j,:)=[100*(1-previous.results.oat.Kappa(match)), ...
        100*(1-c.blind(d,2)/c.blind(d,1)),100*(c.meanPulses(d,2)/c.meanPulses(d,1)-1)];
end
fig=figure('Color','w','Units','centimeters','Position',[2 2 22 16],'Visible',options.Visible,'Renderer','painters');
layout=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
wideTitles={'Energy (\muJ)','Reflectivity','Speed (m/s)','Extinction (10^{-5}/m)'};
for row=1:3
    lo=min([0;wide(:,row)],[],'omitnan');hi=max([5;wide(:,row)],[],'omitnan');span=hi-lo;
    for p=1:4
        ax=nexttile(layout);hold(ax,'on');ids=find([old.manifest.cases.parameterIndex]==p);
        x=[old.manifest.cases(ids).value]*old.manifest.displayScale(p);y=wide(ids,row);
        plot(ax,x,y,'-o','Color',colors(1,:),'MarkerSize',4,'LineWidth',1.1);yline(ax,0,'--','Color',[.4 .4 .4]);
        if p==1
            ax.XScale='log';xlim(ax,[min(x)/1.12 max(x)*1.12]);
        else
            xpad=.07*(max(x)-min(x));xlim(ax,[min(x)-xpad,max(x)+xpad]);
        end
        ylim(ax,[lo-.12*span hi+.2*span]);xticks(ax,x);xticklabels(ax,compose('%g',x));grid(ax,'on');ax.FontSize=7.5;
        for z=find(~isfinite(y)).',text(ax,x(z),hi+.08*span,'NA','HorizontalAlignment','center','FontSize',7);end
        if row==1,title(ax,wideTitles{p},'FontWeight','normal');end
        if p==1,ylabel(ax,labels{row});end
        if row==3,xlabel(ax,wideTitles{p});end
    end
end
title(layout,'Operating boundaries: all earlier wide-range configurations (N = 3000)','FontSize',10,'FontWeight','normal');
xlabel(layout,'Descriptive estimates; NA = requirement unattained, not missing simulation data. Reference width = 125 mrad.','FontSize',8);
exportFigure(fig,fullfile(out,'wide_range_boundaries'));figures(end+1)=fig;
wideTable=table(string({old.manifest.cases.id}).',wide(:,1),wide(:,2),wide(:,3), ...
    'VariableNames',{'Case','AngularReduction','BlindReduction','PulseIncrease'});
writetable(wideTable,fullfile(options.Directory,'wide_range_boundaries.csv'));
end

function exportFigure(fig,path)
exportgraphics(fig,path+".png",'Resolution',400);
exportgraphics(fig,path+".pdf",'ContentType','vector');savefig(fig,path+".fig");
end

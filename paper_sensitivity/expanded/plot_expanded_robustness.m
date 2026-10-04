function handles = plot_expanded_robustness()
%PLOT_EXPANDED_ROBUSTNESS Unified three-row display and transparent screening.
folder=fileparts(mfilename('fullpath'));saved=load(fullfile(folder,'results','analysis.mat'),'analysis');a=saved.analysis;
out=fullfile(folder,'figures');if ~isfolder(out),mkdir(out);end
titles={'Pulse energy','Target reflectivity','Atmospheric extinction','Scan angular speed'};
labels={'Angular-width reduction (%)','Blind-zone reduction (%)','Effective-pulse increase (%)'};
color=[0 .45 .36];active=find(a.halfRange>0);np=numel(active);
assert(np>=3,'Expanded:Domain','Fewer than three parameter domains qualify.');
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 6*np 16],'Renderer','painters');
layout=tiledlayout(fig,3,np,'TileSpacing','compact','Padding','compact');
for metric=1:3
    vals=reshape(a.ci(a.chosen,metric,:),[],1);lo=min([0;vals]);hi=max([5;vals]);span=hi-lo;
    limits=[lo-.08*span hi+.12*span];
    for col=1:np
        p=active(col);ax=nexttile(layout);hold(ax,'on');
        ids=find(a.chosen & (a.design.factorIndex==p | a.design.factorIndex==0));
        [scale,order]=sort(a.design.scale(ids,p));ids=ids(order);x=round(100*(scale-1));h=a.halfRange(p);xlimits=[-h-3 h+3];
        patch(ax,[xlimits fliplr(xlimits)],[0 0 limits(2) limits(2)],[.94 .98 .96],'EdgeColor','none','HandleVisibility','off');
        y=a.point(ids,metric);low=a.ci(ids,metric,1);high=a.ci(ids,metric,2);
        errorbar(ax,x,y,y-low,high-y,'-o','Color',color,'MarkerFaceColor','w','MarkerSize',3.4,'LineWidth',1.15,'CapSize',3);
        base=find(x==0);plot(ax,0,y(base),'d','Color',color,'MarkerFaceColor',color,'MarkerSize',5);
        yline(ax,0,'--','Color',[.4 .4 .4]);xline(ax,0,':','Color',[.6 .6 .6]);
        xlim(ax,xlimits);ylim(ax,limits);xticks(ax,-h:10:h);if h>=60,xticks(ax,unique([-h -60:20:60 h]));end
        ax.FontName='Arial';ax.FontSize=8;grid(ax,'on');ax.Layer='top';
        letter=char('a'+(metric-1)*np+col-1);
        if metric==1,title(ax,sprintf('(%s) %s (\\pm%d%%)',letter,titles{p},h),'FontWeight','normal','FontSize',9);
        else,title(ax,['(' letter ')'],'FontWeight','normal');end
        if col==1,ylabel(ax,labels{metric});end
        if metric==3,xlabel(ax,'Parameter change (%)');end
    end
end
title(layout,'Three annular-beam benefits over expanded symmetric perturbations','FontWeight','normal','FontSize',11);
subtitle(layout,'Angular requirement: P_D(15 s) \geq 60%; blind zone and pulses: w_D = 125 mrad','FontSize',9);
xlabel(layout,{'Positive = annular advantage; paired 95% pointwise intervals. Moving N = 512; static N = 2048.', ...
    'E_0 = 150 \muJ, \rho_0 = 0.5, \alpha_0 = 1.5 \times 10^{-5} m^{-1}, \Omega_0 = 2\pi rad s^{-1}; static coverage: one cycle.'},'FontSize',8);
exportFigure(fig,fullfile(out,'expanded_three_benefits'));handles=fig;
% Absolute responses establish physical directions without asserting that
% relative gain itself must be monotone in every parameter.
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 6*np 16],'Renderer','painters');
layout=tiledlayout(fig,3,np,'TileSpacing','compact','Padding','compact');
absoluteLabels={'Required divergence (mrad)','Blind fraction (%)','Pulses / successful window'};
for metric=1:3
    mult=1;if metric==2,mult=100;end
    vals=reshape(a.absolute(a.chosen,:,metric),[],1)*mult;pad=.15*max(max(vals)-min(vals),1);lim=[max(0,min(vals)-pad),max(vals)+pad];
    for col=1:np
        p=active(col);ax=nexttile(layout);hold(ax,'on');ids=find(a.chosen & (a.design.factorIndex==p | a.design.factorIndex==0));
        [scale,order]=sort(a.design.scale(ids,p));ids=ids(order);x=round(100*(scale-1));
        plot(ax,x,a.absolute(ids,1,metric)*mult,'-o','Color',[0 .447 .698],'LineWidth',1.2,'MarkerSize',3.5);
        plot(ax,x,a.absolute(ids,2,metric)*mult,'-^','Color',[.835 .369 0],'LineWidth',1.2,'MarkerSize',3.5);
        ylim(ax,lim);h=a.halfRange(p);xlim(ax,[-h-3 h+3]);xticks(ax,unique([-h 0 h]));grid(ax,'on');ax.FontSize=8;
        if metric==1,title(ax,titles{p},'FontWeight','normal');end
        if col==1,ylabel(ax,absoluteLabels{metric});end
        if metric==3,xlabel(ax,'Parameter change (%)');end
        if metric==1 && col==1,legend(ax,{'Line','Annular'},'Location','best','Box','off');end
    end
end
title(layout,'Absolute responses in the selected expanded operating ranges','FontSize',11,'FontWeight','normal');
exportFigure(fig,fullfile(out,'expanded_absolute_responses'));handles(end+1)=fig;
% All candidate values remain visible, even when not selected for the main.
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 24 16],'Renderer','painters');
layout=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
for metric=1:3
    vals=a.point(:,metric);vals=vals(isfinite(vals));lo=min([0;vals]);hi=max([5;vals]);span=hi-lo;
    for p=1:4
        ax=nexttile(layout);hold(ax,'on');ids=find(a.design.factorIndex==p | a.design.factorIndex==0);
        [scale,order]=sort(a.design.scale(ids,p));ids=ids(order);x=round(100*(scale-1));y=a.point(ids,metric);
        plot(ax,x,y,'-o','Color',[.5 .5 .5],'LineWidth',1,'MarkerSize',3);
        selected=a.chosen(ids);plot(ax,x(selected),y(selected),'o','Color',color,'MarkerFaceColor',color,'MarkerSize',4);
        yline(ax,0,'--','Color',[.4 .4 .4]);xlim(ax,[min(x)-4 max(x)+4]);ylim(ax,[lo-.1*span hi+.23*span]);
        for z=find(~isfinite(y)).',text(ax,x(z),hi+.1*span,'NA','FontSize',7,'HorizontalAlignment','center');end
        xticks(ax,unique([min(x) 0 max(x)]));grid(ax,'on');ax.FontSize=8;
        if metric==1,title(ax,titles{p},'FontWeight','normal');end
        if p==1,ylabel(ax,labels{metric});end
        if metric==3,xlabel(ax,'Parameter change (%)');end
    end
end
title(layout,'All tested candidate levels, including excluded operating conditions','FontSize',11,'FontWeight','normal');
xlabel(layout,{'Green markers: selected contiguous symmetric domain. Gray: other candidates; point estimates only.', ...
    'NA: consult all_candidates.csv for unattained requirements or reference-screen exclusion; no interpolation of missing outcomes.'},'FontSize',8);
exportFigure(fig,fullfile(out,'candidate_boundaries'));handles(end+1)=fig;
end

function exportFigure(fig,name)
exportgraphics(fig,name+".png",'Resolution',400);exportgraphics(fig,name+".pdf",'ContentType','vector');savefig(fig,name+".fig");
end

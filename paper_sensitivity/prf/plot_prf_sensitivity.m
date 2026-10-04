function plot_prf_sensitivity(referenceOnly,resultDirectory,exportDirectory)
%PLOT_PRF_SENSITIVITY Show all PRF levels, including undefined ratios.
if nargin<1,referenceOnly=false;end
folder=fileparts(mfilename('fullpath'));
if nargin<2,resultDirectory=fullfile(folder,'results');end
if nargin<3,exportDirectory=fullfile(folder,'figures');end
s=load(fullfile(resultDirectory,'analysis.mat'),'a');a=s.a;
deadline=15;target=.6;
if isfield(a,'criterion'),deadline=a.criterion.DeadlineSeconds;target=a.criterion.TargetProbability;end
out=exportDirectory;if ~isfolder(out),mkdir(out);end
f=a.design.frequencies/1000;color=[0 .45 .36];
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 24 14]);
layout=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
labels={'Required divergence (mrad)','Static blind fraction (%)','Pulses / successful window'};
indices=1:3;if referenceOnly,indices=[4 2 3];labels{1}=sprintf('Detection within %g s (%%)',deadline);end
for metric=1:3
    ax=nexttile(layout);hold(ax,'on');mult=1;if metric==2 || indices(metric)==4,mult=100;end
    plot(ax,f,a.absolute(:,1,indices(metric))*mult,'-o','LineWidth',1.2,'MarkerSize',4);
    plot(ax,f,a.absolute(:,2,indices(metric))*mult,'-^','LineWidth',1.2,'MarkerSize',4);
    ylabel(ax,labels{metric});xlabel(ax,'Pulse repetition frequency (kHz)');grid(ax,'on');xlim(ax,[.7 10.3]);xticks(ax,1:10);ax.FontSize=9;
    if metric==1,legend(ax,{'Line','Annular'},'Location','best');end
    title(ax,char('a'+metric-1),'FontWeight','normal');
    if ~referenceOnly && metric==1
        lim=ylim(ax);missing=all(~isfinite(a.absolute(:,:,1)),2);
        ylim(ax,lim);
        for j=find(missing).',text(ax,f(j),lim(1)+.8*diff(lim),'NA','HorizontalAlignment','center','FontSize',8);end
    end
end
gainLabels={'Required-divergence reduction (%)','Blind fraction: line - annular (pp)','Effective-pulse increase (%)'};
if referenceOnly,gainLabels{1}='Detection: annular - line (pp)';end
for metric=1:3
    ax=nexttile(layout);hold(ax,'on');
    if metric==2,y=a.blindDifference(:,1);lo=a.blindDifference(:,2);hi=a.blindDifference(:,3);
    else,y=a.point(:,indices(metric));lo=a.ci(:,indices(metric),1);hi=a.ci(:,indices(metric),2);end
    plot(ax,f,y,'-o','Color',color,'LineWidth',1.2,'MarkerSize',4);
    good=isfinite(lo)&isfinite(hi);errorbar(ax,f(good),y(good),y(good)-lo(good),hi(good)-y(good),'LineStyle','none','Color',color,'CapSize',4);
    yline(ax,0,'--','Color',[.4 .4 .4]);xlim(ax,[.7 10.3]);xticks(ax,1:10);grid(ax,'on');ax.FontSize=9;
    ylabel(ax,gainLabels{metric});xlabel(ax,'Pulse repetition frequency (kHz)');title(ax,char('d'+metric-1),'FontWeight','normal');
    lim=ylim(ax);for j=find(~isfinite(y)).',text(ax,f(j),lim(1)+.1*diff(lim),'NA','HorizontalAlignment','center','FontSize',8);end
end
title(layout,'PRF sensitivity at fixed pulse energy and scan angular speed','FontSize',12,'FontWeight','normal');
if referenceOnly
xlabel(layout,{'All metrics at 125 mrad. E = 150 \muJ; \Omega = 2\pi rad/s; moving N = 512; static N = 2048.', ...
    'Paired 95% pointwise intervals; positive lower-row values favor annular. Static coverage: one complete scan cycle.', ...
    'pp = percentage points; successful-window pulse counts are conditional statistics.'},'FontSize',8);
name='prf_reference';
else
xlabel(layout,{'E = 150 \muJ; \Omega = 2\pi rad/s; moving N = 512; static N = 2048; paired 95% pointwise intervals.', ...
    sprintf('Angular criterion: P_D(%g s) \\geq %g%%; other metrics: 125 mrad. Blind difference uses percentage points to avoid zero denominators.',deadline,100*target), ...
    'NA angular width: requirement unattained in searched grid; all PRF levels retained.'},'FontSize',8);
name='prf_sensitivity';
end
exportgraphics(fig,fullfile(out,[name '.png']),'Resolution',300);
exportgraphics(fig,fullfile(out,[name '.pdf']),'ContentType','vector');savefig(fig,fullfile(out,[name '.fig']));close(fig);
end

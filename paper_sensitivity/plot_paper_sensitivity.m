function figureHandle = plot_paper_sensitivity(options)
%PLOT_PAPER_SENSITIVITY Six panels with raw outcomes and paired uncertainty.
arguments
    options.ResultsFile (1,1) string = ""
    options.OutputDirectory (1,1) string = ""
    options.Visible (1,1) string {mustBeMember(options.Visible,["on","off"])} = "off"
end
folder=fileparts(mfilename('fullpath'));
if options.ResultsFile=="",options.ResultsFile=fullfile(folder,'results','paper_sensitivity.mat');end
if options.OutputDirectory=="",options.OutputDirectory=fullfile(folder,'figures_final');end
if ~isfolder(options.OutputDirectory),mkdir(options.OutputDirectory);end
saved=load(options.ResultsFile,'results');r=saved.results;
assert(~isempty(r.joint) && r.jointMetadata.AllComplete,'Sensitivity:JointMissing','Complete joint results required.');
lineColor=[0 .447 .698];ringColor=[.835 .369 0];jointColor=[0 .48 .38];
figureHandle=figure('Color','w','Units','centimeters','Position',[2 2 19 13.8], ...
    'Visible',options.Visible,'Renderer','painters','Name','Parameter sensitivity and joint robustness');
layout=tiledlayout(figureHandle,2,3,'TileSpacing','compact','Padding','compact');
panels=[1 2 4 5];letters={'(a)','(b)','(c)','(d)'};
titles={'Pulse energy','Target reflectivity','Target speed','Atmospheric extinction'};
labels={'E_p (\muJ)','\rho','v (m s^{-1})','\alpha (10^{-5} m^{-1})'};
for p=1:4
    ax=nexttile(layout,panels(p));hold(ax,'on');
    t=r.oat(r.oat.ParameterIndex==p,:);
    x=t.Value*r.manifest.displayScale(p);
    baseline=r.manifest.baseline(p)*r.manifest.displayScale(p);
    h1=errorbar(ax,x,t.LineWidth,t.LineWidth-t.LineLow,t.LineHigh-t.LineWidth, ...
        '-o','Color',lineColor,'MarkerFaceColor','w','LineWidth',1.1,'MarkerSize',3.6,'CapSize',3);
    h2=errorbar(ax,x,t.RingWidth,t.RingWidth-t.RingLow,t.RingHigh-t.RingWidth, ...
        '-^','Color',ringColor,'MarkerFaceColor','w','LineWidth',1.1,'MarkerSize',4,'CapSize',3);
    xline(ax,baseline,':','Color',[.5 .5 .5],'HandleVisibility','off');
    if p==1
        ax.XScale='log';xlim(ax,[23 950]);
        legend(ax,[h1 h2],{'Line','Annular'},'Location','northeast','Box','off','FontSize',7.5);
    else
        dx=max(x)-min(x);xlim(ax,[min(x)-.09*dx,max(x)+.09*dx]);
    end
    for j=1:height(t)
        names=strings(0,1);
        if isnan(t.LineWidth(j)),names(end+1)="L";end %#ok<AGROW>
        if isnan(t.RingWidth(j)),names(end+1)="R";end %#ok<AGROW>
        if ~isempty(names)
            text(ax,x(j),118,strjoin(names,'/')+" NA",'FontSize',7,'HorizontalAlignment','center','Color',[.45 .45 .45]);
        end
    end
    ylim(ax,[45 125]);yticks(ax,50:25:125);xticks(ax,x);xticklabels(ax,compose('%g',x));
    xlabel(ax,labels{p});
    if p==1 || p==3,ylabel(ax,'Minimum w_D (mrad)');end
    title(ax,[letters{p} ' ' titles{p}],'FontWeight','normal');
    style(ax);
end
for panel=1:2
    ax=nexttile(layout,3*panel);hold(ax,'on');t=r.joint;x=t.Scenario;
    if panel==1
        value=100*t.PdGain;low=100*t.PdSimLow;high=100*t.PdSimHigh;
        label='P_{D,R} - P_{D,L} (pp)';ttl='(e) Joint: timely detection';
    else
        value=100*t.BlindReduction;low=100*t.BlindSimLow;high=100*t.BlindSimHigh;
        label='B_L - B_R (pp)';ttl='(f) Joint: blind-zone reduction';
    end
    bounds=[min(-2,min(low)-1),max(high)+2];
    patch(ax,[-1 max(x)+1 max(x)+1 -1],[0 0 bounds(2) bounds(2)],[.94 .98 .96], ...
        'EdgeColor','none','HandleVisibility','off');
    errorbar(ax,x(2:end),value(2:end),value(2:end)-low(2:end),high(2:end)-value(2:end), ...
        'o','Color',jointColor,'MarkerFaceColor',jointColor,'LineWidth',.75,'MarkerSize',2.8,'CapSize',2);
    errorbar(ax,x(1),value(1),value(1)-low(1),high(1)-value(1), ...
        'd','Color',[.2 .2 .2],'MarkerFaceColor',[.2 .2 .2],'LineWidth',1,'MarkerSize',4,'CapSize',3);
    yline(ax,0,'--','Color',[.4 .4 .4],'LineWidth',.8);
    xlim(ax,[-1 max(x)+1]);ylim(ax,bounds);xticks(ax,unique([0 6:6:max(x) max(x)]));
    xlabel(ax,'Joint scenario (0 = baseline)');ylabel(ax,label);title(ax,ttl,'FontWeight','normal');
    style(ax);
end
title(layout,sprintf('Equal requirement: P_D(%g s) \\geq %g%%   |   Joint test: w_D = %g mrad', ...
    r.options.Deadline,100*r.options.RequiredProbability,r.jointMetadata.WidthMrad), ...
    'FontSize',9,'FontWeight','normal');
xlabel(layout,{sprintf('OAT: N = %d; joint: N_{moving} = %d, N_{static} = %d; four factors within +/- %g%%.', ...
    r.manifest.N,r.jointMetadata.N,r.jointMetadata.StaticN,100*r.jointMetadata.HalfRange), ...
    'Paired bootstrap: pointwise 95% intervals in (a-d); simultaneous 95% bands in (e-f).'},'FontSize',7);
drawnow;
file=fullfile(options.OutputDirectory,'paper_parameter_sensitivity');
exportgraphics(figureHandle,file+".png",'Resolution',400);
exportgraphics(figureHandle,file+".pdf",'ContentType','vector');
savefig(figureHandle,file+".fig");
end

function style(ax)
set(ax,'FontName','Arial','FontSize',8,'LineWidth',.7,'TickDir','out','Box','off', ...
    'XGrid','off','YGrid','on','GridAlpha',.13,'Layer','top');
end

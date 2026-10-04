function figures = plotBackground(result)
%PLOTBACKGROUND Background-faceted Pd(wD), with explicit receiver controls.
arguments
    result (1,1) struct
end
cfg = result.config;
colors = [0,0.447,0.698;0.835,0.369,0];
[fig1,layout] = discussion.figure("Background-dependent angular-width response",2,2);
fig1.Position = [2,2,22,15];
for j = 1:numel(result.backgroundScales)
    ax = nexttile(layout);
    curve = result.curves{j,1};
    h = discussion.probabilityPlot(ax,curve.widths_mrad,isfinite(curve.times),colors);
    yline(ax,100*cfg.targetProbability,':','80%','HandleVisibility','off');
    iso = result.minima{j,1};
    labels = ["Line","Annular"];
    for k = 1:2
        if isfinite(iso.widths_mrad(k))
            labels(k) = labels(k)+": w_D^*="+iso.widths_mrad(k)+" mrad";
        else
            labels(k) = labels(k)+": 80% unattained";
        end
    end
    location = 'southeast';
    if j==4, location = 'northeast'; end
    legend(ax,h,labels,'Box','off','Location',location,'FontSize',8);
    discussion.styleAxes(ax,"w_D (mrad)","Overall P_D (%)", ...
        sprintf('(%c) L_{sky}/L_0 = %g','a'+j-1,result.backgroundScales(j)));
    xlim(ax,[min(cfg.widthGrid_mrad),max(cfg.widthGrid_mrad)]);
    ylim(ax,[0,100]);
end
title(layout,sprintf('Shape-matched receiver; L_0 = %.1g W m^{-2} sr^{-1} nm^{-1}; N=%d', ...
    cfg.system.skyRadiance_W_m2_sr_nm,cfg.N),'FontWeight','normal','FontSize',10);

[fig2,layout2] = discussion.figure("Receiver-control width response",2,3);
fig2.Position = [2,2,27,15];
receiverLabels = ["Shape matched","Own circular envelope","Common fixed circular FOV"];
for row = 1:2
    bgIndex = row+1; % 10 and 100, declared before computing.
    for col = 1:3
        ax = nexttile(layout2);
        curve = result.curves{bgIndex,col};
        h = discussion.probabilityPlot(ax,curve.widths_mrad,isfinite(curve.times),colors);
        yline(ax,80,':','HandleVisibility','off');
        discussion.styleAxes(ax,"w_D (mrad)","Overall P_D (%)", ...
            sprintf('(%c) %s; sky x%g','a'+(row-1)*3+col-1, ...
            receiverLabels(col),result.backgroundScales(bgIndex)));
        ylim(ax,[0,100]);
        xlim(ax,[min(cfg.widthGrid_mrad),max(cfg.widthGrid_mrad)]);
        if row==1 && col==1
            legend(ax,h,["Line","Annular"],'Box','off','Location','southeast');
        end
    end
end
commonHalf = max(cfg.widthGrid_mrad)/2+3*cfg.system.beamWidth_rad*1e3;
title(layout2,sprintf('Common receiver is fixed at %.1f mrad half field across the width sweep',commonHalf), ...
    'FontWeight','normal','FontSize',10);
figures = [fig1,fig2];

% Preserve the previous jitter/gap/phase checks in a separate, explicitly
% finite-deadline figure. These archived experiments use T=15 s and N=128.
oldPath = fullfile(cfg.root,"discussion_results",cfg.mode,"tmp_4","tmp_4.mat");
if isfile(oldPath)
    saved = load(oldPath,'result');
    old = saved.result;
    [fig3,layout3] = discussion.figure("Other engineering checks: 15 s task",1,3);
    fig3.Position = [2,2,26,8];
    ax = nexttile(layout3);
    h = discussion.probabilityPlot(ax,old.jitterScales,old.jitterHit,colors);
    legend(ax,h,["Line","Annular"],'Box','off','Location','southwest','FontSize',8);
    discussion.styleAxes(ax,"\sigma_{point}/w_d","P_{15 s} (%)","(a) Common pointing jitter");
    ax = nexttile(layout3);
    plot(ax,old.gaps_deg,old.gapProbability*100,'-o','LineWidth',1.2,'MarkerSize',3);
    legend(ax,["Lost energy","Fixed total energy"],'Box','off','FontSize',8,'Location','southwest');
    discussion.styleAxes(ax,"Missing sector (deg)","Annular P_{15 s} (%)","(b) Closure loss");
    ylim(ax,[0,100]);
    ax = nexttile(layout3);
    bar(ax,old.phaseWindowProbability*100);
    xticks(ax,1:3);
    xticklabels(ax,["Original phase","Full cycle","30 ms window"]);
    xtickangle(ax,20);
    discussion.styleAxes(ax,"Model setting","P_{15 s} (%)","(c) Phase and window");
    ylim(ax,[0,100]);
    title(layout3,sprintf('Archived controls, unchanged data: w_D=120 mrad, N=%d, deadline 15 s', ...
        old.config.N),'FontWeight','normal','FontSize',10);
    figures = [figures,fig3];
end
end

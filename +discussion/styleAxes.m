function styleAxes(ax,xText,yText,panel)
%STYLEAXES Readable labels and unobtrusive journal-style axes.
arguments
    ax (1,1) matlab.graphics.axis.Axes
    xText (1,1) string
    yText (1,1) string
    panel (1,1) string = ""
end
set(ax,'FontName','Arial','FontSize',9,'LineWidth',0.7, ...
    'TickDir','out','Box','off','XGrid','off','YGrid','on', ...
    'GridAlpha',0.12,'Layer','top');
xlabel(ax,xText,'FontSize',10);
ylabel(ax,yText,'FontSize',10);
if strlength(panel)>0
    title(ax,panel,'FontSize',10,'FontWeight','normal');
end
end

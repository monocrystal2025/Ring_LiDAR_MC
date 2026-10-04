function handles = probabilityPlot(ax,x,hits,colors)
%PROBABILITYPLOT Probability curves with pointwise 95% Wilson intervals.
% hits is N x number-of-x-values x number-of-beams.
arguments
    ax (1,1) matlab.graphics.axis.Axes
    x double
    hits
    colors double
end
hold(ax,'on');
n = size(hits,1);
handles = gobjects(1,size(hits,3));
styles = ["--","-"];
markers = ["s","o"];
for k = 1:size(hits,3)
    p = reshape(mean(hits(:,:,k),1),1,[]);
    [lo,hi] = discussion.wilson(p*n,n);
    patch(ax,[x(:).',fliplr(x(:).')],[lo,fliplr(hi)]*100, ...
        colors(k,:),'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
    handles(k) = plot(ax,x,p*100,'Color',colors(k,:), ...
        'LineStyle',styles(k),'Marker',markers(k),'MarkerSize',3.5,'LineWidth',1.3);
end
ylim(ax,[0,100]);
end

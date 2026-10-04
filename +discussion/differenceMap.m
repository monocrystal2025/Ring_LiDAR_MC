function differenceMap(ax,values)
%DIFFERENCEMAP Symmetric blue-white-orange color scale centered at zero.
arguments
    ax (1,1) matlab.graphics.axis.Axes
    values double
end
blue = [0,0.447,0.698];
orange = [0.835,0.369,0];
fraction = linspace(0,1,128).';
map = [blue+(1-blue).*fraction;1+(orange-1).*fraction];
colormap(ax,map);
limit = max(abs(values(isfinite(values))));
if isempty(limit) || limit==0
    limit = 1;
end
clim(ax,[-limit,limit]);
set(ax,'XGrid','off','YGrid','off');
end

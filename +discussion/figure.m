function [fig,layout,colors] = figure(name,rows,columns)
%FIGURE Consistent journal-size figures with colorblind-friendly styling.
arguments
    name (1,1) string
    rows (1,1) double {mustBeInteger,mustBePositive}
    columns (1,1) double {mustBeInteger,mustBePositive}
end
fig = figure('Color','w','Name',name,'NumberTitle','off', ...
    'Units','centimeters','Position',[2,2,18,max(6.5,6.1*rows)]);
layout = tiledlayout(fig,rows,columns,'TileSpacing','compact','Padding','compact');
colors = [0,0.447,0.698;0.835,0.369,0;0,0.62,0.45;0.45,0.45,0.45];
end

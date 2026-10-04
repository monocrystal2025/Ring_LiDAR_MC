function fig = directionDiagram()
%DIRECTIONDIAGRAM Distinguish entry heading, field rotation and relative track.
[fig,layout] = discussion.figure("Direction definitions",1,2);
fig.Position = [2,2,25,12];
ax = nexttile(layout);
hold(ax,'on');
[az,el] = meshgrid(linspace(0,2*pi,61),linspace(0,pi/2,21));
surf(ax,cos(el).*cos(az),cos(el).*sin(az),sin(el), ...
    'FaceColor',[0.75,0.82,0.88],'FaceAlpha',0.12,'EdgeColor','none');
theta = linspace(0,2*pi,200);
plot3(ax,cos(theta),sin(theta),zeros(size(theta)),'-','Color',[0.5,0.5,0.5]);
n = [sqrt(3)/2,0,0.5];
east = [0,1,0];
north = cross(n,east);
p = n;
headings = [0,90,180];
colors = [0,0.447,0.698;0,0.62,0.45;0.835,0.369,0];
handles = gobjects(1,3);
for j = 1:3
    v = -0.5*n+sqrt(3)/2*(cosd(headings(j))*north+sind(headings(j))*east);
    sphereExit = -2*dot(p,v);
    groundExit = Inf;
    if v(3)<0
        groundExit = -p(3)/v(3);
    end
    endpoint = min(sphereExit,groundExit);
    t = linspace(0,endpoint,100).';
    path = p+t.*v;
    handles(j) = plot3(ax,path(:,1),path(:,2),path(:,3), ...
        'Color',colors(j,:),'LineWidth',2);
    finish = path(end,:);
    plot3(ax,finish(1),finish(2),finish(3),'o','Color',colors(j,:),'MarkerFaceColor',colors(j,:));
end
plot3(ax,p(1),p(2),p(3),'ko','MarkerFaceColor','k');
text(ax,p(1)+0.05,p(2),p(3),'Entry point','FontSize',9);
legend(ax,handles,["Heading 0 deg","Heading 90 deg","Heading 180 deg"], ...
    'Box','off','Location','northwest','FontSize',8);
axis(ax,'equal');
view(ax,38,23);
xlabel(ax,'x/R');ylabel(ax,'y/R');zlabel(ax,'z/R');
title(ax,'(a) Same radial incidence; ground clips the chord','FontWeight','normal','FontSize',10);
ax.FontName = 'Arial';
ax.FontSize = 9;
grid(ax,'on');
ax = nexttile(layout);
hold(ax,'on');
radius = 60;
plot(ax,radius*cos(theta),radius*sin(theta),'-','Color',[0.6,0.6,0.6],'LineWidth',1.3);
rotations = [0,30,90];
handles = gobjects(1,3);
for j = 1:3
    direction = [sind(rotations(j)),cosd(rotations(j))];
    handles(j) = plot(ax,[-radius,radius]*direction(1), ...
        [-radius,radius]*direction(2),'Color',colors(j,:),'LineWidth',2);
end
quiver(ax,-72,-20,140,0,0,'k','LineWidth',1,'MaxHeadSize',0.12);
text(ax,-72,-28,'Illustrative local relative track','FontSize',8);
arc = linspace(0,pi/6,60);
plot(ax,22*sin(arc),22*cos(arc),'k-');
text(ax,12,28,'\phi=30 deg','FontSize',9);
legend(ax,handles,["Original long axis","Long axis rotated 30 deg","Long axis rotated 90 deg"], ...
    'Box','off','Location','southoutside','FontSize',8);
axis(ax,'equal');
xlim(ax,[-80,80]);ylim(ax,[-72,72]);
discussion.styleAxes(ax,'Original local short coordinate (mrad)', ...
    'Original local long coordinate (mrad)','(b) Field rotation; w_D=120 mrad');
title(layout,'Entry heading changes the target path; field rotation keeps the path fixed', ...
    'FontWeight','normal','FontSize',11);
end

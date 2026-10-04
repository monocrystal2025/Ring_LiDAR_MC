function files = draw_fig1(outputDir)
%DRAW_FIG1 Draw three manuscript panels with native MATLAB graphics.
%   draw_fig1 exports three PNGs and editable FIGs under output/fig1_matlab.
%   draw_fig1(outputDir) selects another output directory.
%
% Geometry sources:
%   MC_EA -> MC_func -> NPY written by path_init/lidarpath.
%   MC_ROI_NEW -> MC_func_ROI_NEW -> build_ROI_scan_cycle.
% Phi is the polar angle from +Z. ROI theta is a half-angle.
% The line-beam long axis is perpendicular to the local trajectory tangent
% (user-requested convention, distinct from the old ROE meridional frame).
% All footprints are linear homothetic expansions from the LiDAR origin.
% Line thickness and annular thickness are enlarged for visibility only.
%
% The orthographic projection is used only to compose the illustration;
% all ROE scan and beam geometry is calculated first in 3-D coordinates.
% ROI panel c uses the inverse angular-plane map of the actual spherical
% samples, without clipping or replacing any boundary connector.
if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(fileparts(mfilename('fullpath')),'output','fig1_matlab');
end
if ~isfolder(outputDir), mkdir(outputDir); end
cfg.ink = [0.045 0.16 0.28];
cfg.path = [0.10 0.12 0.14];
cfg.blue = [0.25 0.48 0.67];
cfg.green = [0.05 0.49 0.30];
cfg.gray = [0.46 0.50 0.55];
cfg.f = 5000;
cfg.omega = 2*pi;
cfg.roeWidth = 0.325;
cfg.ringInnerRatio = 0.83; % Shared visible thickness: line width = Ro - Ri.
cfg.roiWidth = 0.1;
cfg.roiAngle = 30;
cfg.dpi = 300;
cfg.origin = [590 640];
cfg.scale = 430;
az = 0.55;
cfg.projection = [cos(az),sin(az),0; ...
    -0.34*sin(az),0.34*cos(az),-sqrt(1-0.34^2)];
cfg.npy = 'G:\BeamVEC_NEW\beam_vector_fast_5000Hz_h3000m_D325000u.npy';
assert(isfile(cfg.npy),'draw_fig1:MissingNpy', ...
    'ROE NPY file not found: %s. Update cfg.npy.',cfg.npy);
beam = readPointingNpy(cfg.npy);
assert(max(abs(vecnorm(beam,2,2)-1))<1e-10, ...
    'draw_fig1:NotUnit','NPY rows must be unit direction vectors.');
[~,ra] = build_ROI_scan_cycle(cfg.roiAngle,cfg.roiWidth, ...
    cfg.f,cfg.omega,'raster');
[~,sp] = build_ROI_scan_cycle(cfg.roiAngle,cfg.roiWidth, ...
    cfg.f,cfg.omega,'spiral');
figs = gobjects(1,3);
[figs(1),audit] = drawRoe(beam,cfg);
figs(2) = drawRoi(cfg);
figs(3) = drawTrajectories(ra,sp,cfg);
names = {'Fig1a_ROE','Fig1b_ROI','Fig1c_Trajectories'};
files = strings(1,3);
for j = 1:3
    drawnow;
    files(j) = fullfile(outputDir,[names{j},'.png']);
    set(figs(j),'PaperPositionMode','auto','InvertHardcopy','off');
    print(figs(j),char(files(j)),'-dpng',['-r',num2str(cfg.dpi)]);
    savefig(figs(j),fullfile(outputDir,[names{j},'.fig']));
end
audit.phiDegrees = rad2deg(acos(beam([1,end],3)));
audit.roeNpy = cfg.npy;
audit.roiRadius = deg2rad(cfg.roiAngle);
audit.scanRadius = deg2rad(cfg.roiAngle)+cfg.roiWidth/2;
audit.raMaxStep = max(vecnorm(diff(ra),2,2));
audit.spMaxStep = max(vecnorm(diff(sp),2,2));
assert(audit.raMaxStep <= 1.01*cfg.omega/cfg.f, ...
    'draw_fig1:RaContinuity','Unexpected jump in raster samples.');
save(fullfile(outputDir,'geometry_audit.mat'),'audit','ra','sp','cfg');
fprintf('Created three MATLAB figures:\n%s\n%s\n%s\n',files);
fprintf('ROE endpoint polar angles: %.6f, %.6f deg\n',audit.phiDegrees);
fprintf('Line long-axis dot tangent: %.3g\n',audit.lineOrthogonality);
end

function [fig,ax] = canvas(tag,~)
fig = figure('Name',['Fig. 1 ',tag],'Color','w', ...
    'Units','pixels','Position',[90 80 800 650], ...
    'NumberTitle','off');
ax = axes(fig,'Position',[0 0 1 1],'XLim',[0 1200], ...
    'YLim',[0 850],'YDir','reverse','DataAspectRatio',[1 1 1]);
hold(ax,'on'); axis(ax,'off');
% Keep coplanar vector elements above translucent surfaces in draw order.
ax.SortMethod = 'childorder';

end

function p = project(v,cfg)
p = cfg.origin + cfg.scale*(v*cfg.projection.');
end

function h = line2(ax,p,color,width,varargin)
h = plot(ax,p(:,1),p(:,2),'-','Color',color,'LineWidth',width,varargin{:});
end

function flat(ax,p,color,alpha,edge,width)
patch(ax,'XData',p(:,1),'YData',p(:,2),'FaceColor',color, ...
    'FaceAlpha',alpha,'EdgeColor',edge,'LineWidth',width);
end

function arrow2(ax,p,q,color,width)
d = q-p; d=d/norm(d); n=[-d(2),d(1)];
line2(ax,[p;q],color,width);
flat(ax,[q;q-10*d+4*n;q-10*d-4*n],color,1,'none',0.5);
end

function v = unit(v)
v = v/norm(v);
end

function points = footprint(u,radius)
t = linspace(0,2*pi,181).';
v = unit(cross(u,[0 0 1])); w=cross(u,v);
points = u+radius*(cos(t).*v+sin(t).*w);
end

function p = endpointCircle(u,cfg)
r = hypot(u(1),u(2));
t0 = atan2(u(2),u(1));
t = (t0:cfg.omega/cfg.f/r:t0+2*pi).';
p = [r*cos(t),r*sin(t),repmat(u(3),numel(t),1)];
% The tiny last-to-first pulse interval is displayed as a closed circle.
p = [p;p(1,:)];
end

function [fig,audit] = drawRoe(beam,cfg)
[fig,ax] = canvas('(a)',cfg);

% Ground and shaded hemisphere, matching entry_heading_explanation.png.
ground = [-.98 -.98 0;.98 -.98 0;.98 .98 0;-.98 .98 0];
flat(ax,project(ground,cfg),[.68 .70 .74],.20,'none',.5);
[az,el] = meshgrid(linspace(0,2*pi,100),linspace(0,pi/2,42));
x=cos(el).*cos(az); y=cos(el).*sin(az); z=sin(el);
pp=project([x(:),y(:),z(:)],cfg);
shade=.70+.12*z-.08*x+.04*y;
surf(ax,reshape(pp(:,1),size(x)),reshape(pp(:,2),size(x)), ...
    -5*ones(size(x)),shade,'FaceColor','interp','EdgeColor','none','FaceAlpha',.43);
colormap(ax,[linspace(.44,.84,128).',linspace(.66,.92,128).',linspace(.83,.99,128).']);
% Light geographic guides do not obscure the black scan path.
for a=0:pi/3:5*pi/3
    e=linspace(0,pi/2,100).';
    p=project([cos(e)*cos(a),cos(e)*sin(a),sin(e)],cfg);
    line2(ax,p,[.66 .74 .80],.55,'LineStyle','--');
end
t=linspace(0,2*pi,401).';
line2(ax,project([cos(t),sin(t),zeros(size(t))],cfg),cfg.blue,.85);
line2(ax,[cfg.origin;project([0 0 1.15],cfg)],cfg.gray,.85,'LineStyle','--');
line2(ax,project(beam,cfg),cfg.path,1.15);
top=endpointCircle(beam(1,:),cfg); bottom=endpointCircle(beam(end,:),cfg);
line2(ax,project(top,cfg),cfg.path,1.5);
line2(ax,project(bottom,cfg),cfg.path,1.5);
targets = [330 445;545 315;850 465];
pp=project(beam,cfg);
% Avoid grazing projections: choose visible front-facing footprints so that
% their area centroid, rather than an edge, visibly meets the scan path.
towardViewer=unit(cross(cfg.projection(1,:),cfg.projection(2,:)));
frontCandidates=find(beam*towardViewer.' > .55);
for j=1:3
    [~,nearest]=min(vecnorm(pp(frontCandidates,:)-targets(j,:),2,2));
    idx=frontCandidates(nearest);
    u=beam(idx,:); center=pp(idx,:);
    tangent=beam(min(idx+1,end),:)-beam(max(idx-1,1),:);
    tangent=unit(tangent-dot(tangent,u)*u);
    long=unit(cross(tangent,u));
    % At unit propagation range, all maximum dimensions equal 2*r.
    r=tan(cfg.roeWidth/2);
    ringInnerRadius=cfg.ringInnerRatio*r;
    sharedFullWidth=r-ringInnerRadius; % Radial ring thickness, not twice it.
    if j==2
        shortWidth=sharedFullWidth/2;
        ends=[u-r*long-shortWidth*tangent;u+r*long-shortWidth*tangent; ...
            u+r*long+shortWidth*tangent;u-r*long+shortWidth*tangent];
        audit.lineOrthogonality=dot(long,tangent);
    else
        ends=footprint(u,r);
    end
    ep=project(ends,cfg);
    if j==2
        halfShort=shortWidth;
    else
        halfShort=0;
    end
    % Project actual 3-D propagation support, not a filled 2-D hull.
    % A ray is in the beam only when 0<=s<=1 and its transverse
    % coordinates satisfy the solid / rectangular / annular cross-section.
    drawBeamVolume(ax,u,long,tangent,r,ringInnerRadius,halfShort,j,cfg);
    if j==3
        inner=project(footprint(u,ringInnerRadius),cfg);
        vertices=[ep;inner]; n=size(ep,1);
        faces=[(1:n-1).',(2:n).',(n+2:2*n).',(n+1:2*n-1).'];
        patch(ax,'Faces',faces,'Vertices',vertices,'FaceColor',cfg.green, ...
            'FaceAlpha',.85,'EdgeColor','none');
        line2(ax,ep,cfg.green,1.2);line2(ax,inner,cfg.green,.8);
    else
        flat(ax,ep,cfg.green,.75,cfg.green,.9);
    end
    % Compute the projected area centroid independently from the path point.
    areaShape=polyshape(ep(:,1),ep(:,2),'Simplify',false);
    if j==3
        areaShape=subtract(areaShape,polyshape(inner(:,1),inner(:,2),'Simplify',false));
    end
    [cx,cy]=centroid(areaShape);
    audit.projectedCentroidError(j)=norm([cx,cy]-pp(idx,:));
    assert(audit.projectedCentroidError(j)<1e-8, ...
        'draw_fig1:OffPathCenter','Footprint centroid must be on the scan.');
    % Repaint the local trajectory on top of the translucent illumination.
    local=max(1,idx-100):min(size(beam,1),idx+100);
    line2(ax,pp(local,:),cfg.path,1.25);
    plot(ax,center(1),center(2),'.','Color',cfg.path,'MarkerSize',12);
    drone(ax,center+[0 -56],.8,cfg);
    audit.maximumDimension(j)=2*r;
    if j==2
        audit.lineFullLength=norm(ends(2,:)-ends(1,:));
        audit.lineFullWidth=norm(ends(3,:)-ends(2,:));
    elseif j==3
        audit.ringOuterDiameter=2*r;
        audit.ringRadialWidth=r-ringInnerRadius;
    else
        audit.spotDiameter=2*r;
    end
    audit.beamSample(j)=idx;
    audit.centerError(j)=norm(mean([ends(1,:);ends(1+floor(size(ends,1)/2),:)],1)-u);
end
assert(abs(audit.lineFullWidth-audit.ringRadialWidth)<1e-12, ...
    'draw_fig1:WidthMismatch','Line width must equal radial ring width.');
assert(max(abs([audit.spotDiameter,audit.lineFullLength, ...
    audit.ringOuterDiameter]-2*r))<1e-12, ...
    'draw_fig1:DiameterMismatch','The three maximum dimensions must match.');
% Polar angle from +Z, drawn to the actual annular beam center direction.
u=beam(audit.beamSample(3),:);
phi=acos(u(3)); a=atan2(u(2),u(1)); tt=linspace(0,phi,80).';
arc=.25*[sin(tt)*cos(a),sin(tt)*sin(a),cos(tt)];
line2(ax,project(arc,cfg),cfg.ink,1.05);
plot(ax,cfg.origin(1),cfg.origin(2),'o','Color',cfg.ink, ...
    'MarkerFaceColor',cfg.ink,'MarkerSize',5);


end

function fig = drawRoi(cfg)
[fig,ax]=canvas('(b)',cfg);

o=[135 510]; axisVector=[790 -45]; t=linspace(0,2*pi,241).';
base=sin(deg2rad(cfg.roiAngle))*[66*cos(t),450*sin(t)];
center=o+axisVector; rim=center+base;
points=[o;rim]; h=convhull(points);
flat(ax,points(h,:),[.54 .76 .91],.28,cfg.blue,1.1);
flat(ax,rim,[.54 .76 .91],.16,cfg.blue,1);
line2(ax,[o;center],cfg.gray,.85,'LineStyle','--');
for z=[.43 .71]
    line2(ax,o+z*(rim-o),[.52 .63 .71],.75,'LineStyle','--');
end
bc=center+[0 -100]; edge=bc+[13*cos(t),23*sin(t)];
hull=convhull([o;edge]); pts=[o;edge];
flat(ax,pts(hull,:),cfg.green,.14,'none',.5);
for z=[.43 .71 1]
    flat(ax,o+z*(edge-o),cfg.green,.75,cfg.green,.75);
end
dronePos=o+.69*axisVector+[0 39];
drone(ax,dronePos,1,cfg);
arrow2(ax,dronePos+[38 -19],dronePos+[137 -117],cfg.green,1.5);
plot(ax,o(1),o(2),'o','MarkerFaceColor',cfg.ink,'Color',cfg.ink,'MarkerSize',5);
angles=atan2(rim(:,2)-o(2),rim(:,1)-o(1));
q=linspace(min(angles),max(angles),70).';
line2(ax,o+118*[cos(q),sin(q)],cfg.ink,.9);

end

function fig = drawTrajectories(ra,sp,cfg)
[fig,ax]=canvas('(c)',cfg);

theta=deg2rad(cfg.roiAngle); scale=360;
t=linspace(0,2*pi,361).'; paths={ra,sp};
centers=[305 421;895 421];
for j=1:2
    c=centers(j,:);
    disc=c+theta*scale*[cos(t),sin(t)];
    flat(ax,disc,[.61 .84 .75],.26,[.62 .75 .70],.6);
    b=paths{j};
    polar=acos(max(-1,min(1,b(:,3)))); az=atan2(b(:,2),b(:,1));
    angular=[polar.*cos(az),polar.*sin(az)];
    pixels=c+scale*angular.*[1 -1];
    % No scan-envelope circle: avoid disguising the real boundary arcs.
    line2(ax,pixels,cfg.path,2);
end


end

function drawBeamVolume(ax,u,long,short,r,ri,halfShort,kind,cfg)
%DRAWBEAMVOLUME Orthographic ray projection of finite 3-D beam support.
% Green opacity encodes geometric path length through the schematic beam
% volume. It is NOT a radiometric prediction or visible laser photograph.
% Unlike an outer convex-hull fill, this preserves the annular empty core.
camera=unit(cross(cfg.projection(1,:),cfg.projection(2,:)));
if kind==2
    rim=[u-r*long-halfShort*short;u+r*long-halfShort*short; ...
        u+r*long+halfShort*short;u-r*long+halfShort*short];
else
    rim=footprint(u,r);
end
p=[cfg.origin;project(rim,cfg)];
lo=floor(min(p,[],1))-2; hi=ceil(max(p,[],1))+2;
xx=linspace(lo(1),hi(1),max(20,ceil(hi(1)-lo(1))+1));
yy=linspace(lo(2),hi(2),max(20,ceil(hi(2)-lo(2))+1));
[x,y]=meshgrid(xx,yy);
screen=[(x(:)-cfg.origin(1))/cfg.scale, ...
    (y(:)-cfg.origin(2))/cfg.scale];
q=screen*cfg.projection;
s0=q*u.'; ns=dot(camera,u);
% s = dot(q+t*camera,u), and each section lies in a plane normal to u.
tLower=(0-s0)/ns; tUpper=(1-s0)/ns;
assert(ns>0,'draw_fig1:Camera','Selected directions must face viewer.');
if kind==2
    [tLower,tUpper]=clipLinear(tLower,tUpper,q*long.'-r*s0, ...
        dot(camera,long)-r*ns);
    [tLower,tUpper]=clipLinear(tLower,tUpper,-q*long.'-r*s0, ...
        -dot(camera,long)-r*ns);
    [tLower,tUpper]=clipLinear(tLower,tUpper,q*short.'-halfShort*s0, ...
        dot(camera,short)-halfShort*ns);
    [tLower,tUpper]=clipLinear(tLower,tUpper,-q*short.'-halfShort*s0, ...
        -dot(camera,short)-halfShort*ns);
    distance=max(0,tUpper-tLower);
else
    distance=coneChord(q,s0,ns,r,tLower,tUpper);
    if kind==3
        distance=max(0,distance-coneChord(q,s0,ns,ri,tLower,tUpper));
    end
end
opacity=reshape(.18*distance/max(max(distance),eps),size(x));
rgb=repmat(reshape(cfg.green,1,1,3),size(x,1),size(x,2));
image(ax,'XData',[xx(1) xx(end)],'YData',[yy(1) yy(end)], ...
    'CData',rgb,'AlphaData',opacity,'AlphaDataMapping','none');
% No internal generators or arbitrary intermediate cross-section lines.
end

function [lo,hi]=clipLinear(lo,hi,b,a)
% Intersect the ray interval with b+a*t <= 0.
if abs(a)<1e-12
    hi(b>0)=-Inf;
elseif a>0
    hi=min(hi,-b/a);
else
    lo=max(lo,-b/a);
end
end

function length=coneChord(q,s0,ns,r,lo,hi)
% |q+t*n|^2 - (1+r^2)*(s0+t*ns)^2 <= 0.
% q lies in the screen plane, hence dot(q,n)=0.
a=1-(1+r^2)*ns^2;
b=-2*(1+r^2)*s0*ns;
c=sum(q.^2,2)-(1+r^2)*s0.^2;
disc=b.^2-4*a*c;
if abs(a)<1e-12
    % Tangential viewing case; this figure deliberately avoids it.
    error('draw_fig1:TangentCamera','Choose a non-tangent beam direction.');
end
root=sqrt(max(disc,0));
r1=(-b-root)/(2*a);r2=(-b+root)/(2*a);
left=min(r1,r2);right=max(r1,r2);
inside=max(0,min(hi,right)-max(lo,left));
inside(disc<0)=0;
if a>0
    length=inside;
else
    length=max(0,hi-lo)-inside;
end
length=max(0,length);
end

function drone(ax,center,scale,~)
% Shared mechanical quadcopter vector: actual two-bladed propellers.
% No enclosed rotor ellipses or central dots that resemble eyes.
tr=@(p) center+scale*p;
body=[.27 .32 .36]; arm=[.41 .46 .50]; blade=[.24 .29 .33];
motors=[-27 -15;27 -15;-27 15;27 15];
for k=1:4
    q=motors(k,:);
    d=unit(q); n=[-d(2) d(1)];
    flat(ax,tr([4*d+2*n;q+2*n;q-2*n;4*d-2*n]),arm,1,'none',.5);
    flat(ax,tr(q+[-3 -2;3 -2;3 2;-3 2]),body,1,'none',.5);
end
% Faceted fuselage rather than a cartoon face.
flat(ax,tr([-7 -16;7 -16;12 0;7 15;-7 15;-12 0]),body,1,'none',.5);
flat(ax,tr([-5 -12;5 -12;7 -2;-7 -2]),[.60 .67 .71],1,'none',.5);
line2(ax,tr([-6 7;6 7]),[.57 .64 .68],.65);
% Two slim paddles per rotor; four distinct fixed orientations.
angles=[18 -22 -22 18];
paddle=[1 -1.1;8 -2.7;16 -2;17 -.1;7 1.4;1 1.1];
for k=1:4
    a=deg2rad(angles(k)); rot=[cos(a) -sin(a);sin(a) cos(a)];
    for sign=[-1 1]
        points=(sign*paddle)*rot.';
        % A slight vertical flattening keeps the same oblique icon style.
        points(:,2)=.78*points(:,2);
        flat(ax,tr(points+motors(k,:)),blade,1,'none',.5);
    end
end
end

function data = readPointingNpy(path)
% Minimal NPY reader for real float64 N-by-3 pointing-vector files.
fid=fopen(path,'r','ieee-le');
assert(fid>=0,'draw_fig1:NpyOpen','Cannot open %s.',path);
cleanup=onCleanup(@() fclose(fid));
magic=fread(fid,6,'*uint8').';
assert(isequal(magic,[147 uint8('NUMPY')]),'draw_fig1:NpyMagic','Invalid NPY file.');
ver=fread(fid,2,'*uint8');
if ver(1)==1, n=fread(fid,1,'uint16'); else, n=fread(fid,1,'uint32'); end
header=char(fread(fid,n,'*uint8').');
assert(contains(header,'''<f8''') || contains(header,'''|f8'''), ...
    'draw_fig1:NpyType','Expected little-endian float64 NPY.');
dims=regexp(header,'''shape''\s*:\s*\((\d+),\s*(\d+)','tokens','once');
assert(~isempty(dims),'draw_fig1:NpyShape','Expected two-dimensional NPY.');
rows=str2double(dims{1}); cols=str2double(dims{2});
assert(cols==3,'draw_fig1:NpyShape','Expected N-by-3 NPY.');
raw=fread(fid,rows*cols,'*double');
assert(numel(raw)==rows*cols,'draw_fig1:NpyLength','Incomplete NPY payload.');
if contains(header,'''fortran_order'': True')
    data=reshape(raw,rows,cols);
else
    data=reshape(raw,cols,rows).';
end
end




















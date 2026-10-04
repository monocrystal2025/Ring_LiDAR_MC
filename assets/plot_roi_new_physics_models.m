clear; close all; clc;
set(groot, 'defaultTextInterpreter', 'none');
set(groot, 'defaultAxesTickLabelInterpreter', 'none');
set(groot, 'defaultLegendInterpreter', 'none');

scriptDir = fileparts(mfilename('fullpath'));

R = 2000;
JIAODU = 25;
fasan_D = 50e-3;
fasan_d = 0.8e-3;
f = 5e3;
omiga = 2*pi;

drawConeAndMapping(scriptDir, R, JIAODU, fasan_D, f, omiga);
drawLineModel(scriptDir, JIAODU, fasan_D, fasan_d, f, omiga);
drawPointModel(scriptDir, JIAODU, fasan_D, f, omiga);

disp(fullfile(scriptDir, 'roi_new_cone_and_scan_mapping.png'));
disp(fullfile(scriptDir, 'roi_new_line_physical_model.png'));
disp(fullfile(scriptDir, 'roi_new_point_physical_model.png'));

function drawConeAndMapping(scriptDir, R, jiaodu, fasan_D, f, omiga)
alpha = deg2rad(jiaodu);
scanRadius = tan(alpha);
z = linspace(0, R, 160);
rho = z*tan(alpha);

fig = figure('Color','w','Position',[100 100 1500 650]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

nexttile;
hold on;
fill([z fliplr(z)], [rho -fliplr(rho)], [0.72 0.84 1.00], ...
    'FaceAlpha',0.25,'EdgeColor','none');
plot(z, rho, 'Color',[0.10 0.35 0.85], 'LineWidth',2);
plot(z, -rho, 'Color',[0.10 0.35 0.85], 'LineWidth',2);
plot([0 R], [0 0], 'k-', 'LineWidth',1.2);
scatter(0,0,55,'k','filled');
text(35,40,'LiDAR / cone apex','FontSize',10);

s0 = 1450;
p0 = [s0*cos(alpha), s0*sin(alpha)];
quiver(p0(1), p0(2), 260*0.94, -260*0.34, 0, ...
    'Color',[0.80 0.10 0.10], 'LineWidth',1.8, 'MaxHeadSize',0.7);
scatter(p0(1),p0(2),60,[0.80 0.10 0.10],'filled');
text(p0(1)+35,p0(2)-135,{'initial UAV on cone side', ...
    'slant range 100-2000 m', 'velocity points inward'},'FontSize',9);
text(780, R*tan(alpha)*0.80, 'ROI condition: range <= R, rho <= z*tan(JIAODU)', 'FontSize',10);
xlabel('z along cone axis (m)');
ylabel('radial distance rho (m)');
title('Cone ROI used by NEW line/point code');
axis equal; grid on; box off;
xlim([0 R*1.04]);
ylim([-R*tan(alpha)*1.05 R*tan(alpha)*1.05]);

nexttile;
hold on;
theta = linspace(0,2*pi,300);
plot(scanRadius*cos(theta), scanRadius*sin(theta), ...
    'Color',[0.10 0.35 0.85], 'LineWidth',2);
beamD = 2*tan(fasan_D/2);
stepSize = omiga/f;
raster = generateRasterPath(beamD, stepSize, scanRadius);
spiral = generateSpiralPath(scanRadius, beamD, stepSize);
plot(raster(:,1), raster(:,2), 'Color',[0.00 0.48 0.75], 'LineWidth',1.1);
plot(spiral(:,1), spiral(:,2), 'Color',[0.45 0.20 0.75], 'LineWidth',1.1);
q = [0.23 -0.13];
B = [q 1] ./ norm([q 1]);
scatter(q(1),q(2),55,'k','filled');
text(q(1)+0.02,q(2)-0.04,'q=[x,y] in tangent plane','FontSize',9);
text(-0.43,0.43,'beam direction B = normalize([x, y, 1])','FontSize',10);
text(-0.43,0.37,sprintf('example B=[%.2f, %.2f, %.2f]',B(1),B(2),B(3)),'FontSize',9);
xlabel('x/z'); ylabel('y/z');
title('Scan paths generated in tangent-plane coordinates');
legend({'cone boundary','raster path','spiral path'},'Location','northeast');
axis equal; grid on; box off;
xlim(scanRadius*[-1.13 1.13]); ylim(scanRadius*[-1.13 1.13]);

exportgraphics(fig, fullfile(scriptDir, 'roi_new_cone_and_scan_mapping.png'), 'Resolution', 180);
close(fig);
end

function drawLineModel(scriptDir, jiaodu, fasan_D, fasan_d, f, omiga)
scanRadius = tand(jiaodu);
beamD = 2*tan(fasan_D/2);
stepSize = omiga/f;
raster = generateRasterPath(beamD, stepSize, scanRadius);
spiral = generateSpiralPath(scanRadius, beamD, stepSize);
[eLong, eShort] = localFrameSpiralXY(spiral, beamD);

fig = figure('Color','w','Position',[90 80 1500 1050]);
tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

nexttile;
hold on;
drawBoundary(scanRadius);
plot(raster(:,1), raster(:,2), 'Color',[0.00 0.48 0.75], 'LineWidth',1.0);
q = raster(round(size(raster,1)/2),:);
drawRotRect(q, [1 0], [0 1], 0.18, 0.035, [0.80 0.10 0.10]);
quiver(q(1),q(2),0.14,0,0,'Color',[0.80 0.10 0.10],'LineWidth',1.7,'MaxHeadSize',0.7);
quiver(q(1),q(2),0,0.09,0,'Color',[0.05 0.55 0.20],'LineWidth',1.7,'MaxHeadSize',0.7);
text(q(1)+0.07,q(2)+0.02,'long','Color',[0.80 0.10 0.10],'FontWeight','bold');
text(q(1)+0.02,q(2)+0.075,'short','Color',[0.05 0.55 0.20],'FontWeight','bold');
text(-0.44,-0.47,'raster line frame in code: e_long=[1,0], e_short=[0,1]','FontSize',9);
title('Line beam, raster: long side fixed along x');
xlabel('x/z'); ylabel('y/z'); finishPlaneAxes(scanRadius);

nexttile;
hold on;
drawBoundary(scanRadius);
plot(spiral(:,1), spiral(:,2), 'Color',[0.45 0.20 0.75], 'LineWidth',1.2);
idx = round(size(spiral,1)*0.65);
q = spiral(idx,:);
el = eLong(idx,:);
es = eShort(idx,:);
drawRotRect(q, el, es, 0.32, 0.18, [0.80 0.10 0.10]);
quiver(q(1),q(2),0.16*el(1),0.16*el(2),0,'Color',[0.80 0.10 0.10],'LineWidth',1.7,'MaxHeadSize',0.7);
quiver(q(1),q(2),0.13*es(1),0.13*es(2),0,'Color',[0.05 0.55 0.20],'LineWidth',1.7,'MaxHeadSize',0.7);
text(q(1)+0.12*el(1),q(2)+0.12*el(2),'long=normal','Color',[0.80 0.10 0.10],'FontWeight','bold');
text(q(1)+0.12*es(1),q(2)+0.12*es(2),'short=tangent','Color',[0.05 0.55 0.20],'FontWeight','bold');
text(-0.44,-0.47,{'spiral line frame in code:', ...
    'e_short is analytic path tangent; e_long=[-e_short_y, e_short_x]'},'FontSize',9);
title('Line beam, spiral: legacy long/short physical orientation');
xlabel('x/z'); ylabel('y/z'); finishPlaneAxes(scanRadius);

nexttile;
long = linspace(-fasan_D, fasan_D, 360);
short = linspace(-4*fasan_d, 4*fasan_d, 360);
[LL,SS] = meshgrid(long, short);
gaussianWidth = fasan_d/2;
I = (abs(LL) <= fasan_D/2) .* exp(-2*(SS/gaussianWidth).^2);
imagesc(long*1e3, short*1e3, I);
set(gca,'YDir','normal');
colormap(gca,'hot'); colorbar;
xline(-fasan_D*500,'w--','LineWidth',1);
xline(fasan_D*500,'w--','LineWidth',1);
title('Line illumination used for SNR');
xlabel('long_offset (mrad), hard edge at +/- fasan_D/2');
ylabel('short_offset (mrad), Gaussian width fasan_d/2');

nexttile;
axis off;
text(0.02,0.98,{ ...
    'Line-beam detection equations from NEW code', ...
    '', ...
    'long_half_width = fasan_D/2', ...
    'gaussian_width = fasan_d/2', ...
    'candidate if:', ...
    '  |long_offset| <= long_half_width + target_half_theta', ...
    '  |short_offset| <= 3*fasan_d + target_half_theta', ...
    '', ...
    'illumination factor:', ...
    '  inside_long * exp(-2*(sample_short/gaussian_width)^2)', ...
    '', ...
    'SNR decision:', ...
    '  sliding window, threshold = 2', ...
    '  WINDOW_PULSES = ceil(f*fasan_D/omiga)'}, ...
    'VerticalAlignment','top','FontName','Consolas','FontSize',10);

exportgraphics(fig, fullfile(scriptDir, 'roi_new_line_physical_model.png'), 'Resolution', 180);
close(fig);
end

function drawPointModel(scriptDir, jiaodu, fasan_D, f, omiga)
scanRadius = tand(jiaodu);
beamD = 2*tan(fasan_D/2);
stepSize = omiga/f;
raster = generateRasterPath(beamD, stepSize, scanRadius);
spiral = generateSpiralPath(scanRadius, beamD, stepSize);

fig = figure('Color','w','Position',[90 80 1500 1050]);
tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

drawPointPathTile(raster, scanRadius, 'Point beam, raster: circular Gaussian spot follows path');
drawPointPathTile(spiral, scanRadius, 'Point beam, spiral: circular Gaussian spot follows path');

nexttile;
lim = 4*(fasan_D/2);
x = linspace(-lim, lim, 360);
y = linspace(-lim, lim, 360);
[X,Y] = meshgrid(x,y);
theta = hypot(X,Y);
fovHalf = 3*(fasan_D/2);
I = (theta <= fovHalf) .* exp(-2*(theta/(fasan_D/2)).^2);
imagesc(x*1e3, y*1e3, I);
set(gca,'YDir','normal');
axis equal tight;
hold on;
drawCircle(0,0,(fasan_D/2)*1e3,'LineStyle','--','Color','w','LineWidth',1.1);
drawCircle(0,0,fovHalf*1e3,'LineStyle','-','Color','w','LineWidth',1.0);
colormap(gca,'parula'); colorbar;
title('Point-beam illumination used for SNR');
xlabel('angular x offset (mrad)');
ylabel('angular y offset (mrad)');

nexttile;
axis off;
text(0.02,0.98,{ ...
    'Point-beam detection equations from NEW code', ...
    '', ...
    'GAUSSIAN_HALF_ANGLE = fasan_D/2', ...
    'FOV_HALF_ANGLE = 3*GAUSSIAN_HALF_ANGLE', ...
    '', ...
    'candidate if target square intersects FOV:', ...
    '  cos_theta >= cos(FOV_HALF_ANGLE + target_half_diag_theta)', ...
    '', ...
    'illumination factor:', ...
    '  inside_fov * exp(-2*(sample_theta/GAUSSIAN_HALF_ANGLE)^2)', ...
    '', ...
    'SNR decision:', ...
    '  sliding window, threshold = 2', ...
    '  WINDOW_PULSES = ceil(f*(3*fasan_D)/omiga)'}, ...
    'VerticalAlignment','top','FontName','Consolas','FontSize',10);

exportgraphics(fig, fullfile(scriptDir, 'roi_new_point_physical_model.png'), 'Resolution', 180);
close(fig);
end

function drawPointPathTile(path, scanRadius, titleText)
nexttile;
hold on;
drawBoundary(scanRadius);
plot(path(:,1), path(:,2), 'Color',[0.35 0.20 0.75], 'LineWidth',1.1);
q = path(round(size(path,1)/2),:);
drawCircle(q(1), q(2), 0.08, 'LineStyle','-','Color','r','LineWidth',1.8);
drawCircle(q(1), q(2), 0.15, 'Color',[0.95 0.45 0.05], 'LineStyle','--', 'LineWidth',1.3);
scatter(q(1), q(2), 35, 'k', 'filled');
text(q(1)+0.035,q(2)+0.075,'Gaussian center B=normalize([x,y,1])','FontSize',8);
text(-0.44,-0.47,'solid: Gaussian half-angle fasan_D/2; dashed: FOV half-angle 3*fasan_D/2','FontSize',8);
title(titleText);
xlabel('x/z'); ylabel('y/z'); finishPlaneAxes(scanRadius);
end

function path = generateSpiralPath(rEnd, pitch, stepSize)
a = pitch/(2*pi);
thetaMax = rEnd/a;
L = 0.5*a*(thetaMax*sqrt(thetaMax^2 + 1) + log(thetaMax + sqrt(thetaMax^2 + 1)));
nPoints = ceil(L/stepSize) + 1;
path = zeros(nPoints,2);
path(1,:) = [0 0];
pointCount = 1;
theta = 0;
while true
    dsDtheta = a*sqrt(theta^2 + 1);
    thetaNew = theta + stepSize/dsDtheta;
    rNew = a*thetaNew;
    if rNew >= rEnd
        thetaEnd = rEnd/a;
        pointCount = pointCount + 1;
        path(pointCount,:) = [rEnd*cos(thetaEnd), rEnd*sin(thetaEnd)];
        break
    end
    pointCount = pointCount + 1;
    path(pointCount,:) = [rNew*cos(thetaNew), rNew*sin(thetaNew)];
    theta = thetaNew;
end
path = path(1:pointCount,:);
end

function path = generateRasterPath(D, d, a)
if D >= 2*a
    y = (-a-d:d:a+d).';
    x = zeros(size(y));
    path = [x y];
    return
end
lineCount = max(1, ceil(2*a/D));
path = zeros(0,2);
for i = 1:lineCount
    x0 = -a + D/2 + D*(i-1);
    ySamples = linspace(-a,a,max(2,ceil(2*a/d))).';
    if mod(i,2) == 0
        ySamples = flipud(ySamples);
    end
    xSamples = x0*ones(size(ySamples));
    path = [path; xSamples ySamples]; %#ok<AGROW>
    if i ~= lineCount
        xNext = -a + D/2 + D*i;
        connectorCount = max(2,ceil(abs(xNext-x0)/d));
        xConn = linspace(x0,xNext,connectorCount).';
        yConn = ySamples(end)*ones(size(xConn));
        path = [path; xConn(2:end) yConn(2:end)]; %#ok<AGROW>
    end
end
inside = hypot(path(:,1),path(:,2)) <= a + eps(a);
path = path(inside,:);
if isempty(path)
    path = [0 0];
end
end

function [eLong, eShort] = localFrameSpiralXY(P, pitch)
x = P(:,1);
y = P(:,2);
r = hypot(x,y);
b = pitch/(2*pi);
tx = zeros(size(x));
ty = zeros(size(y));
mask = r > eps;
tx(mask) = (b./r(mask)).*x(mask) - y(mask);
ty(mask) = (b./r(mask)).*y(mask) + x(mask);
tx(~mask) = 1;
ty(~mask) = 0;
nt = hypot(tx,ty);
eShort = [tx./nt ty./nt];
eLong = [-eShort(:,2) eShort(:,1)];
end

function drawBoundary(scanRadius)
theta = linspace(0,2*pi,300);
plot(scanRadius*cos(theta), scanRadius*sin(theta), ...
    'Color',[0.10 0.35 0.85], 'LineWidth',1.7);
end

function finishPlaneAxes(scanRadius)
axis equal; grid on; box off;
xlim(scanRadius*[-1.13 1.13]);
ylim(scanRadius*[-1.13 1.13]);
end

function drawRotRect(center, eLong, eShort, longLen, shortLen, color)
corners = [ ...
    center - longLen/2*eLong - shortLen/2*eShort; ...
    center + longLen/2*eLong - shortLen/2*eShort; ...
    center + longLen/2*eLong + shortLen/2*eShort; ...
    center - longLen/2*eLong + shortLen/2*eShort];
patch(corners(:,1), corners(:,2), color, 'FaceColor','none', ...
    'EdgeColor',color, 'LineWidth',2);
end

function drawCircle(x0,y0,r,varargin)
theta = linspace(0,2*pi,200);
plot(x0 + r*cos(theta), y0 + r*sin(theta), varargin{:});
end

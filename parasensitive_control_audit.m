function control = parasensitive_control_audit()
%PARASENSITIVE_CONTROL_AUDIT Paired moving-target controls for baseline changes.
% Reuses all 64 saved trajectories and phases; no static simulations.
root = fileparts(mfilename('fullpath'));
loaded = load(fullfile(root,'parasensitive_data.mat'),'data');
d = loaded.data;
pop = d.population;
settings = [150e-6 6000 5e-5;200e-6 5000 5e-5; ...
    200e-6 6000 1.5e-5;150e-6 5000 1.5e-5];
names = ["Energy only: 150 uJ","PRF only: 5 kHz", ...
    "Alpha only: 1.5e-5","Paper-entry parameters"];
widths = d.config.WidthsMrad;
n = d.config.N;
times = nan(n,numel(widths),2,4);
detect = false(size(times));
if isempty(gcp('nocreate'))
    parpool('threads');
end
timer = tic;
for c = 1:4
    f = settings(c,2);
    physics = [settings(c,1),0.5,settings(c,3),settings(c,3)/50];
    for k = 1:numel(widths)
        scan = makeScan(widths(k)*1e-3,f,d.config.omega);
        idx = floor(pop.phaseFraction*scan.phaseLength)+1;
        for b = 1:2
            if b==1
                kernel = @MC_line_snr;
            else
                kernel = @MC_ring_snr;
            end
            hit = false(n,1);
            pulse = nan(n,1);
            P = pop.entryPositions;
            V = pop.unitVelocity*30;
            width = widths(k)*1e-3;
            block = mc_ea_block_size(width,char(d.config.beams(b)));
            parfor i = 1:n
                [hit(i),pulse(i)] = kernel(width,1e-3,scan.path1,scan.beam, ...
                    scan.path2,f,2*pi,P(i,:),V(i,:),idx(i),2000,block,Inf,'fixed',physics);
            end
            detect(:,k,b,c) = hit;
            times(:,k,b,c) = pulse/f;
        end
    end
    fprintf('Completed paired control %s in %.1f s\n',names(c),toc(timer));
end
requirements = [15 .6;30 .7;60 .9;19 .88];
minimumWidths = nan(4,2,4);
reduction = nan(4,4);
for c = 1:4
    for r = 1:4
        for b = 1:2
            t = requirements(r,1);
            p = cummax(mean(detect(:,:,b,c) & isfinite(times(:,:,b,c)) ...
                & times(:,:,b,c)>=0 & times(:,:,b,c)<t+eps(t),1));
            k = find(p>=requirements(r,2),1);
            if ~isempty(k)
                if k==1
                    minimumWidths(r,b,c) = widths(1);
                elseif abs(diff(p(k-1:k)))<=eps(max(abs(p(k-1:k))))
                    minimumWidths(r,b,c) = widths(k);
                else
                    minimumWidths(r,b,c) = interp1(p(k-1:k),widths(k-1:k),requirements(r,2));
                end
            end
        end
        reduction(r,c) = 100*(1-minimumWidths(r,2,c)/minimumWidths(r,1,c));
    end
end
control = struct('settings',settings,'names',names,'requirements',requirements, ...
    'minimumWidths',minimumWidths,'reductionPct',reduction,'detect',detect, ...
    'firstTime',times,'widthsMrad',widths,'N',n);
save(fullfile(root,'parasensitive_control_audit.mat'),'control');
disp(names);
disp([requirements reduction]);
end
function scan = makeScan(width,f,omega)
% Same spherical spiral and ordering as flipud(lidarpath(...)), without UI
% or external NPY files. Regenerate on the actual 1/f pulse-time grid.
spacing = 2*tan(width/2);
endpoint = min(1,spacing*sqrt(1-spacing^2/4));
phiLower = acos(endpoint)+0.5*asin(endpoint);
phiUpper = asin(endpoint)/2;
k = spacing/(2*pi);
arc = integral(@(phi) sqrt(k^2+sin(phi).^2),phiUpper,phiLower, ...
    'RelTol',1e-11,'AbsTol',1e-12)/k;
pulseTimes = (0:floor(arc/omega*f)).'/f;
solution = ode45(@(~,phi) -k./sqrt(k^2+sin(phi).^2),[0 arc], ...
    phiLower,odeset('RelTol',1e-9,'AbsTol',1e-11));
phi = deval(solution,omega*pulseTimes).';
phi = min(phiLower,max(phiUpper,phi));
theta = (phiLower-phi)/k;
beam = [sin(phi).*cos(theta),sin(phi).*sin(theta),cos(phi)];
scan.beam = flipud(beam./vecnorm(beam,2,2));
scan.path1 = boundaryCircle(scan.beam(1,:),f,omega);
scan.path2 = boundaryCircle(scan.beam(end,:),f,omega);
scan.phaseLength = size(scan.path1,1)+size(scan.beam,1)+size(scan.path2,1);
scan.cycleLength = scan.phaseLength+size(scan.beam,1);
scan.f = f;
end

function circle = boundaryCircle(point,f,omega)
theta0 = atan2(point(2),point(1));
radius = hypot(point(1),point(2));
theta = (theta0:omega/f/radius:theta0+2*pi).';
circle = [radius*cos(theta),radius*sin(theta),repmat(point(3),numel(theta),1)];
end


function parasensitive_add_clusters()
%PARASENSITIVE_ADD_CLUSTERS Replay successful saved trials to recover clusters.
% Counts consecutive positive-signal groups exactly as advantage_data.m.
% Original samples, detections and pulse counts must match before saving.
root = fileparts(mfilename('fullpath'));
file = fullfile(root,'parasensitive_data.mat');
s = load(file,'data');
data = s.data;
cfg = data.config;
assert(data.complete,'Complete original simulation first.');
for k = 2:numel(cfg.sourceNames)
    assert(strcmp(fileread(fullfile(root,cfg.sourceNames(k))),cfg.sourceText{k}), ...
        'Simulation dependency changed: %s',cfg.sourceNames(k));
end
if isempty(gcp('nocreate'))
    parpool('threads');
end
pop = data.population;
timer = tic;
for d = numel(cfg.WidthsMrad):-1:1
    width = cfg.WidthsMrad(d)*1e-3;
    for f = unique(data.cases(:,2)).'
        members = find(data.cases(:,2)==f).';
        scan = makeScan(width,f,cfg.omega);
        idx = floor(pop.phaseFraction*scan.phaseLength)+1;
        for j = members
            setting = data.cases(j,:);
            physics = [setting(1),0.5,setting(4),setting(4)/cfg.lidarRatio];
            for b = 1:2
                raw = data.dynamic{d,j,b};
                if isfield(raw,'effectiveClusterCount')
                    continue
                end
                selected = find(raw.detect);
                P = pop.entryPositions(selected,:);
                V = pop.unitVelocity(selected,:)*setting(3);
                starts = idx(selected);
                if b==1
                    kernel = @MC_line_snr;
                else
                    kernel = @MC_ring_snr;
                end
                block = mc_ea_block_size(width,char(cfg.beams(b)));
                clusters = nan(numel(selected),1);
                counts = clusters;
                firstTime = clusters;
                hit = false(numel(selected),1);
                parfor i = 1:numel(selected)
                    [hit(i),pulse,~,window] = kernel(width,cfg.smallWidthRad, ...
                        scan.path1,scan.beam,scan.path2,f,cfg.omega, ...
                        P(i,:),V(i,:),starts(i),cfg.R,block,Inf,'fixed',physics);
                    positive = isfinite(window.signal_photons(:)) & window.signal_photons(:)>0;
                    counts(i) = nnz(positive);
                    clusters(i) = nnz(positive & [true;~positive(1:end-1)]);
                    firstTime(i) = pulse/f;
                end
                assert(all(hit) && isequal(counts,raw.effectivePulseCount(selected)) ...
                    && isequal(firstTime,raw.firstTime(selected)), ...
                    'Replay mismatch: width=%g, case=%d, beam=%d',width,j,b);
                raw.effectiveClusterCount = nan(cfg.N,1);
                raw.effectiveClusterCount(selected) = clusters;
                data.dynamic{d,j,b} = raw;
            end
        end
    end
    save(file,'data','-v7.3');
    fprintf('Cluster checkpoint D=%g mrad (%.1f s).\n',width*1e3,toc(timer));
end
data.clusterMetadata = struct('method',"consecutive positive signal samples", ...
    'replayVerified',true,'seconds',toc(timer));
save(file,'data','-v7.3');
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


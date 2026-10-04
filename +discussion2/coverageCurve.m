function output = coverageCurve(cfg,bank)
%COVERAGECURVE Geometric blind fraction on uniform hemispherical directions.
arguments
    cfg (1,1) struct
    bank (1,1) struct
end
u = bank.u;
az = 2*pi*u(:,1);
cz = u(:,2);
directions = [sqrt(1-cz.^2).*cos(az),sqrt(1-cz.^2).*sin(az),cz];
widths = cfg.widthGrid_mrad;
covered = false(bank.N,numel(widths),2);
source = [string(fileread(fullfile(cfg.root,"+discussion2/coverageMask.m"))), ...
    string(fileread(fullfile(cfg.root,"+discussion/scan.m")))];
folder = fullfile(cfg.cacheRoot,"coverage");
if ~isfolder(folder)
    mkdir(folder);
end
for j = 1:numel(widths)
    inputs = struct('version',cfg.version,'directions',directions, ...
        'width_mrad',widths(j),'thickness',cfg.system.beamWidth_rad, ...
        'angularStep',cfg.system.scanRate_rad_s/cfg.system.prf_Hz, ...
        'scanSource',cfg.scanSource,'scanLibrary',cfg.scanLibrary,'source',source);
    filename = fullfile(folder,discussion2.hash(inputs)+".mat");
    if isfile(filename)
        saved = load(filename,'hit','inputs');
        if isequaln(saved.inputs,inputs)
            covered(:,j,:) = reshape(saved.hit,bank.N,1,2);
            continue
        end
    end
    scan = discussion.scan(cfg,widths(j)*1e-3);
    hit = false(bank.N,2);
    for k = 1:2
        names = ["line","ring"];
        beam = discussion.beam(names(k),widths(j),cfg.system.beamWidth_rad);
        hit(:,k) = discussion2.coverageMask(directions,scan,beam);
    end
    covered(:,j,:) = reshape(hit,bank.N,1,2);
    save(filename,'inputs','hit','-v7');
    fprintf('Geometric coverage wD=%g step=%g: L %.4f R %.4f\n',widths(j), ...
        inputs.angularStep,mean(hit(:,1)),mean(hit(:,2)));
end
output = struct('widths_mrad',widths,'covered',covered,'N',bank.N, ...
    'blind',1-reshape(mean(covered,1),numel(widths),2));
end

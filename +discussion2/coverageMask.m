function covered = coverageMask(directions,scan,beam)
%COVERAGEMASK Exact point-centre 1/e^2 footprint, as static_MC_EA_NEW.m.
% This purely geometric test has no SNR, target-size, or range dependence.
arguments
    directions (:,3) double
    scan (1,1) struct
    beam (1,1) struct
end
covered = false(size(directions,1),1);
beamBlock = 10000;
targetBlock = 128;
for first = 1:targetBlock:size(directions,1)
    indices = first:min(size(directions,1),first+targetBlock-1);
    targets = directions(indices,:);
    hit = false(numel(indices),1);
    for start = 1:beamBlock:size(scan.B,1)
        active = find(~hit);
        if isempty(active)
            break
        end
        rows = start:min(size(scan.B,1),start+beamBlock-1);
        along = targets(active,:)*scan.B(rows,:).';
        if beam.name == "ring"
            outer = beam.wD_rad/2+beam.wd_rad/2;
            inner = max(0,beam.wD_rad/2-beam.wd_rad/2);
            inside = along>=cos(outer) & along<=cos(inner);
        else
            long = targets(active,:)*scan.long(rows,:).';
            short = targets(active,:)*scan.short(rows,:).';
            inside = along>0 & abs(long)<=along*tan(beam.wD_rad/2) & ...
                abs(short)<=along*tan(beam.wd_rad/2);
        end
        hit(active(any(inside,2))) = true;
    end
    covered(indices) = hit;
end
end

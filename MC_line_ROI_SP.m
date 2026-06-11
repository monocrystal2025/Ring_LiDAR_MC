function [detect,first_detect]=MC_line_ROI_SP(fasan_D,fasan_d,beam_vec,path2,f,omiga,UP,UV,init_beam_idx,z,a)
%MC_line_ROI_SP
% Spiral-path MC detection for a line beam.
% The line-beam long side is perpendicular to the local path tangent.
%
% Compared with the previous version, the local tangent is NOT estimated by
% neighbor differences. Instead, it is computed analytically from the point
% coordinates with respect to the origin:
%   1) beam_vec / flipud(beam_vec): Archimedean spiral tangent
%   2) path2                  : circular tangent
%
% Inputs are kept compatible with MC_ring_ROI_SP.

    %#ok<*NASGU>
    detect = 0;
    first_detect = nan;
    N_KERNAL = 100000;

    D = 2*z*tan(fasan_D/2);   % long-side footprint size
    d = 2*z*tan(fasan_d/2);   % short-side footprint size

    % One-period precomputation -------------------------------------------------
    % beam_vec is the forward Archimedean spiral.
    [eL_beam, eS_beam] = local_frame_spiral_xy(beam_vec, D);

    if isempty(path2)
        P_cycle  = [beam_vec; flipud(beam_vec)];
        eL_cycle = [eL_beam; flipud(eL_beam)];
        eS_cycle = [eS_beam; flipud(eS_beam)];
    else
        % path2 is the circular segment.
        [eL_path2, eS_path2] = local_frame_circle_xy(path2);

        P_cycle  = [beam_vec; path2; flipud(beam_vec)];
        eL_cycle = [eL_beam; eL_path2; flipud(eL_beam)];
        eS_cycle = [eS_beam; eS_path2; flipud(eS_beam)];
    end

    % First pulse ---------------------------------------------------------------
    L = size(P_cycle,1);
    idx0 = mod(init_beam_idx-1, L) + 1;

    B0  = P_cycle(idx0,:);
    eL0 = eL_cycle(idx0,:);
    eS0 = eS_cycle(idx0,:);

    rel0 = UP - B0;
    u0 = rel0 * eL0.';
    v0 = rel0 * eS0.';
    if (abs(u0) <= D/2) && (abs(v0) <= d/2)
        detect = 1;
        first_detect = 0;
        return
    end

    % Batch propagation ---------------------------------------------------------
    P_now = [ones(N_KERNAL,1)*UP(1), ones(N_KERNAL,1)*UP(2)];
    dt = (0:1/f:N_KERNAL*1/f).';
    dt(1) = [];
    P_now = P_now + [UV(1)*dt, UV(2)*dt];

    i = 1;
    beam_idx = init_beam_idx + 1;

    B_now  = cut_loop_simple_anyk(P_cycle,  beam_idx, N_KERNAL);
    EL_now = cut_loop_simple_anyk(eL_cycle, beam_idx, N_KERNAL);
    ES_now = cut_loop_simple_anyk(eS_cycle, beam_idx, N_KERNAL);

    while true
        norm_P_now = vecnorm(P_now, 2, 2);

        rel = P_now - B_now;
        u = sum(rel .* EL_now, 2);
        v = sum(rel .* ES_now, 2);
        judge = (abs(u) <= D/2) & (abs(v) <= d/2);
        judge(norm_P_now > a) = 0;

        if any(judge)
            detect = 1;
            first_detect = (i-1)*N_KERNAL + find(judge,1,'first');
            return
        end

        norm_P_now(norm_P_now > a) = nan;
        i = i + 1;
        if anynan(norm_P_now)
            break
        end

        P_now(:,1) = P_now(:,1) + N_KERNAL*UV(1)/f;
        P_now(:,2) = P_now(:,2) + N_KERNAL*UV(2)/f;

        beam_idx = beam_idx + N_KERNAL;
        B_now  = cut_loop_simple_anyk(P_cycle,  beam_idx, N_KERNAL);
        EL_now = cut_loop_simple_anyk(eL_cycle, beam_idx, N_KERNAL);
        ES_now = cut_loop_simple_anyk(eS_cycle, beam_idx, N_KERNAL);
    end
end


function [e_long, e_short] = local_frame_spiral_xy(P, pitch)
% For Archimedean spiral r = b*theta, b = pitch/(2*pi).
% Tangent in polar basis:
%   t = dr/dtheta * e_r + r * e_theta = b*e_r + r*e_theta
% In Cartesian with point (x,y), r = hypot(x,y):
%   t = (b/r) * [x, y] + [-y, x]
%
% e_short : along tangent
% e_long  : perpendicular to tangent

    x = P(:,1);
    y = P(:,2);
    r = hypot(x,y);
    b = pitch/(2*pi);

    tx = zeros(size(x));
    ty = zeros(size(y));

    mask = (r > eps);
    tx(mask) = (b ./ r(mask)) .* x(mask) - y(mask);
    ty(mask) = (b ./ r(mask)) .* y(mask) + x(mask);

    % Limit at r -> 0 for Archimedean spiral: tangent points along +x.
    tx(~mask) = 1;
    ty(~mask) = 0;

    nt = hypot(tx,ty);
    e_short = [tx ./ nt, ty ./ nt];
    e_long  = [-e_short(:,2), e_short(:,1)];
end


function [e_long, e_short] = local_frame_circle_xy(P)
% For circular path, tangent is perpendicular to the radial direction.
% Using CCW tangent:
%   t = [-y, x]
% For the rectangle test, +/- tangent gives the same orientation, so the
% sign is immaterial here.

    x = P(:,1);
    y = P(:,2);
    r = hypot(x,y);

    tx = -y;
    ty =  x;

    mask = (r > eps);
    tx(~mask) = 0;
    ty(~mask) = 1;

    nt = hypot(tx,ty);
    e_short = [tx ./ nt, ty ./ nt];
    e_long  = [-e_short(:,2), e_short(:,1)];
end


function seg = cut_loop_simple_anyk(P,k,m)
% k can be any positive absolute index in the infinite repeated sequence.
% Returns m consecutive rows starting from index k with loop wrap-around.

    if m <= 0
        seg = zeros(0,size(P,2));
        return;
    end

    L = size(P,1);
    k0 = mod(k-1, L) + 1;
    idx = mod((k0-1) + (0:m-1), L) + 1;
    seg = P(idx, :);
end

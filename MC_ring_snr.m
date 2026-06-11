function [detect, first_detect] = MC_ring(fasan_D, fasan_d, path1, beam_vec, path2, f, omiga, UP, UV, init_beam_idx, R)
%MC_RING  Monte Carlo detection judgment for annular beam.
%
% The geometric hit condition is unchanged: the target is still simplified as
% a point, and a hit is judged only by the angular relation between the target
% direction and the annular footprint.
%
% New detection rule:
%   1) SNR >= 2.80: one hit is accepted as successful detection.
%   2) 1.64 <= SNR < 2.80: the hit is stored as a candidate. Detection is
%      accepted only if another hit with SNR >= 1.64 occurs before the LiDAR
%      pointing direction moves farther than 2*wD from the candidate hit.
%   3) SNR < 1.64: the hit is ignored.
%
% first_detect is the pulse index at which the warning is actually generated.
% Therefore, for a weak two-hit confirmation it is the second-hit index.

detect = 0;
first_detect = nan;
N_KERNAL = 100000;

SNR_WEAK = 1.64;
SNR_STRONG = 2.80;
CONFIRM_ANGLE = 2 * fasan_D;      % angular interval between two hit events

pending_valid = false;
pending_beam = [nan, nan, nan];

%% One complete beam cycle used by the original code
Beam_all = [beam_vec; path2; flipud(beam_vec); path1];
L_cycle = size(Beam_all, 1);

%% Check the initial pulse, which corresponds to first_detect = 0
idx0 = mod(init_beam_idx - 1, L_cycle) + 1;
B0 = Beam_all(idx0, :);
UP_norm = norm(UP);
UP_N = UP ./ UP_norm;

dtheta0 = safe_acos(dot(UP_N, unitvec(B0)));
if ((fasan_D/2 - fasan_d) <= dtheta0) && (dtheta0 <= fasan_D/2)
    snr0 = MC_single_pulse_snr('ring', fasan_D, fasan_d, UP_norm);
    if snr0 >= SNR_STRONG
        detect = 1;
        first_detect = 0;
        return;
    elseif snr0 >= SNR_WEAK
        pending_valid = true;
        pending_beam = B0;
    end
end

%% Precompute the first block of target positions and beam pointings
P_now = [ones(N_KERNAL,1)*UP(1), ones(N_KERNAL,1)*UP(2), ones(N_KERNAL,1)*UP(3)];
dt = (1:N_KERNAL)' ./ f;
dp = [UV(1)*dt, UV(2)*dt, UV(3)*dt];
P_now = P_now + dp;

B_now = cut_loop_simple_anyk(path1, beam_vec, path2, init_beam_idx + 1, N_KERNAL);
beam_idx = init_beam_idx + 1;
block_id = 1;

while true
    norm_P_now = vecnorm(P_now, 2, 2);

    dtheta = safe_acos_vec(dot(P_now, B_now, 2) ./ norm_P_now);
    judge = (((fasan_D/2 - fasan_d) <= dtheta) & (dtheta <= fasan_D/2));

    % Keep the original range-exit logic, but prevent invalid out-of-region
    % samples from being accepted inside the current vectorized block.
    judge(norm_P_now > R) = false;
    judge(P_now(:,3) < 0) = false;

    hit_idx = find(judge);
    if ~isempty(hit_idx)
        z_hit = norm_P_now(hit_idx);
        snr_hit = MC_single_pulse_snr('ring', fasan_D, fasan_d, z_hit);

        for kk = 1:numel(hit_idx)
            this_snr = snr_hit(kk);
            this_beam = B_now(hit_idx(kk), :);
            this_step = (block_id - 1) * N_KERNAL + hit_idx(kk);

            if this_snr >= SNR_STRONG
                detect = 1;
                first_detect = this_step;
                return;
            elseif this_snr >= SNR_WEAK
                if pending_valid
                    beam_sep = angular_distance(pending_beam, this_beam);
                    if beam_sep <= CONFIRM_ANGLE
                        detect = 1;
                        first_detect = this_step;
                        return;
                    else
                        % The old weak hit is outside the allowed short-time
                        % angular window. Discard it and start a new candidate.
                        pending_beam = this_beam;
                    end
                else
                    pending_valid = true;
                    pending_beam = this_beam;
                end
            end
        end
    end

    % Original termination condition: target has left the detection region.
    norm_P_tmp = norm_P_now;
    norm_P_tmp(norm_P_tmp > R) = nan;
    if any(isnan(norm_P_tmp)) || any(P_now(:,3) < 0)
        break;
    end

    block_id = block_id + 1;
    P_now(:,1) = P_now(:,1) + N_KERNAL * UV(1) / f;
    P_now(:,2) = P_now(:,2) + N_KERNAL * UV(2) / f;
    P_now(:,3) = P_now(:,3) + N_KERNAL * UV(3) / f;

    beam_idx = beam_idx + N_KERNAL;
    B_now = cut_loop_simple_anyk(path1, beam_vec, path2, beam_idx, N_KERNAL);
end
end

function seg = cut_loop_simple_anyk(a1, a2, a3, k, m)
% k: arbitrary start index in the infinite periodic beam sequence.
% seg: m-by-3 segment starting from k.

if m <= 0
    seg = zeros(0, 3);
    return;
end

P = [a2; a3; flipud(a2); a1];
L = size(P, 1);
k0 = mod(k - 1, L) + 1;
idx = mod((k0 - 1) + (0:m-1), L) + 1;
seg = P(idx, :);
end

function u = unitvec(v)
nv = norm(v);
if nv == 0
    u = v;
else
    u = v ./ nv;
end
end

function a = angular_distance(v1, v2)
a = safe_acos(dot(unitvec(v1), unitvec(v2)));
end

function y = safe_acos(x)
y = acos(max(min(x, 1), -1));
end

function y = safe_acos_vec(x)
y = acos(max(min(x, 1), -1));
end

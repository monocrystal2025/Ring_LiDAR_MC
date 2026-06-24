function [detect, first_detect] = MC_ring_snr(fasan_D, fasan_d, path1, beam_vec, path2, f, ~, UP, UV, init_beam_idx, R)
%MC_RING  Monte Carlo detection judgment for annular beam.
%
% The annular beam uses a radial Gaussian profile centered at fasan_D/2.
% fasan_d is the radial 1/e^2 angular width. The target is no longer treated
% as a point for the contact judgment: the code averages the normalized
% Gaussian illumination over the finite target area before calculating SNR.
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
TARGET_W = 0.297;                 % m, kept consistent with MC_single_pulse_snr
GAUSSIAN_CUTOFF_WIDTH = 3;         % ignore tails below exp(-18) before SNR calc
RING_RADIUS = fasan_D / 2;

if fasan_d <= 0
    error('MC_ring_snr:InvalidGaussianWidth', ...
        'fasan_d must be positive for the annular Gaussian profile.');
end

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
if is_ring_gaussian_candidate(dtheta0, UP_norm, TARGET_W, ...
        RING_RADIUS, fasan_d, GAUSSIAN_CUTOFF_WIDTH)
    illum0 = ring_gaussian_target_factor(dtheta0, UP_norm, ...
        TARGET_W, RING_RADIUS, fasan_d);
    snr0 = MC_single_pulse_snr('ring', fasan_D, fasan_d, UP_norm, ...
        illum0, 'gaussian');
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
    judge = is_ring_gaussian_candidate(dtheta, norm_P_now, TARGET_W, ...
        RING_RADIUS, fasan_d, GAUSSIAN_CUTOFF_WIDTH);

    % Keep the original range-exit logic, but prevent invalid out-of-region
    % samples from being accepted inside the current vectorized block.
    judge(norm_P_now > R) = false;
    judge(P_now(:,3) < 0) = false;

    hit_idx = find(judge);
    if ~isempty(hit_idx)
        z_hit = norm_P_now(hit_idx);
        illum_hit = ring_gaussian_target_factor(dtheta(hit_idx), z_hit, ...
            TARGET_W, RING_RADIUS, fasan_d);
        snr_hit = MC_single_pulse_snr('ring', fasan_D, fasan_d, z_hit, ...
            illum_hit, 'gaussian');

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

function candidate = is_ring_gaussian_candidate(center_theta, z, target_w, ring_radius, gaussian_width, cutoff_width)
% Fast geometric prefilter. The exact SNR still comes from target-area
% averaging; this only skips Gaussian tails that are numerically negligible.

z = z(:);
center_theta = center_theta(:);
target_half_diag_theta = atan((sqrt(2) * target_w / 2) ./ z);
candidate = abs(center_theta - ring_radius) <= ...
    (target_half_diag_theta + cutoff_width * gaussian_width);
candidate = candidate & isfinite(center_theta) & isfinite(z) & z > 0;
end

function factor = ring_gaussian_target_factor(center_theta, z, target_w, ring_radius, gaussian_width)
% Average normalized annular-Gaussian intensity over a square target.
%
% The target side length is target_w, matching the area used in the SNR
% model. The square is sampled in the local angular tangent plane; because
% the annular beam is radially symmetric, the square's in-plane orientation
% has only a second-order effect for these small target angles.

persistent sample_x sample_y
if isempty(sample_x)
    grid_n = 7;
    nodes = ((1:grid_n) - 0.5) ./ grid_n - 0.5;
    [xx, yy] = meshgrid(nodes, nodes);
    sample_x = xx(:).';
    sample_y = yy(:).';
end

center_theta = center_theta(:);
z = z(:);
factor = zeros(size(center_theta));

valid = isfinite(center_theta) & isfinite(z) & z > 0 & gaussian_width > 0;
if ~any(valid)
    return;
end

angular_side = 2 .* atan((target_w / 2) ./ z(valid));
dx = angular_side .* sample_x;
dy = angular_side .* sample_y;
sample_theta = sqrt((center_theta(valid) + dx).^2 + dy.^2);
sample_intensity = exp(-2 .* ((sample_theta - ring_radius) ./ gaussian_width).^2);
factor(valid) = mean(sample_intensity, 2);
end

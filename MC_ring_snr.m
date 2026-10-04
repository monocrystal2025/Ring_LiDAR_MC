function [detect, first_detect, first_encounter, pulse_score, ...
    post_detection_photons] = MC_ring_snr( ...
    fasan_D, fasan_d, path1, beam_vec, path2, f, omiga, ...
    UP, UV, init_beam_idx, R, N_KERNAL, max_step_index, window_mode, physics)
%MC_RING  Monte Carlo detection judgment for annular beam.
%
% The annular beam uses a radial Gaussian profile centered at fasan_D/2.
% fasan_d is the full radial 1/e^2 angular width. The target is no longer
% treated as a point for the contact judgment: the code averages the
% normalized Gaussian illumination over the finite target area before
% calculating SNR.
%
% New detection rule:
%   For every emitted pulse, compute N_sig and N_noise at the target range.
%   Missed pulses have N_sig = 0 and still contribute N_noise.  A causal
%   sliding window with M = ceil(f * fasan_D / omiga) pulses is evaluated as
%       D = sqrt(sum(N_sig.^2 ./ (N_sig + N_noise))).
%   Detection is accepted when D >= 2, and first_detect is the tail pulse
%   index of the first successful window.
%
% first_detect is the pulse index at which the warning is actually generated.
% first_encounter is the pulse-index difference between the first pulse that
% passes the geometric prefilter and entry into the detection region. The
% current Monte Carlo trajectories enter at pulse 0.
% In variable mode, pulse_score stores the shortest successful causal
% window. In fixed mode, it stores the single startup-or-M-length window.
% When requested, post_detection_photons stores the following pulse sequence.

if nargin < 12 || isempty(N_KERNAL)
    N_KERNAL = 25000;
end
validateattributes(N_KERNAL, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'positive'}, ...
    mfilename, 'N_KERNAL');
if nargin < 13 || isempty(max_step_index)
    max_step_index = inf;
elseif isfinite(max_step_index)
    validateattributes(max_step_index, {'numeric'}, ...
        {'scalar', 'real', 'integer', 'nonnegative'}, ...
        mfilename, 'max_step_index');
end
if nargin < 14 || isempty(window_mode)
    window_mode = 'variable';
end
window_mode = validatestring(window_mode, {'variable', 'fixed'}, ...
    mfilename, 'window_mode');
collect_post_detection = nargout >= 5;
% Optional [pulse energy (J), reflectivity, alpha (1/m), beta (1/(m sr))].
if nargin < 15 || isempty(physics)
    physics = [150e-6, 0.5, 1.5e-5, 0.3e-6];
end
validateattributes(physics, {'double'}, ...
    {'vector','numel',4,'real','finite','positive'}, mfilename, 'physics');
assert(physics(2) <= 1, 'MC_ring_snr:Reflectivity', 'Reflectivity must be <= 1.');

SNR_THRESHOLD = 2;
TARGET_W = 0.297;                 % m, kept consistent with MC_single_pulse_snr
GAUSSIAN_CUTOFF_WIDTH = 3;         % ignore tails below exp(-18) before SNR calc
RING_RADIUS = fasan_D / 2;
RING_HALF_WIDTH = fasan_d / 2;

if fasan_d <= 0
    error('MC_ring_snr:InvalidGaussianWidth', ...
        'fasan_d must be positive for the annular Gaussian profile.');
end
if omiga <= 0
    error('MC_ring_snr:InvalidAngularSpeed', ...
        'omiga must be positive for sliding-window pulse accumulation.');
end

WINDOW_PULSES = max(1, ceil(f * fasan_D / omiga));
POST_DETECTION_BASE_ANGLE = fasan_D + 3 * fasan_d;
photon_model = get_ring_photon_model(fasan_D, fasan_d, R, physics);
carry_sig = zeros(0, 1);
carry_backscatter = zeros(0, 1);
carry_step = zeros(0, 1);
first_encounter = nan;
post_detection_photons = empty_detection_window_photons();
capture_active = false;
first_post_step = nan;
last_post_step = nan;

%% One complete beam cycle used by the original code
Beam_all = [beam_vec; path2; flipud(beam_vec); path1];
L_cycle = size(Beam_all, 1);

%% Include the initial pulse, which corresponds to first_detect = 0
idx0 = mod(init_beam_idx - 1, L_cycle) + 1;
B0 = Beam_all(idx0, :);
UP_norm = norm(UP);
UP_N = UP ./ UP_norm;

dtheta0 = safe_acos(dot(UP_N, unitvec(B0)));
illum0 = 0;
initial_candidate = is_ring_gaussian_candidate(dtheta0, UP_norm, TARGET_W, ...
    RING_RADIUS, RING_HALF_WIDTH, GAUSSIAN_CUTOFF_WIDTH);
if initial_candidate
    first_encounter = 0;
    illum0 = ring_gaussian_target_factor(dtheta0, UP_norm, ...
        TARGET_W, RING_RADIUS, RING_HALF_WIDTH);
end
[~, backscatter0] = lookup_ring_noise(photon_model, UP_norm);
sig0 = calc_ring_signal_photons(photon_model, UP_norm, illum0);
[detect, first_detect, pulse_score, carry_sig, carry_backscatter, carry_step] = ...
    append_and_check_photon_window( ...
    carry_sig, carry_backscatter, carry_step, sig0, backscatter0, 0, ...
    photon_model.N_bg, photon_model.N_dark, WINDOW_PULSES, SNR_THRESHOLD, ...
    window_mode);
if detect
    if ~collect_post_detection
        return;
    end
    post_pulse_count = post_detection_pulse_count( ...
        f, POST_DETECTION_BASE_ANGLE, omiga, UP_norm, norm(UV));
    first_post_step = first_detect + 1;
    last_post_step = first_detect + post_pulse_count;
    capture_active = true;
end

% Optional finite-horizon mode: step 0 above is the initial pulse.
if max_step_index == 0
    return;
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
    out_of_region = (norm_P_now > R) | (P_now(:,3) < 0);
    exit_idx = find(out_of_region, 1, 'first');
    if isempty(exit_idx)
        block_len = N_KERNAL;
    else
        block_len = exit_idx - 1;
    end
    step_offset = (block_id - 1) * N_KERNAL;
    remaining_steps = max_step_index - step_offset;
    block_len = min(block_len, remaining_steps);
    reached_step_limit = isfinite(max_step_index) && ...
        (step_offset + block_len >= max_step_index);

    if block_len > 0
        P_block = P_now(1:block_len, :);
        B_block = B_now(1:block_len, :);
        z_block = norm_P_now(1:block_len);

        dtheta = safe_acos_vec(dot(P_block, B_block, 2) ./ z_block);
        judge = is_ring_gaussian_candidate(dtheta, z_block, TARGET_W, ...
            RING_RADIUS, RING_HALF_WIDTH, GAUSSIAN_CUTOFF_WIDTH);

        illum_arr = zeros(block_len, 1);
        hit_idx = find(judge);
        step_arr = (block_id - 1) * N_KERNAL + (1:block_len)';
        if isnan(first_encounter) && ~isempty(hit_idx)
            first_encounter = step_arr(hit_idx(1));
        end
        if ~isempty(hit_idx)
            illum_arr(hit_idx) = ring_gaussian_target_factor( ...
                dtheta(hit_idx), z_block(hit_idx), ...
                TARGET_W, RING_RADIUS, RING_HALF_WIDTH);
        end

        % Keep the executed weighted per-pulse score criterion unchanged.
        [~, backscatter_arr] = lookup_ring_noise(photon_model, z_block);
        sig_arr = zeros(block_len, 1);
        if ~isempty(hit_idx)
            sig_arr(hit_idx) = calc_ring_signal_photons( ...
                photon_model, z_block(hit_idx), illum_arr(hit_idx));
        end

        if capture_active
            [post_detection_photons, capture_complete] = ...
                append_post_detection_photons( ...
                post_detection_photons, sig_arr, backscatter_arr, step_arr, ...
                first_post_step, last_post_step, ...
                photon_model.N_bg, photon_model.N_dark);
            if capture_complete
                return;
            end
        else
            [detect, first_detect, pulse_score, carry_sig, ...
                carry_backscatter, carry_step] = ...
                append_and_check_photon_window( ...
                carry_sig, carry_backscatter, carry_step, ...
                sig_arr, backscatter_arr, step_arr, ...
                photon_model.N_bg, photon_model.N_dark, ...
                WINDOW_PULSES, SNR_THRESHOLD, window_mode);
            if detect
                if ~collect_post_detection
                    return;
                end
                detection_pos = find(step_arr == first_detect, 1, 'first');
                detection_range = z_block(detection_pos);
                post_pulse_count = post_detection_pulse_count( ...
                    f, POST_DETECTION_BASE_ANGLE, omiga, ...
                    detection_range, norm(UV));
                first_post_step = first_detect + 1;
                last_post_step = first_detect + post_pulse_count;
                capture_active = true;

                [post_detection_photons, capture_complete] = ...
                    append_post_detection_photons( ...
                    post_detection_photons, sig_arr, backscatter_arr, ...
                    step_arr, first_post_step, last_post_step, ...
                    photon_model.N_bg, photon_model.N_dark);
                if capture_complete
                    return;
                end
            end
        end
    end

    % Original termination condition: target has left the detection region.
    if ~isempty(exit_idx) || reached_step_limit
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

function model = get_ring_photon_model(fasan_D, fasan_d, R, physics)
% Build and cache deterministic photon-count quantities for annular beams.
% The cache avoids recomputing range-gate backscatter integrals for every
% Monte Carlo trajectory.  In a parallel pool, each worker keeps its own
% persistent copy and rebuilds only when fasan_D/fasan_d/R changes.

persistent cached_model cached_key

key = sprintf('%.17g_', [fasan_D, fasan_d, R, physics(:).']);
if ~isempty(cached_model) && strcmp(cached_key, key)
    model = cached_model;
    return;
end

%% Physical constants
h = 6.62607015e-34;               % J*s
c = 299792458;                    % m/s

%% System parameters, identical to MC_single_pulse_snr.m
lambda = 1550e-9;                 % m
E_pulse = physics(1);             % J
D_rx = 0.0508;                    % m
A_rx = pi * (D_rx/2)^2;           % m^2
eta_opt = 0.80;
eta_det = 0.80;
eta_sys = eta_opt * eta_det;
tau_gate = 50e-9;                 % s
alpha = physics(3);               % 1/m
beta = physics(4);                % 1/(m*sr)
L_sky_nm = 5e-9;                  % W/(m^2*sr*nm)
delta_lambda_nm = 10.0;           % nm
dark_count_rate = 400;            % counts/s

%% Target parameters, identical to MC_single_pulse_snr.m
target_w = 0.297;                 % m
A_target = target_w * target_w;   % m^2
rho = physics(2);

E_photon = h * c / lambda;
N_tx = E_pulse / E_photon;
N_dark = dark_count_rate * tau_gate;

ring_radius = fasan_D / 2;
gaussian_width = fasan_d / 2;
Omega_beam = annular_gaussian_effective_solid_angle( ...
    ring_radius, gaussian_width);

fov_half_width = 3 * fasan_d;
theta_fov_inner = max(ring_radius - fov_half_width, 0);
theta_fov_outer = min(ring_radius + fov_half_width, pi);
Omega_fov = 2 * pi * (cos(theta_fov_inner) - cos(theta_fov_outer));

N_bg = calc_background_photons( ...
    L_sky_nm, delta_lambda_nm, A_rx, Omega_fov, tau_gate, eta_sys, E_photon);

z_step = 0.5;
z_lut = (0:z_step:R).';
if isempty(z_lut) || z_lut(end) < R
    z_lut = [z_lut; R];
end
N_bs_lut = calc_backscatter_photons_vec( ...
    N_tx, eta_sys, A_rx, beta, alpha, z_lut, tau_gate, c);

model = struct( ...
    'z_step', z_step, ...
    'z_lut', z_lut, ...
    'N_bs_lut', N_bs_lut, ...
    'N_bg', N_bg, ...
    'N_dark', N_dark, ...
    'N_tx', N_tx, ...
    'A_target', A_target, ...
    'rho', rho, ...
    'A_rx', A_rx, ...
    'eta_sys', eta_sys, ...
    'alpha', alpha, ...
    'Omega_beam', Omega_beam);

cached_key = key;
cached_model = model;
end

function [N_noise, N_bs] = lookup_ring_noise(model, z)
% Nearest-neighbor lookup on the 0.5 m deterministic backscatter table.

z = z(:);
idx = round(z ./ model.z_step) + 1;
idx(~isfinite(idx)) = 1;
idx = max(1, min(numel(model.N_bs_lut), idx));
N_bs = model.N_bs_lut(idx);
N_noise = N_bs + model.N_bg + model.N_dark;
end

function N_sig = calc_ring_signal_photons(model, z, illumination_factor)
% Deterministic target-return photons for annular Gaussian illumination.

z = z(:);
illumination_factor = illumination_factor(:);
N_sig = zeros(size(z));

valid = isfinite(z) & z > 0 & isfinite(illumination_factor) & ...
    illumination_factor > 0 & model.Omega_beam > 0;
if ~any(valid)
    return;
end

beam_area = model.Omega_beam .* z(valid).^2;
T2_target = exp(-2 * model.alpha .* z(valid));
N_sig(valid) = model.N_tx .* ...
    (model.A_target .* illumination_factor(valid) ./ beam_area) .* ...
    model.rho .* (model.A_rx ./ (pi .* z(valid).^2)) .* ...
    T2_target .* model.eta_sys;
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

function y = safe_acos(x)
y = acos(max(min(x, 1), -1));
end

function y = safe_acos_vec(x)
y = acos(max(min(x, 1), -1));
end

function Omega_eff = annular_gaussian_effective_solid_angle(ring_radius, gaussian_width)
% Integral of exp(-2*((theta-ring_radius)/gaussian_width)^2) over angle.

if gaussian_width <= 0
    Omega_eff = nan;
    return;
end

u0 = -sqrt(2) * ring_radius / gaussian_width;
radial_integral = ring_radius * gaussian_width * sqrt(pi) / (2 * sqrt(2)) .* ...
    erfc(u0) + gaussian_width.^2 / 4 .* exp(-2 * (ring_radius ./ gaussian_width).^2);
Omega_eff = 2 * pi .* radial_integral;
end

function N_bs = calc_backscatter_photons_vec(N_tx, eta_sys, A_rx, beta, alpha, z, tau_gate, c)
% Vectorized deterministic backscatter expectation over the range gate.

dr = c * tau_gate / 2;
offsets = linspace(-dr/2, dr/2, 101);
r = z(:) + offsets;
r(r <= 0) = 1e-6;
integrand = beta .* exp(-2 .* alpha .* r) ./ (r.^2);
integ = trapz(offsets, integrand, 2);
N_bs = N_tx .* eta_sys .* A_rx .* integ;
end

function N_bg = calc_background_photons(L_sky_nm, delta_lambda_nm, A_rx, Omega_fov, tau_gate, eta_sys, E_photon)
P_bg = L_sky_nm .* delta_lambda_nm .* A_rx .* Omega_fov;
E_bg = P_bg .* tau_gate;
N_bg = E_bg .* eta_sys ./ E_photon;
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

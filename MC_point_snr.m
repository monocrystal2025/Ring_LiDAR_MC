function [detect, first_detect, first_encounter, pulse_score, ...
    post_detection_photons] = MC_point_snr( ...
    fasan_D, ~, path1, beam_vec, path2, f, omiga, ...
    UP, UV, init_beam_idx, R, N_KERNAL, max_step_index, window_mode)
%MC_POINT_SNR  Monte Carlo detection judgment for a Gaussian point beam.
%
% fasan_D is the full divergence angle at the exp(-2) intensity contour.
% The Gaussian half-angle width is therefore fasan_D/2.  The receiver FOV is
% a hard circular cone with full angle 3*fasan_D.
% first_encounter is the pulse-index difference from region entry (pulse 0)
% to the first geometric-prefilter hit. pulse_score follows the selected
% window mode; post_detection_photons stores following pulses when requested.

if nargin < 12 || isempty(N_KERNAL)
    N_KERNAL = 100000;
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

SNR_THRESHOLD = 2;
TARGET_W = 0.297;             % m, kept consistent with ring/line SNR code
GAUSSIAN_HALF_ANGLE = fasan_D / 2;
FOV_HALF_ANGLE = 3 * GAUSSIAN_HALF_ANGLE;

if fasan_D <= 0
    error('MC_point_snr:InvalidGaussianWidth', ...
        'fasan_D must be positive for the point-beam Gaussian profile.');
end
if omiga <= 0
    error('MC_point_snr:InvalidAngularSpeed', ...
        'omiga must be positive for sliding-window pulse accumulation.');
end

WINDOW_PULSES = max(1, ceil(f * (3 * fasan_D) / omiga));
fasan_d = 0;
POST_DETECTION_BASE_ANGLE = fasan_D + 3 * fasan_d;
photon_model = get_point_photon_model(fasan_D, R);
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

cos_theta0 = dot(UP_N, unitvec(B0));
illum0 = 0;
initial_candidate = is_point_gaussian_candidate_from_cos(cos_theta0, ...
    UP_norm, TARGET_W, FOV_HALF_ANGLE);
if initial_candidate
    first_encounter = 0;
    dtheta0 = safe_acos(cos_theta0);
    illum0 = point_gaussian_target_factor(dtheta0, UP_norm, ...
        TARGET_W, GAUSSIAN_HALF_ANGLE, FOV_HALF_ANGLE);
end
[~, backscatter0] = lookup_point_noise(photon_model, UP_norm);
sig0 = calc_point_signal_photons(photon_model, UP_norm, illum0);
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

        cos_theta = dot(P_block, B_block, 2) ./ z_block;
        judge = is_point_gaussian_candidate_from_cos(cos_theta, z_block, ...
            TARGET_W, FOV_HALF_ANGLE);

        illum_arr = zeros(block_len, 1);
        hit_idx = find(judge);
        step_arr = (block_id - 1) * N_KERNAL + (1:block_len)';
        if isnan(first_encounter) && ~isempty(hit_idx)
            first_encounter = step_arr(hit_idx(1));
        end
        if ~isempty(hit_idx)
            dtheta = safe_acos_vec(cos_theta(hit_idx));
            illum_arr(hit_idx) = point_gaussian_target_factor( ...
                dtheta, z_block(hit_idx), TARGET_W, ...
                GAUSSIAN_HALF_ANGLE, FOV_HALF_ANGLE);
        end

        [~, backscatter_arr] = lookup_point_noise(photon_model, z_block);
        sig_arr = zeros(block_len, 1);
        if ~isempty(hit_idx)
            sig_arr(hit_idx) = calc_point_signal_photons( ...
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

function model = get_point_photon_model(fasan_D, R)
% Build and cache deterministic photon-count quantities for point beams.

persistent cached_model cached_key

key = sprintf('%.17g_%.17g', fasan_D, R);
if ~isempty(cached_model) && strcmp(cached_key, key)
    model = cached_model;
    return;
end

%% Physical constants
h = 6.62607015e-34;               % J*s
c = 299792458;                    % m/s

%% System parameters, identical to ring/line SNR code
lambda = 1550e-9;                 % m
E_pulse = 150e-6;                 % J
D_rx = 0.0508;                    % m
A_rx = pi * (D_rx/2)^2;           % m^2
eta_opt = 0.80;
eta_det = 0.80;
eta_sys = eta_opt * eta_det;
tau_gate = 50e-9;                 % s
alpha = 1.5e-5;                   % 1/m
beta = 0.3e-6;                    % 1/(m*sr)
L_sky_nm = 5e-9;                  % W/(m^2*sr*nm)
delta_lambda_nm = 10.0;           % nm
dark_count_rate = 400;            % counts/s

%% Target parameters, identical to ring/line SNR code
target_w = 0.297;                 % m
A_target = target_w * target_w;   % m^2
rho = 0.5;

E_photon = h * c / lambda;
N_tx = E_pulse / E_photon;
N_dark = dark_count_rate * tau_gate;

gaussian_half_angle = fasan_D / 2;
fov_half_angle = min(3 * gaussian_half_angle, pi);
Omega_beam = point_gaussian_effective_solid_angle(gaussian_half_angle);
Omega_fov = 2 * pi * (1 - cos(fov_half_angle));

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

function [N_noise, N_bs] = lookup_point_noise(model, z)
% Nearest-neighbor lookup on the 0.5 m deterministic backscatter table.

z = z(:);
idx = round(z ./ model.z_step) + 1;
idx(~isfinite(idx)) = 1;
idx = max(1, min(numel(model.N_bs_lut), idx));
N_bs = model.N_bs_lut(idx);
N_noise = N_bs + model.N_bg + model.N_dark;
end

function N_sig = calc_point_signal_photons(model, z, illumination_factor)
% Deterministic target-return photons for circular Gaussian point light.

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

function Omega_eff = point_gaussian_effective_solid_angle(gaussian_half_angle)
% Exact spherical solid angle of exp(-2*(theta/theta_w)^2).

if gaussian_half_angle <= 0
    Omega_eff = nan;
    return;
end

integrand = @(theta) exp(-2 .* (theta ./ gaussian_half_angle).^2) .* sin(theta);
Omega_eff = 2 * pi * integral(integrand, 0, pi, ...
    'RelTol', 1e-10, 'AbsTol', 1e-14);
end

function candidate = is_point_gaussian_candidate_from_cos(cos_theta, z, ...
    target_w, fov_half_angle)
% Fast FOV prefilter.  It avoids acos for pulses that are clearly outside
% the receiver cone plus the target's angular half-diagonal.

z = z(:);
cos_theta = cos_theta(:);
target_half_diag_theta = atan((sqrt(2) * target_w / 2) ./ z);
theta_limit = min(fov_half_angle + target_half_diag_theta, pi);

candidate = cos_theta >= cos(theta_limit);
candidate = candidate & isfinite(cos_theta) & isfinite(z) & z > 0;
end

function factor = point_gaussian_target_factor(center_theta, z, target_w, ...
    gaussian_half_angle, fov_half_angle)
% Average normalized circular-Gaussian intensity over a square target.

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

valid = isfinite(center_theta) & isfinite(z) & z > 0 & ...
    gaussian_half_angle > 0;
if ~any(valid)
    return;
end

angular_side = 2 .* atan((target_w / 2) ./ z(valid));
dx = angular_side .* sample_x;
dy = angular_side .* sample_y;
sample_theta = sqrt((center_theta(valid) + dx).^2 + dy.^2);
inside_fov = sample_theta <= fov_half_angle;
sample_intensity = inside_fov .* ...
    exp(-2 .* (sample_theta ./ gaussian_half_angle).^2);
factor(valid) = mean(sample_intensity, 2);
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

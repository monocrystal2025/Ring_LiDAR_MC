function [detect, first_detect, first_encounter, pulse_score] = ...
    MC_point_RA_snr_NEW(fasan_D, ~, ...
    beam_vec, f, omiga, UP, UV, init_beam_idx, R, jiaodu, window_mode, mc_block_size)
%MC_POINT_RA_SNR_NEW Gaussian point-beam detection in a cone ROI.

if nargin < 11 || isempty(window_mode)
    window_mode = 'variable';
end
if nargin < 12 || isempty(mc_block_size)
    mc_block_size = 100000;
end
validateattributes(mc_block_size, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'positive'}, ...
    mfilename, 'mc_block_size');
N_KERNAL = mc_block_size;
SNR_THRESHOLD = 2;
TARGET_W = 0.297;
GAUSSIAN_HALF_ANGLE = fasan_D / 2;
FOV_HALF_ANGLE = 3 * GAUSSIAN_HALF_ANGLE;

if fasan_D <= 0
    error('MC_point_RA_snr_NEW:InvalidGaussianWidth', ...
        'fasan_D must be positive.');
end
if omiga <= 0
    error('MC_point_RA_snr_NEW:InvalidAngularSpeed', ...
        'omiga must be positive.');
end

WINDOW_PULSES = max(1, ceil(f * (3 * fasan_D) / omiga));
photon_model = get_point_photon_model(fasan_D, R);
carry_sig = zeros(0, 1);
carry_backscatter = zeros(0, 1);
carry_step = zeros(0, 1);
first_encounter = nan;
pulse_score = empty_detection_window_photons();

Beam_all = beam_vec;
idx0 = mod(init_beam_idx - 1, size(Beam_all, 1)) + 1;
B0 = Beam_all(idx0, :);

if is_inside_cone(UP, R, jiaodu)
    range0 = norm(UP);
    U0 = UP ./ range0;
    cos_theta0 = dot(U0, B0);
    illum0 = 0;
    if is_point_gaussian_candidate_from_cos(cos_theta0, range0, ...
            TARGET_W, FOV_HALF_ANGLE)
        first_encounter = 0;
        dtheta0 = safe_acos(cos_theta0);
        illum0 = point_gaussian_target_factor(dtheta0, range0, ...
            TARGET_W, GAUSSIAN_HALF_ANGLE, FOV_HALF_ANGLE);
    end
    [~, backscatter0] = lookup_point_noise(photon_model, range0);
    sig0 = calc_point_signal_photons(photon_model, range0, illum0);
    [detect, first_detect, pulse_score, carry_sig, ...
        carry_backscatter, carry_step] = ...
        append_and_check_photon_window( ...
        carry_sig, carry_backscatter, carry_step, ...
        sig0, backscatter0, 0, photon_model.N_bg, ...
        photon_model.N_dark, WINDOW_PULSES, SNR_THRESHOLD, window_mode);
    if detect
        return;
    end
else
    detect = 0;
    first_detect = nan;
    return;
end

P_now = repmat(UP, N_KERNAL, 1);
dt = (1:N_KERNAL)' ./ f;
P_now = P_now + [UV(1) * dt, UV(2) * dt, UV(3) * dt];

B_now = tangent_path_segment(beam_vec, init_beam_idx + 1, N_KERNAL);
beam_idx = init_beam_idx + 1;
block_id = 1;

while true
    inside = is_inside_cone(P_now, R, jiaodu);
    exit_idx = find(~inside, 1, 'first');
    if isempty(exit_idx)
        block_len = N_KERNAL;
    else
        block_len = exit_idx - 1;
    end

    if block_len > 0
        P_block = P_now(1:block_len, :);
        B_block = B_now(1:block_len, :);
        z_block = vecnorm(P_block, 2, 2);
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

        [detect, first_detect, pulse_score, carry_sig, ...
            carry_backscatter, carry_step] = ...
            append_and_check_photon_window( ...
            carry_sig, carry_backscatter, carry_step, ...
            sig_arr, backscatter_arr, step_arr, ...
            photon_model.N_bg, photon_model.N_dark, ...
            WINDOW_PULSES, SNR_THRESHOLD, window_mode);
        if detect
            return;
        end
    end

    if ~isempty(exit_idx)
        break;
    end

    block_id = block_id + 1;
    P_now = P_now + N_KERNAL .* UV ./ f;
    beam_idx = beam_idx + N_KERNAL;
    B_now = tangent_path_segment(beam_vec, beam_idx, N_KERNAL);
end

detect = 0;
first_detect = nan;
end

function inside = is_inside_cone(P, R, jiaodu)
alpha_tan = tand(jiaodu);
if isvector(P)
    P = P(:).';
end
range = vecnorm(P, 2, 2);
rho = hypot(P(:, 1), P(:, 2));
inside = P(:, 3) >= 0 & range <= R & rho <= P(:, 3) .* alpha_tan + 1e-10;
end

function seg = tangent_path_segment(path, k, m)
L = size(path, 1);
idx = mod((k - 1) + (0:m-1), L) + 1;
seg = path(idx, :);
end

function model = get_point_photon_model(fasan_D, R)
persistent cached_model cached_key
key = sprintf('%.17g_%.17g', fasan_D, R);
if ~isempty(cached_model) && strcmp(cached_key, key)
    model = cached_model;
    return;
end

h = 6.62607015e-34;
c = 299792458;
lambda = 1550e-9;
E_pulse = 150e-6;
D_rx = 0.0508;
A_rx = pi * (D_rx / 2)^2;
eta_sys = 0.80 * 0.80;
tau_gate = 50e-9;
alpha = 1.5e-5;
beta = 0.3e-6;
L_sky_nm = 5e-9;
delta_lambda_nm = 10.0;
dark_count_rate = 400;
target_w = 0.297;
A_target = target_w * target_w;
rho = 0.5;

E_photon = h * c / lambda;
N_tx = E_pulse / E_photon;
N_dark = dark_count_rate * tau_gate;
gaussian_half_angle = fasan_D / 2;
fov_half_angle = min(3 * gaussian_half_angle, pi);
Omega_beam = point_gaussian_effective_solid_angle(gaussian_half_angle);
Omega_fov = 2 * pi * (1 - cos(fov_half_angle));
N_bg = calc_background_photons(L_sky_nm, delta_lambda_nm, A_rx, ...
    Omega_fov, tau_gate, eta_sys, E_photon);

z_step = 0.5;
z_lut = (0:z_step:R).';
if isempty(z_lut) || z_lut(end) < R
    z_lut = [z_lut; R];
end
N_bs_lut = calc_backscatter_photons_vec(N_tx, eta_sys, A_rx, beta, ...
    alpha, z_lut, tau_gate, c);

model = struct('z_step', z_step, 'z_lut', z_lut, ...
    'N_bs_lut', N_bs_lut, 'N_bg', N_bg, 'N_dark', N_dark, ...
    'N_tx', N_tx, 'A_target', A_target, 'rho', rho, ...
    'A_rx', A_rx, 'eta_sys', eta_sys, 'alpha', alpha, ...
    'Omega_beam', Omega_beam);

cached_key = key;
cached_model = model;
end

function [N_noise, N_bs] = lookup_point_noise(model, z)
z = z(:);
idx = round(z ./ model.z_step) + 1;
idx(~isfinite(idx)) = 1;
idx = max(1, min(numel(model.N_bs_lut), idx));
N_bs = model.N_bs_lut(idx);
N_noise = N_bs + model.N_bg + model.N_dark;
end

function N_sig = calc_point_signal_photons(model, z, illumination_factor)
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

function candidate = is_point_gaussian_candidate_from_cos(cos_theta, z, ...
    target_w, fov_half_angle)
z = z(:);
cos_theta = cos_theta(:);
target_half_diag_theta = atan((sqrt(2) * target_w / 2) ./ z);
theta_limit = min(fov_half_angle + target_half_diag_theta, pi);
candidate = cos_theta >= cos(theta_limit);
candidate = candidate & isfinite(cos_theta) & isfinite(z) & z > 0;
end

function factor = point_gaussian_target_factor(center_theta, z, target_w, ...
    gaussian_half_angle, fov_half_angle)
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

function Omega_eff = point_gaussian_effective_solid_angle(gaussian_half_angle)
if gaussian_half_angle <= 0
    Omega_eff = nan;
    return;
end
integrand = @(theta) exp(-2 .* (theta ./ gaussian_half_angle).^2) .* sin(theta);
Omega_eff = 2 * pi * integral(integrand, 0, pi, ...
    'RelTol', 1e-10, 'AbsTol', 1e-14);
end

function y = safe_acos(x)
y = acos(max(min(x, 1), -1));
end

function y = safe_acos_vec(x)
y = acos(max(min(x, 1), -1));
end

function N_bs = calc_backscatter_photons_vec(N_tx, eta_sys, A_rx, beta, alpha, z, tau_gate, c)
dr = c * tau_gate / 2;
offsets = linspace(-dr / 2, dr / 2, 101);
r = z(:) + offsets;
r(r <= 0) = 1e-6;
integrand = beta .* exp(-2 .* alpha .* r) ./ (r.^2);
integ = trapz(offsets, integrand, 2);
N_bs = N_tx .* eta_sys .* A_rx .* integ;
end

function N_bg = calc_background_photons(L_sky_nm, delta_lambda_nm, A_rx, ...
    Omega_fov, tau_gate, eta_sys, E_photon)
P_bg = L_sky_nm .* delta_lambda_nm .* A_rx .* Omega_fov;
E_bg = P_bg .* tau_gate;
N_bg = E_bg .* eta_sys ./ E_photon;
end

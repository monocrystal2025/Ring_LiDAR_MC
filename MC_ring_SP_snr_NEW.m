function [detect, first_detect, first_encounter, effective_pulses] = ...
    MC_ring_SP_snr_NEW(fasan_D, fasan_d, ...
    beam_vec, f, omiga, UP, UV, init_beam_idx, R, jiaodu)
%MC_RING_SP_SNR_NEW Annular-beam spiral detection in a cone ROI.

[detect, first_detect, first_encounter, effective_pulses] = ...
    cone_ring_snr_core(fasan_D, fasan_d, ...
    beam_vec, f, omiga, UP, UV, init_beam_idx, R, jiaodu);
end

function [detect, first_detect, first_encounter, effective_pulses] = ...
    cone_ring_snr_core(fasan_D, fasan_d, ...
    beam_vec, f, omiga, UP, UV, init_beam_idx, R, jiaodu)

N_KERNAL = 100000;
SNR_THRESHOLD = 2;
TARGET_W = 0.297;
GAUSSIAN_CUTOFF_WIDTH = 3;
RING_RADIUS = fasan_D / 2;

if fasan_d <= 0
    error('MC_ring_SP_snr_NEW:InvalidGaussianWidth', ...
        'fasan_d must be positive.');
end
if omiga <= 0
    error('MC_ring_SP_snr_NEW:InvalidAngularSpeed', ...
        'omiga must be positive.');
end

gaussian_width = fasan_d / 2;
fov_ring_half_width = 3 * fasan_d;
WINDOW_PULSES = max(1, ceil(f * fasan_D / omiga));
photon_model = get_ring_photon_model(fasan_D, fasan_d, R);
carry_sig = zeros(0, 1);
carry_noise = zeros(0, 1);
carry_step = zeros(0, 1);
first_encounter = nan;
effective_pulses = nan;

Beam_all = beam_vec;
idx0 = mod(init_beam_idx - 1, size(Beam_all, 1)) + 1;
B0 = Beam_all(idx0, :);

if is_inside_cone(UP, R, jiaodu)
    range0 = norm(UP);
    dtheta0 = safe_acos(dot(UP ./ range0, B0));
    illum0 = 0;
    if is_ring_gaussian_candidate(dtheta0, range0, TARGET_W, ...
            RING_RADIUS, gaussian_width, GAUSSIAN_CUTOFF_WIDTH)
        first_encounter = 0;
        illum0 = ring_gaussian_target_factor(dtheta0, range0, ...
            TARGET_W, RING_RADIUS, gaussian_width, fov_ring_half_width);
    end
    noise0 = lookup_ring_noise(photon_model, range0);
    sig0 = calc_ring_signal_photons(photon_model, range0, illum0);
    [detect, first_detect, effective_pulses, ...
        carry_sig, carry_noise, carry_step] = ...
        append_and_check_window(carry_sig, carry_noise, carry_step, ...
        sig0, noise0, 0, WINDOW_PULSES, SNR_THRESHOLD);
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
        P_unit = P_block ./ z_block;

        dtheta = safe_acos_vec(dot(P_unit, B_block, 2));
        judge = is_ring_gaussian_candidate(dtheta, z_block, TARGET_W, ...
            RING_RADIUS, gaussian_width, GAUSSIAN_CUTOFF_WIDTH);

        illum_arr = zeros(block_len, 1);
        hit_idx = find(judge);
        step_arr = (block_id - 1) * N_KERNAL + (1:block_len)';
        if isnan(first_encounter) && ~isempty(hit_idx)
            first_encounter = step_arr(hit_idx(1));
        end
        if ~isempty(hit_idx)
            illum_arr(hit_idx) = ring_gaussian_target_factor( ...
                dtheta(hit_idx), z_block(hit_idx), ...
                TARGET_W, RING_RADIUS, gaussian_width, ...
                fov_ring_half_width);
        end

        noise_arr = lookup_ring_noise(photon_model, z_block);
        sig_arr = zeros(block_len, 1);
        if ~isempty(hit_idx)
            sig_arr(hit_idx) = calc_ring_signal_photons( ...
                photon_model, z_block(hit_idx), illum_arr(hit_idx));
        end

        [detect, first_detect, effective_pulses, ...
            carry_sig, carry_noise, carry_step] = ...
            append_and_check_window(carry_sig, carry_noise, carry_step, ...
            sig_arr, noise_arr, step_arr, WINDOW_PULSES, SNR_THRESHOLD);
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

function model = get_ring_photon_model(fasan_D, fasan_d, R)
persistent cached_model cached_key

key = sprintf('%.17g_%.17g_%.17g', fasan_D, fasan_d, R);
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
eta_opt = 0.80;
eta_det = 0.80;
eta_sys = eta_opt * eta_det;
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

ring_radius = fasan_D / 2;
gaussian_width = fasan_d / 2;
Omega_beam = annular_gaussian_effective_solid_angle(ring_radius, gaussian_width);
fov_half_width = 3 * fasan_d;
theta_fov_inner = max(ring_radius - fov_half_width, 0);
theta_fov_outer = min(ring_radius + fov_half_width, pi);
Omega_fov = 2 * pi * (cos(theta_fov_inner) - cos(theta_fov_outer));
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

function N_noise = lookup_ring_noise(model, z)
z = z(:);
idx = round(z ./ model.z_step) + 1;
idx(~isfinite(idx)) = 1;
idx = max(1, min(numel(model.N_bs_lut), idx));
N_noise = model.N_bs_lut(idx) + model.N_bg + model.N_dark;
end

function N_sig = calc_ring_signal_photons(model, z, illumination_factor)
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

function [detect, first_detect, effective_pulses, ...
    carry_sig, carry_noise, carry_step] = ...
    append_and_check_window(carry_sig, carry_noise, carry_step, ...
    sig_arr, noise_arr, step_arr, window_pulses, threshold)

detect = 0;
first_detect = nan;
effective_pulses = nan;
sig_arr = sig_arr(:);
noise_arr = noise_arr(:);
step_arr = step_arr(:);
if isempty(sig_arr)
    return;
end
sig_arr(~isfinite(sig_arr)) = 0;
noise_arr(~isfinite(noise_arr)) = 0;

n_carry = numel(carry_sig);
all_sig = [carry_sig; sig_arr];
all_noise = [carry_noise; noise_arr];
all_step = [carry_step; step_arr];
n_all = numel(all_sig);
signal_tail_pos = n_carry + find(sig_arr > 0);

if ~isempty(signal_tail_pos)
    tail_pos = signal_tail_pos(:);
    max_len = min(window_pulses, n_all);
    len_vec = 1:max_len;
    start_pos = tail_pos - len_vec + 1;
    valid_len = start_pos >= 1;
    start_pos(~valid_len) = 1;

    pulse_var = all_sig + all_noise;
    pulse_score = zeros(size(all_sig));
    score_valid = all_sig > 0 & pulse_var > 0 & isfinite(pulse_var);
    pulse_score(score_valid) = all_sig(score_valid).^2 ./ pulse_var(score_valid);
    cs_score = [0; cumsum(pulse_score)];
    win_score = cs_score(tail_pos + 1) - cs_score(start_pos);
    win_snr = sqrt(win_score);
    win_snr(~valid_len | win_score <= 0 | ~isfinite(win_snr)) = -inf;

    hit_matrix = win_snr >= threshold;
    ok = find(any(hit_matrix, 2), 1, 'first');
    if ~isempty(ok)
        successful_len = find(hit_matrix(ok, :), 1, 'first');
        successful_start = start_pos(ok, successful_len);
        detect = 1;
        first_detect = all_step(tail_pos(ok));
        effective_pulses = nnz( ...
            pulse_score(successful_start:tail_pos(ok)) > 0);
        return;
    end
end

n_keep = min(window_pulses - 1, n_all);
if n_keep > 0
    carry_sig = all_sig(end - n_keep + 1:end);
    carry_noise = all_noise(end - n_keep + 1:end);
    carry_step = all_step(end - n_keep + 1:end);
else
    carry_sig = zeros(0, 1);
    carry_noise = zeros(0, 1);
    carry_step = zeros(0, 1);
end
end

function y = safe_acos(x)
y = acos(max(min(x, 1), -1));
end

function y = safe_acos_vec(x)
y = acos(max(min(x, 1), -1));
end

function Omega_eff = annular_gaussian_effective_solid_angle(ring_radius, gaussian_width)
if gaussian_width <= 0
    Omega_eff = nan;
    return;
end
u0 = -sqrt(2) * ring_radius / gaussian_width;
radial_integral = ring_radius * gaussian_width * sqrt(pi) / (2 * sqrt(2)) .* ...
    erfc(u0) + gaussian_width.^2 / 4 .* ...
    exp(-2 * (ring_radius ./ gaussian_width).^2);
Omega_eff = 2 * pi .* radial_integral;
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

function candidate = is_ring_gaussian_candidate(center_theta, z, ...
    target_w, ring_radius, gaussian_width, cutoff_width)
z = z(:);
center_theta = center_theta(:);
target_half_diag_theta = atan((sqrt(2) * target_w / 2) ./ z);
candidate = abs(center_theta - ring_radius) <= ...
    (target_half_diag_theta + cutoff_width * gaussian_width);
candidate = candidate & isfinite(center_theta) & isfinite(z) & z > 0;
end

function factor = ring_gaussian_target_factor(center_theta, z, ...
    target_w, ring_radius, gaussian_width, fov_half_width)
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
inside_fov = abs(sample_theta - ring_radius) <= fov_half_width;
sample_intensity = inside_fov .* ...
    exp(-2 .* ((sample_theta - ring_radius) ./ gaussian_width).^2);
factor(valid) = mean(sample_intensity, 2);
end

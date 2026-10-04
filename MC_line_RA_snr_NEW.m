function [detect, first_detect, first_encounter, pulse_score] = ...
    MC_line_RA_snr_NEW(fasan_D, fasan_d, ...
    line_scan, f, omiga, UP, UV, init_beam_idx, R, jiaodu, window_mode, mc_block_size)
%MC_LINE_RA_SNR_NEW Line-beam detection in a cone ROI.

if nargin < 11 || isempty(window_mode)
    window_mode = 'variable';
end
if nargin < 12 || isempty(mc_block_size)
    mc_block_size = 100000;
end
validateattributes(mc_block_size, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'positive'}, ...
    mfilename, 'mc_block_size');
[detect, first_detect, first_encounter, pulse_score] = ...
    cone_line_snr_core(fasan_D, fasan_d, ...
    line_scan, f, omiga, UP, UV, init_beam_idx, R, jiaodu, window_mode, mc_block_size);
end

function [detect, first_detect, first_encounter, pulse_score] = ...
    cone_line_snr_core(fasan_D, fasan_d, ...
    line_scan, f, omiga, UP, UV, init_beam_idx, R, jiaodu, window_mode, mc_block_size)

N_KERNAL = mc_block_size;
SNR_THRESHOLD = 2;
TARGET_W = 0.297;
FOV_HALF_WIDTH_MULTIPLIER = 3;
GAUSSIAN_CUTOFF_WIDTH = 3;

if fasan_D <= 0
    error('MC_line_RA_snr_NEW:InvalidLongWidth', ...
        'fasan_D must be positive.');
end
if fasan_d <= 0
    error('MC_line_RA_snr_NEW:InvalidGaussianWidth', ...
        'fasan_d must be positive.');
end
if omiga <= 0
    error('MC_line_RA_snr_NEW:InvalidAngularSpeed', ...
        'omiga must be positive.');
end

long_half_width = fasan_D / 2;
gaussian_width = fasan_d / 2;
fov_short_half_width = FOV_HALF_WIDTH_MULTIPLIER * fasan_d;
short_prefilter_half_width = GAUSSIAN_CUTOFF_WIDTH * gaussian_width;
WINDOW_PULSES = max(1, ceil(f * fasan_D / omiga));

photon_model = get_line_photon_model(fasan_D, fasan_d, gaussian_width, R);
carry_sig = zeros(0, 1);
carry_backscatter = zeros(0, 1);
carry_step = zeros(0, 1);
first_encounter = nan;
pulse_score = empty_detection_window_photons();

P_cycle = line_scan.Directions;
eL_cycle = line_scan.LongAxes;
eS_cycle = line_scan.ShortAxes;
idx0 = mod(init_beam_idx - 1, size(P_cycle, 1)) + 1;
B0 = P_cycle(idx0, :);
EL0 = eL_cycle(idx0, :);
ES0 = eS_cycle(idx0, :);

if is_inside_cone(UP, R, jiaodu)
    range0 = norm(UP);
    U0 = UP ./ range0;
    [long0, short0, front0] = line_offsets_from_frame(U0, B0, EL0, ES0);
    illum0 = 0;
    if front0 && is_line_gaussian_candidate(long0, short0, range0, ...
            TARGET_W, long_half_width, short_prefilter_half_width)
        first_encounter = 0;
        illum0 = line_gaussian_target_factor(long0, short0, range0, ...
            TARGET_W, long_half_width, gaussian_width, ...
            fov_short_half_width);
    end
    [~, backscatter0] = lookup_line_noise(photon_model, range0);
    sig0 = calc_line_signal_photons(photon_model, range0, illum0);
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

[B_now, EL_now, ES_now] = tangent_frame_segment(P_cycle, eL_cycle, ...
    eS_cycle, init_beam_idx + 1, N_KERNAL);
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
        EL_block = EL_now(1:block_len, :);
        ES_block = ES_now(1:block_len, :);
        z_block = vecnorm(P_block, 2, 2);
        U_block = P_block ./ z_block;

        [long_offset, short_offset, front] = line_offsets_from_frame( ...
            U_block, B_block, EL_block, ES_block);
        judge = front & is_line_gaussian_candidate(long_offset, short_offset, ...
            z_block, TARGET_W, long_half_width, short_prefilter_half_width);

        illum_arr = zeros(block_len, 1);
        hit_idx = find(judge);
        step_arr = (block_id - 1) * N_KERNAL + (1:block_len)';
        if isnan(first_encounter) && ~isempty(hit_idx)
            first_encounter = step_arr(hit_idx(1));
        end
        if ~isempty(hit_idx)
            illum_arr(hit_idx) = line_gaussian_target_factor( ...
                long_offset(hit_idx), short_offset(hit_idx), ...
                z_block(hit_idx), TARGET_W, long_half_width, ...
                gaussian_width, fov_short_half_width);
        end

        [~, backscatter_arr] = lookup_line_noise(photon_model, z_block);
        sig_arr = zeros(block_len, 1);
        if ~isempty(hit_idx)
            sig_arr(hit_idx) = calc_line_signal_photons( ...
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
    [B_now, EL_now, ES_now] = tangent_frame_segment(P_cycle, eL_cycle, ...
        eS_cycle, beam_idx, N_KERNAL);
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

function [B, E_long, E_short] = tangent_frame_segment(P_cycle, eL_cycle, eS_cycle, k, m)
L = size(P_cycle, 1);
idx = mod((k - 1) + (0:m-1), L) + 1;
B = P_cycle(idx, :);
E_long = eL_cycle(idx, :);
E_short = eS_cycle(idx, :);
end

function [long_offset, short_offset, front] = line_offsets_from_frame(U_dir, B_dir, E_long, E_short)
U_dir = U_dir ./ vecnorm(U_dir, 2, 2);
B_dir = B_dir ./ vecnorm(B_dir, 2, 2);
c3 = sum(U_dir .* B_dir, 2);
c_long = sum(U_dir .* E_long, 2);
c_short = sum(U_dir .* E_short, 2);
front = c3 > 0;
long_offset = atan2(c_long, c3);
short_offset = atan2(c_short, c3);
end

function model = get_line_photon_model(fasan_D, fasan_d, gaussian_width, R)
persistent cached_model cached_key
key = sprintf('%.17g_%.17g_%.17g_%.17g', fasan_D, fasan_d, gaussian_width, R);
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
Omega_beam = fasan_D * gaussian_width * sqrt(pi / 2);
fov_short_full_width = 6 * fasan_d;
Omega_fov = fasan_D * fov_short_full_width;
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

function [N_noise, N_bs] = lookup_line_noise(model, z)
z = z(:);
idx = round(z ./ model.z_step) + 1;
idx(~isfinite(idx)) = 1;
idx = max(1, min(numel(model.N_bs_lut), idx));
N_bs = model.N_bs_lut(idx);
N_noise = N_bs + model.N_bg + model.N_dark;
end

function N_sig = calc_line_signal_photons(model, z, illumination_factor)
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

function candidate = is_line_gaussian_candidate(long_offset, short_offset, z, ...
    target_w, long_half_width, short_prefilter_half_width)
z = z(:);
long_offset = long_offset(:);
short_offset = short_offset(:);
target_half_theta = atan((target_w / 2) ./ z);
candidate = abs(long_offset) <= (long_half_width + target_half_theta) & ...
    abs(short_offset) <= (short_prefilter_half_width + target_half_theta);
candidate = candidate & isfinite(long_offset) & isfinite(short_offset) & ...
    isfinite(z) & z > 0;
end

function factor = line_gaussian_target_factor(long_offset, short_offset, z, ...
    target_w, long_half_width, gaussian_width, fov_short_half_width)
persistent sample_long sample_short
if isempty(sample_long)
    grid_n = 7;
    nodes = ((1:grid_n) - 0.5) ./ grid_n - 0.5;
    [ll, ss] = meshgrid(nodes, nodes);
    sample_long = ll(:).';
    sample_short = ss(:).';
end

long_offset = long_offset(:);
short_offset = short_offset(:);
z = z(:);
factor = zeros(size(long_offset));
valid = isfinite(long_offset) & isfinite(short_offset) & ...
    isfinite(z) & z > 0 & gaussian_width > 0;
if ~any(valid)
    return;
end
angular_side = 2 .* atan((target_w / 2) ./ z(valid));
sample_l = long_offset(valid) + angular_side .* sample_long;
sample_s = short_offset(valid) + angular_side .* sample_short;
inside_long = abs(sample_l) <= long_half_width;
inside_fov = abs(sample_s) <= fov_short_half_width;
sample_intensity = inside_long .* inside_fov .* ...
    exp(-2 .* (sample_s ./ gaussian_width).^2);
factor(valid) = mean(sample_intensity, 2);
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

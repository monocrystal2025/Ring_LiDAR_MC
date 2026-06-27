function [detect, first_detect] = MC_line_snr(fasan_D, fasan_d, path1, beam_vec, path2, f, omiga, UP, UV, init_beam_idx, R)
%MC_LINE_SNR  Monte Carlo detection judgment for a line beam.
%
% The long side keeps the legacy hard-edged flat-top profile.  Its full
% angular length is fasan_D and its direction is the same as MC_line.m.
% The short side is Gaussian.  Following the full-angle convention used by
% the legacy line-beam parameters, the intensity is exp(-2) at +/-fasan_d/2.
%
% Detection uses the same weighted variable sliding-window rule as
% MC_ring_snr, with SNR threshold 2.  Target partial illumination is sampled
% on a 7-by-7 grid in the local line-beam angular coordinates.

N_KERNAL = 100000;

SNR_THRESHOLD = 2;
TARGET_W = 0.297;             % m, kept consistent with MC_ring_snr
GAUSSIAN_CUTOFF_WIDTH = 3;     % short-side geometric prefilter: +/-3*fasan_d

if fasan_D <= 0
    error('MC_line_snr:InvalidLongWidth', ...
        'fasan_D must be positive for the line-beam long side.');
end
if fasan_d <= 0
    error('MC_line_snr:InvalidGaussianWidth', ...
        'fasan_d must be positive for the line-beam Gaussian short side.');
end
if omiga <= 0
    error('MC_line_snr:InvalidAngularSpeed', ...
        'omiga must be positive for sliding-window pulse accumulation.');
end

long_half_width = fasan_D / 2;
gaussian_width = fasan_d / 2;
short_prefilter_half_width = GAUSSIAN_CUTOFF_WIDTH * fasan_d;
WINDOW_PULSES = max(1, ceil(f * fasan_D / omiga));

photon_model = get_line_photon_model(fasan_D, fasan_d, gaussian_width, R);
carry_sig = zeros(0, 1);
carry_noise = zeros(0, 1);
carry_step = zeros(0, 1);

%% One complete beam cycle used by the original code
Beam_all = [beam_vec; path2; flipud(beam_vec); path1];
L_cycle = size(Beam_all, 1);

%% Include the initial pulse, which corresponds to first_detect = 0
idx0 = mod(init_beam_idx - 1, L_cycle) + 1;
B0 = Beam_all(idx0, :);
UP_norm = norm(UP);
UP_N = UP ./ UP_norm;

[long0, short0, front0] = line_offsets(UP_N, unitvec(B0));
illum0 = 0;
if front0 && is_line_gaussian_candidate(long0, short0, UP_norm, ...
        TARGET_W, long_half_width, short_prefilter_half_width)
    illum0 = line_gaussian_target_factor(long0, short0, UP_norm, TARGET_W, ...
        long_half_width, gaussian_width);
end
noise0 = lookup_line_noise(photon_model, UP_norm);
sig0 = calc_line_signal_photons(photon_model, UP_norm, illum0);
[detect, first_detect, carry_sig, carry_noise, carry_step] = ...
    append_and_check_window(carry_sig, carry_noise, carry_step, ...
    sig0, noise0, 0, WINDOW_PULSES, SNR_THRESHOLD);
if detect
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

    if block_len > 0
        P_block = P_now(1:block_len, :);
        B_block = B_now(1:block_len, :);
        z_block = norm_P_now(1:block_len);

        U_block = P_block ./ z_block;
        [long_offset, short_offset, front] = line_offsets(U_block, B_block);
        judge = front & is_line_gaussian_candidate(long_offset, short_offset, ...
            z_block, TARGET_W, long_half_width, short_prefilter_half_width);

        illum_arr = zeros(block_len, 1);
        hit_idx = find(judge);
        if ~isempty(hit_idx)
            illum_arr(hit_idx) = line_gaussian_target_factor( ...
                long_offset(hit_idx), short_offset(hit_idx), z_block(hit_idx), ...
                TARGET_W, long_half_width, gaussian_width);
        end

        noise_arr = lookup_line_noise(photon_model, z_block);
        sig_arr = zeros(block_len, 1);
        if ~isempty(hit_idx)
            sig_arr(hit_idx) = calc_line_signal_photons( ...
                photon_model, z_block(hit_idx), illum_arr(hit_idx));
        end

        step_arr = (block_id - 1) * N_KERNAL + (1:block_len)';
        [detect, first_detect, carry_sig, carry_noise, carry_step] = ...
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
    P_now(:,1) = P_now(:,1) + N_KERNAL * UV(1) / f;
    P_now(:,2) = P_now(:,2) + N_KERNAL * UV(2) / f;
    P_now(:,3) = P_now(:,3) + N_KERNAL * UV(3) / f;

    beam_idx = beam_idx + N_KERNAL;
    B_now = cut_loop_simple_anyk(path1, beam_vec, path2, beam_idx, N_KERNAL);
end
end

function model = get_line_photon_model(fasan_D, fasan_d, gaussian_width, R)
% Build and cache deterministic photon-count quantities for line beams.

persistent cached_model cached_key

key = sprintf('%.17g_%.17g_%.17g_%.17g', fasan_D, fasan_d, gaussian_width, R);
if ~isempty(cached_model) && strcmp(cached_key, key)
    model = cached_model;
    return;
end

%% Physical constants
h = 6.62607015e-34;               % J*s
c = 299792458;                    % m/s

%% System parameters, identical to MC_ring_snr.m
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

%% Target parameters, identical to MC_ring_snr.m
target_w = 0.297;                 % m
A_target = target_w * target_w;   % m^2
rho = 0.5;

E_photon = h * c / lambda;
N_tx = E_pulse / E_photon;
N_dark = dark_count_rate * tau_gate;

Omega_beam = fasan_D * gaussian_width * sqrt(pi / 2);
Omega_fov = fasan_D * (6 * fasan_d);

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

function N_noise = lookup_line_noise(model, z)
% Nearest-neighbor lookup on the 0.5 m deterministic backscatter table.

z = z(:);
idx = round(z ./ model.z_step) + 1;
idx(~isfinite(idx)) = 1;
idx = max(1, min(numel(model.N_bs_lut), idx));
N_noise = model.N_bs_lut(idx) + model.N_bg + model.N_dark;
end

function N_sig = calc_line_signal_photons(model, z, illumination_factor)
% Deterministic target-return photons for flat-long/Gaussian-short line light.

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

function [detect, first_detect, carry_sig, carry_noise, carry_step] = ...
    append_and_check_window(carry_sig, carry_noise, carry_step, ...
    sig_arr, noise_arr, step_arr, window_pulses, threshold)

detect = 0;
first_detect = nan;

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
    pulse_score(score_valid) = all_sig(score_valid).^2 ./ ...
        pulse_var(score_valid);
    cs_score = [0; cumsum(pulse_score)];

    win_score = cs_score(tail_pos + 1) - cs_score(start_pos);
    win_snr = sqrt(win_score);
    win_snr(~valid_len | win_score <= 0 | ~isfinite(win_snr)) = -inf;

    hit_matrix = win_snr >= threshold;
    ok = find(any(hit_matrix, 2), 1, 'first');
    if ~isempty(ok)
        detect = 1;
        first_detect = all_step(tail_pos(ok));
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

function [long_offset, short_offset, front] = line_offsets(U_dir, B_dir)
% Offsets in the legacy MC_line local frame.
% short_offset follows e1 = [-by, bx, 0], long_offset follows cross(B,e1).

U_dir = U_dir ./ vecnorm(U_dir, 2, 2);
B_dir = B_dir ./ vecnorm(B_dir, 2, 2);

bx = B_dir(:,1);
by = B_dir(:,2);
bxy = hypot(bx, by);

mask = (bxy < eps);
bxy_safe = bxy;
bxy_safe(mask) = 1;

e_short = [-by ./ bxy_safe, bx ./ bxy_safe, zeros(size(bx))];
if any(mask)
    e_short(mask,:) = repmat([1 0 0], nnz(mask), 1);
end

e_long = cross(B_dir, e_short, 2);

c3 = sum(U_dir .* B_dir, 2);
c_short = sum(U_dir .* e_short, 2);
c_long = sum(U_dir .* e_long, 2);

front = c3 > 0;
short_offset = atan2(c_short, c3);
long_offset = atan2(c_long, c3);
end

function candidate = is_line_gaussian_candidate(long_offset, short_offset, z, ...
    target_w, long_half_width, short_prefilter_half_width)
% Fast geometric prefilter. The exact signal comes from 7-by-7 target-area
% averaging, so the hard long edge is expanded only by half a target sample
% footprint. The Gaussian short side keeps the requested +/-3*fasan_d range.

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
    target_w, long_half_width, gaussian_width)
% Average normalized line-beam intensity over a square target.

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
sample_intensity = inside_long .* exp(-2 .* (sample_s ./ gaussian_width).^2);
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

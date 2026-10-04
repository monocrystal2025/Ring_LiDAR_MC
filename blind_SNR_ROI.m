clc;
clear;
close all;

%% Parameters
% Static blind-search Monte Carlo simulation inside the spherical cone ROI.
% Targets have uniform VOLUME density, remain stationary, and are observed
% for one complete ROI scan only. Detection uses the MC_EA photon/SNR rule.
% Dependencies (unchanged): build_ROI_scan_cycle, build_ROI_line_scan,
% mc_ea_block_size, append_and_check_photon_window, and its empty-window helper.
% The local ROI window search is equivalent to the shared MC_EA criterion;
% the shared implementation is retained as a nonfinite-arithmetic fallback.
N = 50000;
f = 5e3;
R = 2000;
JIAODU = 30;
fasan_D_LIST = 5e-3:5e-3:400e-3;
fasan_d = 1e-3;
omi_type = 'fast';
rng_seed = 1;

% Keep the original output directory and file names. Running this script
% replaces same-name results, just as static_MC_ROI_NEW does.
result_path = fullfile(pwd, 'STATIC_MC_ROI_RESULTS_NEW');
validate_parameters(N, f, R, JIAODU, fasan_D_LIST, fasan_d, ...
    omi_type);
if ~isfolder(result_path)
    mkdir(result_path);
end

switch lower(omi_type)
    case 'fast'
        omiga = 2 * pi;
    case 'slow'
        omiga = pi / 2;
end

% One population of N targets is reused for every beam, path, and angle.
% Each beam/path combination is simulated once; output column 2 is a copy
% for legacy file compatibility, not an independent Monte Carlo repeat.
rng(rng_seed, 'twister');
target_positions = sample_cone_volume(R, JIAODU, N);

detect_R_RA = zeros(N, 2);
detect_R_SP = zeros(N, 2);
detect_L_RA = zeros(N, 2);
detect_L_SP = zeros(N, 2);
detect_P_RA = zeros(N, 2);
detect_P_SP = zeros(N, 2);

% Keep the original N-by-2 double arrays and variable names. As in
% blind_snr_EA, first_time now means first successful SNR-window tail time
% in seconds (pulse 0 is time 0), and NaN denotes no detection in this scan.
first_time_R_RA = nan(N, 2);
first_time_R_SP = nan(N, 2);
first_time_L_RA = nan(N, 2);
first_time_L_SP = nan(N, 2);
first_time_P_RA = nan(N, 2);
first_time_P_SP = nan(N, 2);

path_types = {'raster', 'spiral'};

for D_index = numel(fasan_D_LIST):-1:1
    iteration_timer = tic;
    fasan_D = fasan_D_LIST(D_index);
    photon_models = build_photon_models(fasan_D, fasan_d, R, f, omiga);
    mc_block_size = mc_ea_block_size(fasan_D);

    % Exactly the same full ping-pong scan and spherical line-beam frames
    % as static_MC_ROI_NEW; do not concatenate another cycle.
    scan_cycles = cell(numel(path_types), 1);
    line_scans = cell(numel(path_types), 1);
    for path_index = 1:numel(path_types)
        scan_cycles{path_index} = build_ROI_scan_cycle( ...
            JIAODU, fasan_D, f, omiga, path_types{path_index});
        line_scans{path_index} = build_ROI_line_scan( ...
            scan_cycles{path_index});
    end

    % Exactly N trials per beam/path; duplicate only the saved columns.
    [detected, first_time] = stationary_snr_all_beams( ...
        target_positions, scan_cycles{1}, line_scans{1}, ...
        photon_models, f, mc_block_size);
    detect_R_RA(:,:) = repmat(detected(:, 1), 1, 2);
    detect_L_RA(:,:) = repmat(detected(:, 2), 1, 2);
    detect_P_RA(:,:) = repmat(detected(:, 3), 1, 2);
    first_time_R_RA(:,:) = repmat(first_time(:, 1), 1, 2);
    first_time_L_RA(:,:) = repmat(first_time(:, 2), 1, 2);
    first_time_P_RA(:,:) = repmat(first_time(:, 3), 1, 2);

    [detected, first_time] = stationary_snr_all_beams( ...
        target_positions, scan_cycles{2}, line_scans{2}, ...
        photon_models, f, mc_block_size);
    detect_R_SP(:,:) = repmat(detected(:, 1), 1, 2);
    detect_L_SP(:,:) = repmat(detected(:, 2), 1, 2);
    detect_P_SP(:,:) = repmat(detected(:, 3), 1, 2);
    first_time_R_SP(:,:) = repmat(first_time(:, 1), 1, 2);
    first_time_L_SP(:,:) = repmat(first_time(:, 2), 1, 2);
    first_time_P_SP(:,:) = repmat(first_time(:, 3), 1, 2);

    D_text = format_mrad(fasan_D);
    d_text = format_mrad(fasan_d);

    save(fullfile(result_path, ...
        ['STATIC_MC_1par_ROI_RING_RA_D', D_text, ...
        'd', d_text, 'mrad.mat']), ...
        'detect_R_RA', 'first_time_R_RA');
    save(fullfile(result_path, ...
        ['STATIC_MC_1par_ROI_RING_SP_D', D_text, ...
        'd', d_text, 'mrad.mat']), ...
        'detect_R_SP', 'first_time_R_SP');

    save(fullfile(result_path, ...
        ['STATIC_MC_1par_ROI_LINE_RA_D', D_text, ...
        'd', d_text, 'mrad.mat']), ...
        'detect_L_RA', 'first_time_L_RA');
    save(fullfile(result_path, ...
        ['STATIC_MC_1par_ROI_LINE_SP_D', D_text, ...
        'd', d_text, 'mrad.mat']), ...
        'detect_L_SP', 'first_time_L_SP');

    % Point-beam width depends on fasan_D alone, as in the original script.
    save(fullfile(result_path, ...
        ['STATIC_MC_1par_ROI_POINT_RA_D', D_text, 'mrad.mat']), ...
        'detect_P_RA', 'first_time_P_RA');
    save(fullfile(result_path, ...
        ['STATIC_MC_1par_ROI_POINT_SP_D', D_text, 'mrad.mat']), ...
        'detect_P_SP', 'first_time_P_SP');

    fprintf(['D = %s mrad | ring: RA %.6f, SP %.6f | ', ...
        'line: RA %.6f, SP %.6f | point: RA %.6f, SP %.6f ', ...
        '| %.2f s\n'], ...
        D_text, mean(detect_R_RA, 'all'), mean(detect_R_SP, 'all'), ...
        mean(detect_L_RA, 'all'), mean(detect_L_SP, 'all'), ...
        mean(detect_P_RA, 'all'), mean(detect_P_SP, 'all'), ...
        toc(iteration_timer));
end

function position = sample_cone_volume(R, jiaodu, N)
%SAMPLE_CONE_VOLUME Uniform volume density in 0 < range <= R, theta <= alpha.
% The ROI has a spherical outer cap, not a flat base at z = R.
% dV = r^2 dr d(cos(theta)) d(phi), hence radius = R * U^(1/3).

alpha = deg2rad(jiaodu);
azimuth = 2 * pi * rand(N, 1);
cos_polar = cos(alpha) + (1 - cos(alpha)) * rand(N, 1);
sin_polar = sqrt(max(1 - cos_polar.^2, 0));
radius = R * nthroot(max(rand(N, 1), realmin), 3);
position = radius .* [sin_polar .* cos(azimuth), ...
    sin_polar .* sin(azimuth), cos_polar];
end

function [detect, first_time] = stationary_snr_all_beams( ...
    targets, beam_cycle, line_scan, models, f, block_size)
%STATIONARY_SNR_ALL_BEAMS Observe each fixed target for a single ROI cycle.

N = size(targets, 1);
detect = zeros(N, 3);
first_pulse = nan(N, 3);
use_parallel = license('test', 'Distrib_Computing_Toolbox') && ...
    ~isempty(ver('parallel'));
if use_parallel
    pool = gcp('nocreate');
    if isempty(pool)
        try
            parpool('threads');
        catch
            parpool('local');
        end
    end
    parfor target_index = 1:N
        [detect(target_index, :), first_pulse(target_index, :)] = ...
            scan_one_target(targets(target_index, :), beam_cycle, ...
            line_scan, models, block_size);
    end
else
    for target_index = 1:N
        [detect(target_index, :), first_pulse(target_index, :)] = ...
            scan_one_target(targets(target_index, :), beam_cycle, ...
            line_scan, models, block_size);
    end
end
first_time = first_pulse ./ f;
end

function [detect, first_pulse] = scan_one_target( ...
    position, beam_cycle, line_scan, models, block_size)
%SCAN_ONE_TARGET Causal SNR windows over pulse indices 0 through L-1 only.

detect = zeros(1, 3);
first_pulse = nan(1, 3);
target_range = norm(position);
direction = position ./ target_range;
carry_sig = repmat({zeros(0, 1)}, 1, 3);
carry_backscatter = repmat({zeros(0, 1)}, 1, 3);
carry_step = repmat({zeros(0, 1)}, 1, 3);

% Position never changes: this is exactly velocity [0, 0, 0]. Backscatter
% is looked up on the same nearest-neighbour 0.5 m table as MC_EA.
noise_index = round(target_range / models(1).z_step) + 1;
noise_index = max(1, min(numel(models(1).N_bs_lut), noise_index));
backscatter = models(1).N_bs_lut(noise_index);
beam_count = size(beam_cycle, 1);

% No modulo indexing, second cycle, post-detection extension, or padding.
% The initial pulse is included exactly once as row 1 / pulse index 0.
for first_beam = 1:block_size:beam_count
    last_beam = min(first_beam + block_size - 1, beam_count);
    beam_index = first_beam:last_beam;
    steps = (beam_index - 1).';
    beams = beam_cycle(beam_index, :);
    for beam_type = 1:3
        if detect(beam_type)
            continue;
        end
        model = models(beam_type);
        illumination = target_illumination(direction, target_range, ...
            beams, line_scan, beam_index, model);
        signal = signal_photons(model, target_range, illumination);
        [detect(beam_type), first_pulse(beam_type), ~, ...
            carry_sig{beam_type}, carry_backscatter{beam_type}, ...
            carry_step{beam_type}] = append_and_check_roi_window( ...
            carry_sig{beam_type}, carry_backscatter{beam_type}, ...
            carry_step{beam_type}, signal, ...
            repmat(backscatter, numel(steps), 1), steps, ...
            model.N_bg, model.N_dark, model.window_pulses, 2);
    end
    % Further pulses cannot change an already known first detection.
    if all(detect)
        break;
    end
end
end

function [detect, first_detect, window_photons, ...
    carry_sig, carry_backscatter, carry_step] = ...
    append_and_check_roi_window( ...
    carry_sig, carry_backscatter, carry_step, ...
    sig_arr, backscatter_arr, step_arr, background_photons, ...
    dark_count_photons, window_pulses, threshold)
%APPEND_AND_CHECK_ROI_WINDOW Exact longest-window screening of causal SNR.
% Each per-pulse score is nonnegative. With finite cumulative scores, the
% longest valid window at a given tail passes iff some valid window passes.
% Screen tails in time order, then enumerate lengths only at the first hit.
% Keep the shared function's arithmetic order, sqrt/threshold comparison,
% shortest-window selection, photon outputs, and cross-block carry unchanged.

detect = 0;
first_detect = nan;
window_photons = empty_detection_window_photons();

sig_arr = sig_arr(:);
backscatter_arr = backscatter_arr(:);
step_arr = step_arr(:);
if isempty(sig_arr)
    return;
end

sig_arr(~isfinite(sig_arr)) = 0;
backscatter_arr(~isfinite(backscatter_arr)) = 0;

n_carry = numel(carry_sig);
all_sig = [carry_sig; sig_arr];
all_backscatter = [carry_backscatter; backscatter_arr];
all_step = [carry_step; step_arr];
n_all = numel(all_sig);

signal_tail_pos = n_carry + find(sig_arr > 0);
if ~isempty(signal_tail_pos)
    tail_pos = signal_tail_pos(:);
    max_len = min(window_pulses, n_all);

    all_noise = all_backscatter + background_photons + dark_count_photons;
    pulse_var = all_sig + all_noise;
    per_pulse_score = zeros(size(all_sig));
    score_valid = all_sig > 0 & pulse_var > 0 & isfinite(pulse_var);
    per_pulse_score(score_valid) = all_sig(score_valid).^2 ./ ...
        pulse_var(score_valid);
    cumulative_score = [0; cumsum(per_pulse_score)];

    % Inf/NaN breaks the finite-score argument: the reference rejects
    % nonfinite SNR. Preserve its behavior by using the reference directly.
    if any(~isfinite(cumulative_score))
        [detect, first_detect, window_photons, ...
            carry_sig, carry_backscatter, carry_step] = ...
            append_and_check_photon_window( ...
            carry_sig, carry_backscatter, carry_step, ...
            sig_arr, backscatter_arr, step_arr, background_photons, ...
            dark_count_photons, window_pulses, threshold);
        return;
    end

    longest_start = max(1, tail_pos - max_len + 1);
    longest_score = cumulative_score(tail_pos + 1) - ...
        cumulative_score(longest_start);
    longest_snr = sqrt(longest_score);
    longest_snr(longest_score <= 0 | ~isfinite(longest_snr)) = -inf;
    first_successful_tail = find(longest_snr >= threshold, 1, 'first');

    if ~isempty(first_successful_tail)
        successful_tail = tail_pos(first_successful_tail);
        len_vec = 1:max_len;
        start_pos = successful_tail - len_vec + 1;
        valid_len = start_pos >= 1;
        start_pos(~valid_len) = 1;
        start_score = reshape(cumulative_score(start_pos(:)), ...
            size(start_pos));
        win_score = cumulative_score(successful_tail + 1) - start_score;
        win_snr = sqrt(win_score);
        win_snr(~valid_len | win_score <= 0 | ~isfinite(win_snr)) = -inf;
        successful_len = find(win_snr >= threshold, 1, 'first');
        selected = start_pos(successful_len):successful_tail;

        detect = 1;
        first_detect = all_step(successful_tail);
        window_photons.signal_photons = all_sig(selected);
        window_photons.backscatter_photons = all_backscatter(selected);
        window_photons.background_photons = background_photons;
        window_photons.dark_count_photons = dark_count_photons;
        window_photons.pulse_steps = uint32(all_step(selected));
        return;
    end
end

n_keep = min(window_pulses - 1, n_all);
if n_keep > 0
    carry_sig = all_sig(end - n_keep + 1:end);
    carry_backscatter = all_backscatter(end - n_keep + 1:end);
    carry_step = all_step(end - n_keep + 1:end);
else
    carry_sig = zeros(0, 1);
    carry_backscatter = zeros(0, 1);
    carry_step = zeros(0, 1);
end
end

function factor = target_illumination(direction, target_range, ...
    beams, line_scan, beam_index, model)
%TARGET_ILLUMINATION MC_EA finite-target Gaussian model in ROI scan frames.
% Preserve MC_EA prefilters and 7-by-7 target-area averaging, including its
% untruncated ring/line Gaussian factors and hard point-beam receiver FOV.
% Only the line orientation changes to the ROI spherical tangent frame.

persistent sample_x sample_y
if isempty(sample_x)
    nodes = ((1:7) - 0.5) / 7 - 0.5;
    [xx, yy] = meshgrid(nodes, nodes);
    sample_x = xx(:).';
    sample_y = yy(:).';
end

factor = zeros(size(beams, 1), 1);
cos_theta = beams * direction.';
angular_side = 2 * atan((model.target_w / 2) / target_range);
dx = angular_side .* sample_x;
dy = angular_side .* sample_y;
half_diag = atan((sqrt(2) * model.target_w / 2) / target_range);

switch model.beam_type
    case 'ring'
        theta = acos(max(min(cos_theta, 1), -1));
        candidate = abs(theta - model.long_half_width) <= ...
            half_diag + 3 * model.gaussian_width;
        if ~any(candidate)
            return;
        end
        sample_theta = sqrt((theta(candidate) + dx).^2 + dy.^2);
        intensity = exp(-2 .* ((sample_theta - model.long_half_width) ...
            ./ model.gaussian_width).^2);
    case 'line'
        long_component = line_scan.LongAxes(beam_index, :) * direction.';
        short_component = line_scan.ShortAxes(beam_index, :) * direction.';
        long_offset = atan2(long_component, cos_theta);
        short_offset = atan2(short_component, cos_theta);
        half_side = atan((model.target_w / 2) / target_range);
        candidate = cos_theta > 0 & ...
            abs(long_offset) <= model.long_half_width + half_side & ...
            abs(short_offset) <= model.fov_half_width + half_side;
        if ~any(candidate)
            return;
        end
        sample_long = long_offset(candidate) + dx;
        sample_short = short_offset(candidate) + dy;
        intensity = (abs(sample_long) <= model.long_half_width) .* ...
            exp(-2 .* (sample_short ./ model.gaussian_width).^2);
    case 'point'
        theta_limit = min(model.fov_half_width + half_diag, pi);
        candidate = cos_theta >= cos(theta_limit);
        if ~any(candidate)
            return;
        end
        theta = acos(max(min(cos_theta(candidate), 1), -1));
        sample_theta = sqrt((theta + dx).^2 + dy.^2);
        intensity = (sample_theta <= model.fov_half_width) .* ...
            exp(-2 .* (sample_theta ./ model.gaussian_width).^2);
    otherwise
        error('blind_SNR_ROI:UnknownBeamType', ...
            'Unknown beam type: %s.', model.beam_type);
end
factor(candidate) = mean(intensity, 2);
end

function signal = signal_photons(model, target_range, illumination)
%SIGNAL_PHOTONS Deterministic MC_EA target return, including two-way loss.

signal = zeros(size(illumination));
valid = illumination > 0 & isfinite(illumination);
beam_area = model.Omega_beam * target_range^2;
signal(valid) = model.N_tx .* ...
    (model.A_target .* illumination(valid) ./ beam_area) .* ...
    model.rho .* (model.A_rx ./ (pi .* target_range^2)) .* ...
    exp(-2 * model.alpha * target_range) .* model.eta_sys;
end

function models = build_photon_models(fasan_D, fasan_d, R, f, omiga)
%BUILD_PHOTON_MODELS Same radiometry and window sizes as MC_*_snr in MC_EA.
% SNR = sqrt(sum(N_sig.^2 ./ (N_sig + N_bs + N_bg + N_dark))) >= 2.
% Missed pulses retain their time slots and noise; partial startup windows
% and cross-block carry are handled by append_and_check_roi_window.

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
rho = 0.5;
E_photon = h * c / lambda;
N_tx = E_pulse / E_photon;

z_step = 0.5;
z_lut = (0:z_step:R).';
if z_lut(end) < R
    z_lut = [z_lut; R];
end
dr = c * tau_gate / 2;
offsets = linspace(-dr / 2, dr / 2, 101);
ranges = z_lut + offsets;
ranges(ranges <= 0) = 1e-6;
integrand = beta .* exp(-2 .* alpha .* ranges) ./ ranges.^2;
N_bs_lut = N_tx .* eta_sys .* A_rx .* trapz(offsets, integrand, 2);

common = struct('z_step', z_step, 'N_bs_lut', N_bs_lut, ...
    'N_bg', 0, 'N_dark', dark_count_rate * tau_gate, ...
    'N_tx', N_tx, 'A_target', target_w * target_w, 'rho', rho, ...
    'A_rx', A_rx, 'eta_sys', eta_sys, 'alpha', alpha, ...
    'target_w', target_w, 'Omega_beam', 0, 'beam_type', '', ...
    'long_half_width', fasan_D / 2, 'gaussian_width', fasan_d / 2, ...
    'fov_half_width', 3 * fasan_d, ...
    'window_pulses', max(1, ceil(f * fasan_D / omiga)));
models = repmat(common, 1, 3);

models(1).beam_type = 'ring';
ring_radius = fasan_D / 2;
width = fasan_d / 2;
u0 = -sqrt(2) * ring_radius / width;
radial_integral = ring_radius * width * sqrt(pi) / (2 * sqrt(2)) * ...
    erfc(u0) + width^2 / 4 * exp(-2 * (ring_radius / width)^2);
models(1).Omega_beam = 2 * pi * radial_integral;
inner = max(ring_radius - 3 * fasan_d, 0);
outer = min(ring_radius + 3 * fasan_d, pi);
fov_solid_angle = zeros(1, 3);
fov_solid_angle(1) = 2 * pi * (cos(inner) - cos(outer));

models(2).beam_type = 'line';
models(2).Omega_beam = fasan_D * width * sqrt(pi / 2);
fov_solid_angle(2) = fasan_D * (6 * fasan_d);

models(3).beam_type = 'point';
models(3).gaussian_width = fasan_D / 2;
models(3).fov_half_width = 3 * fasan_D / 2;
models(3).window_pulses = max(1, ceil(f * (3 * fasan_D) / omiga));
point_integrand = @(theta) exp(-2 .* (theta ./ (fasan_D / 2)).^2) ...
    .* sin(theta);
models(3).Omega_beam = 2 * pi * integral(point_integrand, 0, pi, ...
    'RelTol', 1e-10, 'AbsTol', 1e-14);
fov_solid_angle(3) = 2 * pi * (1 - cos(min(3 * fasan_D / 2, pi)));

for beam_type = 1:3
    P_bg = L_sky_nm * delta_lambda_nm * A_rx * fov_solid_angle(beam_type);
    models(beam_type).N_bg = P_bg * tau_gate * eta_sys / E_photon;
end
end

function validate_parameters(N, f, R, jiaodu, fasan_D_list, ...
    fasan_d, omi_type)
%VALIDATE_PARAMETERS Check physical inputs before creating result files.

validateattributes(N, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'positive'}, mfilename, 'N');
validateattributes(f, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, mfilename, 'f');
validateattributes(R, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, mfilename, 'R');
validateattributes(jiaodu, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>', 0, '<', 90}, mfilename, 'JIAODU');
validateattributes(fasan_D_list, {'numeric'}, ...
    {'vector', 'nonempty', 'real', 'finite', '>', 0, '<', pi}, ...
    mfilename, 'fasan_D_LIST');
validateattributes(fasan_d, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>', 0, '<', pi}, mfilename, 'fasan_d');
if ~any(strcmpi(omi_type, {'fast', 'slow'}))
    error('blind_SNR_ROI:UnknownOmiType', ...
        'omi_type must be fast or slow.');
end
end

function text = format_mrad(angle)
%FORMAT_MRAD Match the numeric filename convention of static_MC_ROI_NEW.

text = num2str(angle * 1000, '%g');
end

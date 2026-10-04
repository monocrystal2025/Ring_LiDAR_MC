clc;
clear;
close all;

%% Parameters
% Static blind-search Monte Carlo simulation using the SNR sliding-window
% detection rule implemented by MC_ring_snr, MC_line_snr, and MC_point_snr.
% Targets are uniformly distributed throughout the upper half-ball, remain
% stationary, and each trial is limited to one complete scan trajectory.
N = 50000;
f = 5e3;
omi_type = 'fast';
fasan_D_LIST = 5e-3:5e-3:400e-3;
fasan_d = 1e-3;
R = 2000;

beam_path_dir = 'G:\BeamVEC_NEW';
result_path = 'D:\matlab\STATIC_MC_RESULTS';

if ~isfolder(result_path)
    mkdir(result_path);
end

% One common target population is reused by all beam sizes and beam types.
% A cube-root radial law gives uniform volume density (not uniform radius).
target_position = sample_uniform_upper_half_ball(N, R);
stationary_velocity = [0, 0, 0];

detect_R = zeros(N, 2);
detect_L = zeros(N, 2);
detect_P = zeros(N, 2);
first_time_R = nan(N, 2);
first_time_L = nan(N, 2);
first_time_P = nan(N, 2);

for i = numel(fasan_D_LIST):-1:1
    iteration_timer = tic;
    fasan_D = fasan_D_LIST(i);
    mc_block_size = mc_ea_block_size(fasan_D);

    [path1, beam_vec, path2, omega] = load_complete_scan_parts( ...
        beam_path_dir, f, fasan_D, omi_type);
    scan_pulse_count = size(beam_vec, 1) + size(path2, 1) + ...
        size(beam_vec, 1) + size(path1, 1);

    % Pulse index 0 is evaluated separately inside each MC_*_snr function.
    % Therefore L scan pulses correspond to a maximum step index of L - 1.
    max_step_index = scan_pulse_count - 1;

    [detect_one, first_one] = run_stationary_snr( ...
        target_position, stationary_velocity, 'ring', ...
        fasan_D, fasan_d, path1, beam_vec, path2, f, omega, ...
        R, mc_block_size, max_step_index);
    detect_R(:,:) = repmat(detect_one, 1, 2);
    first_time_R(:,:) = repmat(first_one, 1, 2);
    save(fullfile(result_path, ...
        ['STATIC_MC_1par_EA_RING_D', format_mrad(fasan_D), ...
        'd', format_mrad(fasan_d), 'mrad.mat']), ...
        'detect_R', 'first_time_R');

    [detect_one, first_one] = run_stationary_snr( ...
        target_position, stationary_velocity, 'line', ...
        fasan_D, fasan_d, path1, beam_vec, path2, f, omega, ...
        R, mc_block_size, max_step_index);
    detect_L(:,:) = repmat(detect_one, 1, 2);
    first_time_L(:,:) = repmat(first_one, 1, 2);
    save(fullfile(result_path, ...
        ['STATIC_MC_1par_EA_LINE_D', format_mrad(fasan_D), ...
        'd', format_mrad(fasan_d), 'mrad.mat']), ...
        'detect_L', 'first_time_L');

    [detect_one, first_one] = run_stationary_snr( ...
        target_position, stationary_velocity, 'point', ...
        fasan_D, [], path1, beam_vec, path2, f, omega, ...
        R, mc_block_size, max_step_index);
    detect_P(:,:) = repmat(detect_one, 1, 2);
    first_time_P(:,:) = repmat(first_one, 1, 2);
    save(fullfile(result_path, ...
        ['STATIC_MC_1par_EA_POINT_D', format_mrad(fasan_D), ...
        'mrad.mat']), ...
        'detect_P', 'first_time_P');

    fprintf(['D = %s mrad: ring = %.6f, line = %.6f, ', ...
        'point = %.6f (%.2f s)\n'], ...
        format_mrad(fasan_D), mean(detect_R(:, 1)), ...
        mean(detect_L(:, 1)), mean(detect_P(:, 1)), ...
        toc(iteration_timer));
end

function target_position = sample_uniform_upper_half_ball(N, R)
%SAMPLE_UNIFORM_UPPER_HALF_BALL Uniform points in a radius-R upper half-ball.

azimuth = 2 * pi * rand(N, 1);
cos_polar = rand(N, 1);
sin_polar = sqrt(1 - cos_polar.^2);
radius = R * nthroot(max(rand(N, 1), realmin), 3);

target_position = radius .* [ ...
    sin_polar .* cos(azimuth), ...
    sin_polar .* sin(azimuth), ...
    cos_polar];
end

function [path1, beam_vec, path2, omega] = load_complete_scan_parts( ...
    beam_path_dir, f, fasan_D, omi_type)
%LOAD_COMPLETE_SCAN_PARTS Load and reproduce the MC_EA complete scan.

switch lower(omi_type)
    case 'fast'
        omega = 2 * pi;
    case 'slow'
        omega = pi / 2;
    otherwise
        error('blind_snr_EA:UnknownOmiType', ...
            'Unknown omi_type: %s. Use fast or slow.', omi_type);
end

file_name = ['beam_vector_', omi_type, '_', num2str(f), ...
    'Hz_h3000m_D', num2str(fasan_D * 1e6), 'u.npy'];
file_path = fullfile(beam_path_dir, file_name);
if ~isfile(file_path)
    error('blind_snr_EA:MissingBeamPath', ...
        'Beam-path file does not exist: %s', file_path);
end

beam_vec = readNPY(file_path);
if size(beam_vec, 2) ~= 3 || isempty(beam_vec)
    error('blind_snr_EA:InvalidBeamPath', ...
        'Beam-path data must be a nonempty N-by-3 array: %s', file_path);
end

path1 = complete_boundary_circle(beam_vec(1, :), omega, f);
path2 = complete_boundary_circle(beam_vec(end, :), omega, f);
end

function path = complete_boundary_circle(start_vector, omega, f)
%COMPLETE_BOUNDARY_CIRCLE Complete one boundary turn at scan angular speed.

x0 = start_vector(1);
y0 = start_vector(2);
z0 = start_vector(3);
theta0 = atan2(y0, x0);
radius_xy = hypot(x0, y0);
if radius_xy <= eps
    error('blind_snr_EA:DegenerateBoundary', ...
        'A scan boundary cannot be completed at the polar axis.');
end

theta = (theta0:omega / (f * radius_xy):theta0 + 2 * pi).';
path = [radius_xy * cos(theta), radius_xy * sin(theta), ...
    repmat(z0, numel(theta), 1)];
end

function [detect, first_time] = run_stationary_snr( ...
    target_position, stationary_velocity, beam_type, ...
    fasan_D, fasan_d, path1, beam_vec, path2, f, omega, ...
    R, mc_block_size, max_step_index)
%RUN_STATIONARY_SNR Apply the MC_EA SNR rule for one complete static scan.

N = size(target_position, 1);
detect = zeros(N, 1);
first_pulse = nan(N, 1);

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
end

if use_parallel
    parfor target_index = 1:N
        [detect(target_index), first_pulse(target_index)] = ...
            evaluate_one_target(target_position(target_index, :), ...
            stationary_velocity, beam_type, fasan_D, fasan_d, ...
            path1, beam_vec, path2, f, omega, R, mc_block_size, ...
            max_step_index);
    end
else
    for target_index = 1:N
        [detect(target_index), first_pulse(target_index)] = ...
            evaluate_one_target(target_position(target_index, :), ...
            stationary_velocity, beam_type, fasan_D, fasan_d, ...
            path1, beam_vec, path2, f, omega, R, mc_block_size, ...
            max_step_index);
    end
end

first_time = first_pulse ./ f;
end

function [detect, first_pulse] = evaluate_one_target( ...
    target_position, stationary_velocity, beam_type, ...
    fasan_D, fasan_d, path1, beam_vec, path2, f, omega, ...
    R, mc_block_size, max_step_index)
%EVALUATE_ONE_TARGET Dispatch to the unchanged MC_EA beam/SNR model.

initial_beam_index = 1;
switch lower(beam_type)
    case 'ring'
        [detect, first_pulse] = MC_ring_snr( ...
            fasan_D, fasan_d, path1, beam_vec, path2, f, omega, ...
            target_position, stationary_velocity, initial_beam_index, ...
            R, mc_block_size, max_step_index);
    case 'line'
        [detect, first_pulse] = MC_line_snr( ...
            fasan_D, fasan_d, path1, beam_vec, path2, f, omega, ...
            target_position, stationary_velocity, initial_beam_index, ...
            R, mc_block_size, max_step_index);
    case 'point'
        [detect, first_pulse] = MC_point_snr( ...
            fasan_D, fasan_d, path1, beam_vec, path2, f, omega, ...
            target_position, stationary_velocity, initial_beam_index, ...
            R, mc_block_size, max_step_index);
    otherwise
        error('blind_snr_EA:UnknownBeamType', ...
            'Unknown beam type: %s.', beam_type);
end
end

function text = format_mrad(angle)
%FORMAT_MRAD Match the numeric filename convention of the reference script.

text = num2str(angle * 1000, '%g');
end

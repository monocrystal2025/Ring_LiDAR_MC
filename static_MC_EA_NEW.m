clc;
clear;
close all;

%% Parameters
% Static spatial-coverage Monte Carlo simulation for the three beam types
% defined by MC_EA.m. A target is covered when one complete scan brings it
% inside the exp(-2) beam footprint. No SNR or receiver-FOV test is used.
N = 50000;
f = 5e3;
omi_type = 'fast';
fasan_D_LIST = 5e-3:5e-3:500e-3;
fasan_d = 1e-3;
R = 2000;

beam_path_dir = 'G:\BeamVEC_NEW';
result_path = 'D:\matlab\STATIC_MC_RESULTS';

if ~isfolder(result_path)
    mkdir(result_path);
end

% Use the same uniformly distributed static targets on the R-radius upper
% hemisphere for all cases. Only their directions enter the angular test.
target_position = R * sample_uniform_hemisphere(N);
target_direction = target_position ./ vecnorm(target_position, 2, 2);

detect_R = zeros(N, 2);
detect_L = zeros(N, 2);
detect_P = zeros(N, 2);
first_time_R = zeros(N, 2);
first_time_L = zeros(N, 2);
first_time_P = zeros(N, 2);

for i = numel(fasan_D_LIST):-1:1
    iteration_timer = tic;
    fasan_D = fasan_D_LIST(i);

    beam_cycle = load_complete_scan_cycle( ...
        beam_path_dir, f, fasan_D, omi_type);

    % Ring: Gaussian annulus centred at fasan_D/2. fasan_d is the full
    % radial exp(-2) width, so the included band is +/-fasan_d/2.
    detect_R(:, 1) = static_coverage( ...
        target_direction, beam_cycle, 'ring', fasan_D, fasan_d);
    detect_R(:, 2) = detect_R(:, 1);
    first_time_R(:) = 1;
    save(fullfile(result_path, ...
        ['STATIC_MC_1par_EA_RING_D', format_mrad(fasan_D), ...
        'd', format_mrad(fasan_d), 'mrad.mat']), ...
        'detect_R', 'first_time_R');

    % Line: hard-edged long side fasan_D and exp(-2) short side fasan_d.
    detect_L(:, 1) = static_coverage( ...
        target_direction, beam_cycle, 'line', fasan_D, fasan_d);
    detect_L(:, 2) = detect_L(:, 1);
    first_time_L(:) = 1;
    save(fullfile(result_path, ...
        ['STATIC_MC_1par_EA_LINE_D', format_mrad(fasan_D), ...
        'd', format_mrad(fasan_d), 'mrad.mat']), ...
        'detect_L', 'first_time_L');

    % Point: MC_EA.m defines fasan_D as the circular exp(-2) full angle;
    % fasan_d is intentionally ignored for this beam type.
    detect_P(:, 1) = static_coverage( ...
        target_direction, beam_cycle, 'point', fasan_D, []);
    detect_P(:, 2) = detect_P(:, 1);
    first_time_P(:) = 1;
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

function target_direction = sample_uniform_hemisphere(N)
%SAMPLE_UNIFORM_HEMISPHERE Uniform unit directions over the upper hemisphere.

azimuth = 2 * pi * rand(N, 1);
cos_polar = rand(N, 1);
sin_polar = sqrt(1 - cos_polar.^2);
target_direction = [sin_polar .* cos(azimuth), ...
    sin_polar .* sin(azimuth), cos_polar];
end

function beam_cycle = load_complete_scan_cycle( ...
    beam_path_dir, f, fasan_D, omi_type)
%LOAD_COMPLETE_SCAN_CYCLE Reproduce the complete scan used by MC_func.m.

switch lower(omi_type)
    case 'fast'
        omega = 2 * pi;
    case 'slow'
        omega = pi / 2;
    otherwise
        error('static_MC_EA_NEW:UnknownOmiType', ...
            'Unknown omi_type: %s. Use fast or slow.', omi_type);
end

file_name = ['beam_vector_', omi_type, '_', num2str(f), ...
    'Hz_h3000m_D', num2str(fasan_D * 1e6), 'u.npy'];
file_path = fullfile(beam_path_dir, file_name);
if ~isfile(file_path)
    error('static_MC_EA_NEW:MissingBeamPath', ...
        'Beam-path file does not exist: %s', file_path);
end

beam_vec = readNPY(file_path);
if size(beam_vec, 2) ~= 3 || isempty(beam_vec)
    error('static_MC_EA_NEW:InvalidBeamPath', ...
        'Beam-path data must be a nonempty N-by-3 array: %s', file_path);
end

path1 = complete_boundary_circle(beam_vec(1, :), omega, f);
path2 = complete_boundary_circle(beam_vec(end, :), omega, f);

% This is one complete scan in the same order as MC_func/MC_*.
beam_cycle = [beam_vec; path2; flipud(beam_vec); path1];
beam_cycle = beam_cycle ./ vecnorm(beam_cycle, 2, 2);
end

function path = complete_boundary_circle(start_vector, omega, f)
%COMPLETE_BOUNDARY_CIRCLE Complete one boundary turn at scan angular speed.

x0 = start_vector(1);
y0 = start_vector(2);
z0 = start_vector(3);
theta0 = atan2(y0, x0);
radius_xy = hypot(x0, y0);
if radius_xy <= eps
    error('static_MC_EA_NEW:DegenerateBoundary', ...
        'A scan boundary cannot be completed at the polar axis.');
end

theta = (theta0:omega / (f * radius_xy):theta0 + 2 * pi).';
path = [radius_xy * cos(theta), radius_xy * sin(theta), ...
    repmat(z0, numel(theta), 1)];
end

function detect = static_coverage( ...
    target_direction, beam_cycle, beam_type, fasan_D, fasan_d)
%STATIC_COVERAGE Test static targets against one complete scan geometrically.

N = size(target_direction, 1);
target_block_size = 250;
beam_block_size = 10000;

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

block_starts = 1:target_block_size:N;
target_blocks = cell(numel(block_starts), 1);
detected_blocks = cell(numel(block_starts), 1);
for block_index = 1:numel(block_starts)
    first_target = block_starts(block_index);
    last_target = min(first_target + target_block_size - 1, N);
    target_blocks{block_index} = ...
        target_direction(first_target:last_target, :);
end

if use_parallel
    parfor block_index = 1:numel(block_starts)
        detected_blocks{block_index} = test_target_block( ...
            target_blocks{block_index}, beam_cycle, ...
            beam_type, fasan_D, fasan_d, beam_block_size);
    end
else
    for block_index = 1:numel(block_starts)
        detected_blocks{block_index} = test_target_block( ...
            target_blocks{block_index}, beam_cycle, ...
            beam_type, fasan_D, fasan_d, beam_block_size);
    end
end

detect = double(vertcat(detected_blocks{:}));
end

function detected = test_target_block( ...
    targets, beam_cycle, beam_type, fasan_D, fasan_d, beam_block_size)
%TEST_TARGET_BLOCK Stop testing each target after its first geometric hit.

detected = false(size(targets, 1), 1);
for first_beam = 1:beam_block_size:size(beam_cycle, 1)
    active = find(~detected);
    if isempty(active)
        break;
    end

    last_beam = min(first_beam + beam_block_size - 1, ...
        size(beam_cycle, 1));
    beam_block = beam_cycle(first_beam:last_beam, :);
    target_block = targets(active, :);

    switch lower(beam_type)
        case 'ring'
            cos_angle = target_block * beam_block.';
            cos_angle = max(min(cos_angle, 1), -1);
            angle = acos(cos_angle);
            radial_offset = abs(angle - fasan_D / 2);
            hit = any(radial_offset <= fasan_d / 2, 2);
        case 'point'
            cos_angle = target_block * beam_block.';
            hit = any(cos_angle >= cos(fasan_D / 2), 2);
        case 'line'
            hit = line_beam_hit(target_block, beam_block, ...
                fasan_D, fasan_d);
        otherwise
            error('static_MC_EA_NEW:UnknownBeamType', ...
                'Unknown beam type: %s.', beam_type);
    end

    detected(active(hit)) = true;
end
end

function hit = line_beam_hit(targets, beams, fasan_D, fasan_d)
%LINE_BEAM_HIT Rectangular exp(-2) footprint in the MC_line local frame.

bx = beams(:, 1);
by = beams(:, 2);
horizontal_norm = hypot(bx, by);
near_axis = horizontal_norm < eps;
horizontal_norm(near_axis) = 1;

short_axis = [-by ./ horizontal_norm, bx ./ horizontal_norm, ...
    zeros(size(bx))];
short_axis(near_axis, :) = repmat([1, 0, 0], nnz(near_axis), 1);
long_axis = cross(beams, short_axis, 2);

along_beam = targets * beams.';
short_component = targets * short_axis.';
long_component = targets * long_axis.';

inside = along_beam > 0 & ...
    abs(long_component) <= along_beam * tan(fasan_D / 2) & ...
    abs(short_component) <= along_beam * tan(fasan_d / 2);
hit = any(inside, 2);
end

function text = format_mrad(angle)
%FORMAT_MRAD Match the numeric filename convention of the legacy script.

text = num2str(angle * 1000, '%g');
end

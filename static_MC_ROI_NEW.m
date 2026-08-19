clc;
clear;
close all;

%% Parameters
% Static 3-D spatial-coverage Monte Carlo simulation for the cone ROI.
% A target is covered if one complete scan places it inside the exp(-2)
% footprint of the selected beam.  No SNR or receiver-FOV test is applied.
N = 50000;
f = 5e3;
R = 2000;
JIAODU = 30;
fasan_D_LIST = 5e-3:5e-3:500e-3;
fasan_d = 1e-3;
omi_type = 'fast';
repeat_count = 2;
rng_seed = 1;

result_path = fullfile(pwd, 'STATIC_MC_ROI_RESULTS_NEW');
if ~isfolder(result_path)
    mkdir(result_path);
end

validate_parameters(N, f, R, JIAODU, fasan_D_LIST, fasan_d, ...
    omi_type, repeat_count);

switch lower(omi_type)
    case 'fast'
        omiga = 2 * pi;
    case 'slow'
        omiga = pi / 2;
end

% Use two independent Monte Carlo target populations for the two columns
% used by static_MC_ROI.m.  Reusing each population for all cases gives a
% lower-variance comparison between beam types, paths, and fasan_D values.
rng(rng_seed, 'twister');
target_direction = cell(repeat_count, 1);
for repeat_index = 1:repeat_count
    target_position = sample_cone_end_cap(R, JIAODU, N);
    target_direction{repeat_index} = target_position ./ ...
        vecnorm(target_position, 2, 2);
end

detect_R_RA = zeros(N, repeat_count);
detect_R_SP = zeros(N, repeat_count);
detect_L_RA = zeros(N, repeat_count);
detect_L_SP = zeros(N, repeat_count);
detect_P_RA = zeros(N, repeat_count);
detect_P_SP = zeros(N, repeat_count);

% The legacy static function returns the scalar value 1 for first_time;
% preserve the resulting N-by-2 MAT-file structure exactly.
first_time_R_RA = ones(N, repeat_count);
first_time_R_SP = ones(N, repeat_count);
first_time_L_RA = ones(N, repeat_count);
first_time_L_SP = ones(N, repeat_count);
first_time_P_RA = ones(N, repeat_count);
first_time_P_SP = ones(N, repeat_count);

path_types = {'raster', 'spiral'};

for D_index = numel(fasan_D_LIST):-1:1
    iteration_timer = tic;
    fasan_D = fasan_D_LIST(D_index);

    % Build exactly the spherical scan cycles and line-beam frames used by
    % MC_func_ROI_NEW.  Each returned beam_cycle is one complete ping-pong
    % scan cycle sampled at the physical angular step omiga/f.
    scan_cycles = cell(numel(path_types), 1);
    line_scans = cell(numel(path_types), 1);
    for path_index = 1:numel(path_types)
        scan_cycles{path_index} = build_ROI_scan_cycle( ...
            JIAODU, fasan_D, f, omiga, path_types{path_index});
        line_scans{path_index} = build_ROI_line_scan( ...
            scan_cycles{path_index});
    end

    for repeat_index = 1:repeat_count
        targets = target_direction{repeat_index};

        [detect_R_RA(:, repeat_index), ...
            detect_L_RA(:, repeat_index), ...
            detect_P_RA(:, repeat_index)] = static_coverage_all_beams( ...
            targets, scan_cycles{1}, line_scans{1}, ...
            fasan_D, fasan_d);

        [detect_R_SP(:, repeat_index), ...
            detect_L_SP(:, repeat_index), ...
            detect_P_SP(:, repeat_index)] = static_coverage_all_beams( ...
            targets, scan_cycles{2}, line_scans{2}, ...
            fasan_D, fasan_d);
    end

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

    % MC_ROI_NEW defines a point beam by fasan_D alone; fasan_d is not a
    % point-beam parameter and is therefore omitted from these file names.
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

function position = sample_cone_end_cap(R, jiaodu, N)
%SAMPLE_CONE_END_CAP Uniform targets on the ROI outer spherical end cap.
% MC_ROI_NEW bounds the ROI by range <= R and polar angle <= jiaodu.
% Consequently its outer "bottom" is the spherical cap at range R, not the
% conical side.  Uniform surface sampling is uniform in azimuth and cos(theta).

alpha = deg2rad(jiaodu);
azimuth = 2 * pi * rand(N, 1);
cos_polar = cos(alpha) + (1 - cos(alpha)) * rand(N, 1);
sin_polar = sqrt(max(1 - cos_polar.^2, 0));

position = R .* [ ...
    sin_polar .* cos(azimuth), ...
    sin_polar .* sin(azimuth), ...
    cos_polar];
end

function [detect_ring, detect_line, detect_point] = ...
    static_coverage_all_beams(targets, beam_cycle, line_scan, ...
    fasan_D, fasan_d)
%STATIC_COVERAGE_ALL_BEAMS Geometric coverage over one complete scan.

target_block_size = 128;
beam_block_size = 8192;
N = size(targets, 1);

block_start = 1:target_block_size:N;
target_blocks = cell(numel(block_start), 1);
ring_blocks = cell(numel(block_start), 1);
line_blocks = cell(numel(block_start), 1);
point_blocks = cell(numel(block_start), 1);

for block_index = 1:numel(block_start)
    first_target = block_start(block_index);
    last_target = min(first_target + target_block_size - 1, N);
    target_blocks{block_index} = targets(first_target:last_target, :);
end

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

    parfor block_index = 1:numel(block_start)
        [ring_blocks{block_index}, line_blocks{block_index}, ...
            point_blocks{block_index}] = test_target_block( ...
            target_blocks{block_index}, beam_cycle, line_scan, ...
            fasan_D, fasan_d, beam_block_size);
    end
else
    for block_index = 1:numel(block_start)
        [ring_blocks{block_index}, line_blocks{block_index}, ...
            point_blocks{block_index}] = test_target_block( ...
            target_blocks{block_index}, beam_cycle, line_scan, ...
            fasan_D, fasan_d, beam_block_size);
    end
end

detect_ring = double(vertcat(ring_blocks{:}));
detect_line = double(vertcat(line_blocks{:}));
detect_point = double(vertcat(point_blocks{:}));
end

function [detected_ring, detected_line, detected_point] = ...
    test_target_block(targets, beam_cycle, line_scan, ...
    fasan_D, fasan_d, beam_block_size)
%TEST_TARGET_BLOCK Stop testing a beam type after each target's first hit.

target_count = size(targets, 1);
detected_ring = false(target_count, 1);
detected_line = false(target_count, 1);
detected_point = false(target_count, 1);

% Definitions inherited from the NEW dynamic model:
%   ring  - centre radius fasan_D/2, exp(-2) radial width fasan_d;
%   line  - long full angle fasan_D, exp(-2) short width fasan_d;
%   point - circular exp(-2) full angle fasan_D.
ring_inner = max((fasan_D - fasan_d) / 2, 0);
ring_outer = min((fasan_D + fasan_d) / 2, pi);
ring_cos_upper = cos(ring_inner);
ring_cos_lower = cos(ring_outer);
point_cos_limit = cos(fasan_D / 2);
line_long_tangent = tan(fasan_D / 2);
line_short_tangent = tan(fasan_d / 2);

beam_count = size(beam_cycle, 1);
for first_beam = 1:beam_block_size:beam_count
    if all(detected_ring) && all(detected_line) && all(detected_point)
        break;
    end

    last_beam = min(first_beam + beam_block_size - 1, beam_count);
    beam_index = first_beam:last_beam;
    beam_block = beam_cycle(beam_index, :);

    % Ring and point coverage depend only on the angular separation from
    % the instantaneous optical axis.  Dot-product bounds avoid acos and
    % are exactly equivalent over [0, pi].
    angular_active = find(~detected_ring | ~detected_point);
    if ~isempty(angular_active)
        direction_cosine = targets(angular_active, :) * beam_block.';

        need_ring = ~detected_ring(angular_active);
        if any(need_ring)
            ring_hit = any( ...
                direction_cosine(need_ring, :) >= ring_cos_lower & ...
                direction_cosine(need_ring, :) <= ring_cos_upper, 2);
            ring_target = angular_active(need_ring);
            detected_ring(ring_target(ring_hit)) = true;
        end

        need_point = ~detected_point(angular_active);
        if any(need_point)
            point_hit = any( ...
                direction_cosine(need_point, :) >= point_cos_limit, 2);
            point_target = angular_active(need_point);
            detected_point(point_target(point_hit)) = true;
        end
    end

    % The line footprint uses the same spherical tangent frame as
    % MC_line_RA_snr_NEW: the short axis follows the scan tangent and the
    % long axis is perpendicular to it in the beam tangent plane.
    line_active = find(~detected_line);
    if ~isempty(line_active)
        target_block = targets(line_active, :);
        along_beam = target_block * beam_block.';
        long_component = target_block * ...
            line_scan.LongAxes(beam_index, :).';
        short_component = target_block * ...
            line_scan.ShortAxes(beam_index, :).';

        inside_line = along_beam > 0 & ...
            abs(long_component) <= along_beam .* line_long_tangent & ...
            abs(short_component) <= along_beam .* line_short_tangent;
        line_hit = any(inside_line, 2);
        detected_line(line_active(line_hit)) = true;
    end
end
end

function validate_parameters(N, f, R, jiaodu, fasan_D_list, ...
    fasan_d, omi_type, repeat_count)
%VALIDATE_PARAMETERS Fail early for invalid physical or batch parameters.

if ~isscalar(N) || N <= 0 || N ~= floor(N)
    error('static_MC_ROI_NEW:InvalidTargetCount', ...
        'N must be a positive integer.');
end
if ~isscalar(f) || f <= 0
    error('static_MC_ROI_NEW:InvalidPulseRate', ...
        'f must be positive.');
end
if ~isscalar(R) || R <= 0
    error('static_MC_ROI_NEW:InvalidRange', ...
        'R must be positive.');
end
if ~isscalar(jiaodu) || jiaodu <= 0 || jiaodu >= 90
    error('static_MC_ROI_NEW:InvalidAngle', ...
        'JIAODU must be between 0 and 90 degrees.');
end
if isempty(fasan_D_list) || any(fasan_D_list <= 0 | fasan_D_list >= pi)
    error('static_MC_ROI_NEW:InvalidBeamWidth', ...
        'Every fasan_D value must be between 0 and pi radians.');
end
if ~isscalar(fasan_d) || fasan_d <= 0 || fasan_d >= pi
    error('static_MC_ROI_NEW:InvalidGaussianWidth', ...
        'fasan_d must be between 0 and pi radians.');
end
if ~any(strcmpi(omi_type, {'fast', 'slow'}))
    error('static_MC_ROI_NEW:UnknownOmiType', ...
        'omi_type must be fast or slow.');
end
if repeat_count ~= 2
    error('static_MC_ROI_NEW:InvalidRepeatCount', ...
        'repeat_count must be 2 to preserve the legacy result structure.');
end
end

function text = format_mrad(angle)
%FORMAT_MRAD Match the numeric filename convention of static_MC_ROI.m.

text = num2str(angle * 1000, '%g');
end

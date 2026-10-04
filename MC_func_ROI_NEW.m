function [detect, first_time, first_encounter_time, pulse_score] = ...
    MC_func_ROI_NEW(R, N, f, jiaodu, beam_type, ...
    fasan_D, fasan_d, v_min, v_max, omi_type, path_type, window_mode, mc_block_size)
%MC_FUNC_ROI_NEW Cone ROI Monte Carlo dispatcher for ROI beams.

isRing = strcmpi(beam_type, 'ring');
isLine = strcmpi(beam_type, 'line');
isPoint = strcmpi(beam_type, 'point');
if ~(isRing || isLine || isPoint)
    error('MC_func_ROI_NEW:UnsupportedBeamType', ...
        'Unknown beam_type: %s. Use ring, line, or point.', beam_type);
end

if strcmpi(omi_type, 'fast')
    omiga = 2 * pi;
elseif strcmpi(omi_type, 'slow')
    omiga = pi / 2;
else
    error('MC_func_ROI_NEW:UnknownOmiType', ...
        'Unknown omi_type: %s. Use fast or slow.', omi_type);
end

if jiaodu <= 0 || jiaodu >= 90
    error('MC_func_ROI_NEW:InvalidAngle', ...
        'jiaodu must be between 0 and 90 degrees.');
end
if fasan_D <= 0 || fasan_D >= pi
    error('MC_func_ROI_NEW:InvalidBeamWidth', ...
        'fasan_D must be between 0 and pi radians.');
end
if f <= 0
    error('MC_func_ROI_NEW:InvalidPulseRate', 'f must be positive.');
end
if nargin < 12 || isempty(window_mode)
    window_mode = 'variable';
end
window_mode = validatestring(window_mode, {'variable', 'fixed'}, ...
    mfilename, 'window_mode');
if nargin < 13 || isempty(mc_block_size)
    mc_block_size = 100000;
end
validateattributes(mc_block_size, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'positive'}, ...
    mfilename, 'mc_block_size');


[UAV_P, UAV_V] = init_UAV_ROI_NEW(R, jiaodu, v_min, v_max, N);

beam_cycle = build_ROI_scan_cycle(jiaodu, fasan_D, f, omiga, path_type);
beam_cycle_len = size(beam_cycle, 1);
beam_idx_rand = randi([1, beam_cycle_len], N, 1);
if isLine
    line_scan = build_ROI_line_scan(beam_cycle);
else
    line_scan = [];
end

detect = zeros(N, 1);
first_detect = nan(N, 1);
first_encounter = nan(N, 1);
pulse_score = repmat(empty_detection_window_photons(), N, 1);

UAV_P_x = UAV_P(:, 1); UAV_P_y = UAV_P(:, 2); UAV_P_z = UAV_P(:, 3);
UAV_V_x = UAV_V(:, 1); UAV_V_y = UAV_V(:, 2); UAV_V_z = UAV_V(:, 3);

useParallel = license('test', 'Distrib_Computing_Toolbox') && ...
    ~isempty(ver('parallel'));
if useParallel
    p = gcp('nocreate');
    if isempty(p)
        try
            parpool('threads');
        catch
            parpool('local');
        end
    end
end

tic
if useParallel
    parfor i = 1:N
        UP = [UAV_P_x(i), UAV_P_y(i), UAV_P_z(i)];
        UV = [UAV_V_x(i), UAV_V_y(i), UAV_V_z(i)];
        init_beam_idx = beam_idx_rand(i);
        if strcmpi(path_type, 'raster')
            if isRing
                [detect(i), first_detect(i), first_encounter(i), ...
                    pulse_score(i)] = MC_ring_RA_snr_NEW( ...
                    fasan_D, fasan_d, beam_cycle, f, omiga, UP, UV, ...
                    init_beam_idx, R, jiaodu, window_mode, mc_block_size);
            elseif isLine
                [detect(i), first_detect(i), first_encounter(i), ...
                    pulse_score(i)] = MC_line_RA_snr_NEW( ...
                    fasan_D, fasan_d, line_scan, f, omiga, UP, UV, ...
                    init_beam_idx, R, jiaodu, window_mode, mc_block_size);
            else
                [detect(i), first_detect(i), first_encounter(i), ...
                    pulse_score(i)] = MC_point_RA_snr_NEW( ...
                    fasan_D, fasan_d, beam_cycle, f, omiga, UP, UV, ...
                    init_beam_idx, R, jiaodu, window_mode, mc_block_size);
            end
        else
            if isRing
                [detect(i), first_detect(i), first_encounter(i), ...
                    pulse_score(i)] = MC_ring_SP_snr_NEW( ...
                    fasan_D, fasan_d, beam_cycle, f, omiga, UP, UV, ...
                    init_beam_idx, R, jiaodu, window_mode, mc_block_size);
            elseif isLine
                [detect(i), first_detect(i), first_encounter(i), ...
                    pulse_score(i)] = MC_line_SP_snr_NEW( ...
                    fasan_D, fasan_d, line_scan, f, omiga, UP, UV, ...
                    init_beam_idx, R, jiaodu, window_mode, mc_block_size);
            else
                [detect(i), first_detect(i), first_encounter(i), ...
                    pulse_score(i)] = MC_point_SP_snr_NEW( ...
                    fasan_D, fasan_d, beam_cycle, f, omiga, UP, UV, ...
                    init_beam_idx, R, jiaodu, window_mode, mc_block_size);
            end
        end
    end
else
    for i = 1:N
        UP = [UAV_P_x(i), UAV_P_y(i), UAV_P_z(i)];
        UV = [UAV_V_x(i), UAV_V_y(i), UAV_V_z(i)];
        init_beam_idx = beam_idx_rand(i);
        if strcmpi(path_type, 'raster')
            if isRing
                [detect(i), first_detect(i), first_encounter(i), ...
                    pulse_score(i)] = MC_ring_RA_snr_NEW( ...
                    fasan_D, fasan_d, beam_cycle, f, omiga, UP, UV, ...
                    init_beam_idx, R, jiaodu, window_mode, mc_block_size);
            elseif isLine
                [detect(i), first_detect(i), first_encounter(i), ...
                    pulse_score(i)] = MC_line_RA_snr_NEW( ...
                    fasan_D, fasan_d, line_scan, f, omiga, UP, UV, ...
                    init_beam_idx, R, jiaodu, window_mode, mc_block_size);
            else
                [detect(i), first_detect(i), first_encounter(i), ...
                    pulse_score(i)] = MC_point_RA_snr_NEW( ...
                    fasan_D, fasan_d, beam_cycle, f, omiga, UP, UV, ...
                    init_beam_idx, R, jiaodu, window_mode, mc_block_size);
            end
        else
            if isRing
                [detect(i), first_detect(i), first_encounter(i), ...
                    pulse_score(i)] = MC_ring_SP_snr_NEW( ...
                    fasan_D, fasan_d, beam_cycle, f, omiga, UP, UV, ...
                    init_beam_idx, R, jiaodu, window_mode, mc_block_size);
            elseif isLine
                [detect(i), first_detect(i), first_encounter(i), ...
                    pulse_score(i)] = MC_line_SP_snr_NEW( ...
                    fasan_D, fasan_d, line_scan, f, omiga, UP, UV, ...
                    init_beam_idx, R, jiaodu, window_mode, mc_block_size);
            else
                [detect(i), first_detect(i), first_encounter(i), ...
                    pulse_score(i)] = MC_point_SP_snr_NEW( ...
                    fasan_D, fasan_d, beam_cycle, f, omiga, UP, UV, ...
                    init_beam_idx, R, jiaodu, window_mode, mc_block_size);
            end
        end
    end
end
toc

SUCCESS = sum(detect) / length(detect);
avetime = mean(first_detect, 'omitnan') / f;
disp([omi_type, '_', path_type, '_', num2str(fasan_D), '_', ...
    num2str(fasan_d), '_', num2str(SUCCESS), '_', num2str(avetime)]);

first_time = first_detect ./ f;
first_encounter_time = first_encounter ./ f;
end

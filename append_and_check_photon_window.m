function [detect, first_detect, window_photons, ...
    carry_sig, carry_backscatter, carry_step] = ...
    append_and_check_photon_window( ...
    carry_sig, carry_backscatter, carry_step, ...
    sig_arr, backscatter_arr, step_arr, background_photons, ...
    dark_count_photons, window_pulses, threshold, window_mode)
%APPEND_AND_CHECK_PHOTON_WINDOW Apply a selectable causal-window criterion.
% VARIABLE keeps the legacy shortest-successful-window behavior. FIXED uses
% all available startup pulses, then exactly WINDOW_PULSES for every window.

detect = 0;
first_detect = nan;
window_photons = empty_detection_window_photons();

if nargin < 11 || isempty(window_mode)
    window_mode = 'variable';
end
window_mode = validatestring(window_mode, {'variable', 'fixed'}, ...
    mfilename, 'window_mode');

sig_arr = sig_arr(:);
backscatter_arr = backscatter_arr(:);
step_arr = step_arr(:);
if isempty(sig_arr)
    return;
end

sig_arr(~isfinite(sig_arr)) = 0;
backscatter_arr(~isfinite(backscatter_arr)) = 0;

n_carry = numel(carry_sig);
if strcmp(window_mode, 'fixed')
    % A zero-signal pulse is never a candidate tail in the legacy rule.
    % Skip empty blocks without creating full-block score or tail arrays.
    new_signal_pos = find(sig_arr > 0);
    if isempty(new_signal_pos)
        [carry_sig, carry_backscatter, carry_step] = keep_recent_photons( ...
            carry_sig, carry_backscatter, carry_step, ...
            sig_arr, backscatter_arr, step_arr, window_pulses);
        return;
    end

    all_sig = [carry_sig; sig_arr];
    all_backscatter = [carry_backscatter; backscatter_arr];
    all_step = [carry_step; step_arr];
    all_pulse_var = all_sig + all_backscatter + ...
        background_photons + dark_count_photons;
    score_valid = all_sig > 0 & all_pulse_var > 0 & ...
        isfinite(all_pulse_var);
    score_pos = find(score_valid);
    score_value = all_sig(score_pos).^2 ./ all_pulse_var(score_pos);
    cumulative_score = [0; cumsum(score_value)];

    % Check only newly appended positive-signal tails. Both score-event edges
    % move monotonically, giving a streaming sum over sparse score events.
    candidate_tail_pos = n_carry + new_signal_pos;
    left_event = 1;
    right_event = 0;
    n_score_events = numel(score_pos);
    successful_tail = [];
    successful_start = [];
    for candidate_idx = 1:numel(candidate_tail_pos)
        tail_pos = candidate_tail_pos(candidate_idx);
        start_pos = max(1, tail_pos - window_pulses + 1);
        while right_event < n_score_events && ...
                score_pos(right_event + 1) <= tail_pos
            right_event = right_event + 1;
        end
        while left_event <= right_event && score_pos(left_event) < start_pos
            left_event = left_event + 1;
        end
        if left_event <= right_event
            win_score = cumulative_score(right_event + 1) - ...
                cumulative_score(left_event);
        else
            win_score = 0;
        end
        if all(win_score > 0 & isfinite(win_score) & ...
                sqrt(win_score) >= threshold)
            successful_tail = tail_pos;
            successful_start = start_pos;
            break;
        end
    end

    if ~isempty(successful_tail)
        selected = successful_start:successful_tail;
        detect = 1;
        first_detect = all_step(successful_tail);
        window_photons.signal_photons = all_sig(selected);
        window_photons.backscatter_photons = ...
            all_backscatter(selected);
        window_photons.background_photons = background_photons;
        window_photons.dark_count_photons = dark_count_photons;
        window_photons.pulse_steps = uint32(all_step(selected));
        return;
    end
else
    % Variable mode follows the legacy shortest-successful-window rule. An
    % empty signal block cannot detect, so retain only the carry needed by a
    % future block and avoid allocating score arrays for the full block.
    if ~any(sig_arr > 0)
        [carry_sig, carry_backscatter, carry_step] = keep_recent_photons( ...
            carry_sig, carry_backscatter, carry_step, ...
            sig_arr, backscatter_arr, step_arr, window_pulses);
        return;
    end

    all_sig = [carry_sig; sig_arr];
    all_backscatter = [carry_backscatter; backscatter_arr];
    all_step = [carry_step; step_arr];
    n_all = numel(all_sig);
    all_noise = all_backscatter + background_photons + dark_count_photons;
    pulse_var = all_sig + all_noise;
    per_pulse_score = zeros(size(all_sig));
    score_valid = all_sig > 0 & pulse_var > 0 & isfinite(pulse_var);
    per_pulse_score(score_valid) = all_sig(score_valid).^2 ./ ...
        pulse_var(score_valid);
    cumulative_score = [0; cumsum(per_pulse_score)];

    signal_tail_pos = n_carry + find(sig_arr > 0);
    if ~isempty(signal_tail_pos)
        tail_pos = signal_tail_pos(:);
        max_len = min(window_pulses, n_all);
        len_vec = 1:max_len;
        start_pos = tail_pos - len_vec + 1;
        valid_len = start_pos >= 1;
        start_pos(~valid_len) = 1;

        tail_score = cumulative_score(tail_pos + 1);
        start_score = reshape(cumulative_score(start_pos(:)), ...
            size(start_pos));
        win_score = tail_score - start_score;
        win_snr = sqrt(win_score);
        win_snr(~valid_len | win_score <= 0 | ~isfinite(win_snr)) = -inf;

        hit_matrix = win_snr >= threshold;
        first_successful_tail = find(any(hit_matrix, 2), 1, 'first');
        if ~isempty(first_successful_tail)
            successful_len = find( ...
                hit_matrix(first_successful_tail, :), 1, 'first');
            successful_start = start_pos( ...
                first_successful_tail, successful_len);
            successful_tail = tail_pos(first_successful_tail);
            selected = successful_start:successful_tail;

            detect = 1;
            first_detect = all_step(successful_tail);
            window_photons.signal_photons = all_sig(selected);
            window_photons.backscatter_photons = ...
                all_backscatter(selected);
            window_photons.background_photons = background_photons;
            window_photons.dark_count_photons = dark_count_photons;
            window_photons.pulse_steps = uint32(all_step(selected));
            return;
        end
    end
end

[carry_sig, carry_backscatter, carry_step] = keep_recent_photons( ...
    carry_sig, carry_backscatter, carry_step, ...
    sig_arr, backscatter_arr, step_arr, window_pulses);
end

function [carry_sig, carry_backscatter, carry_step] = keep_recent_photons( ...
    carry_sig, carry_backscatter, carry_step, ...
    sig_arr, backscatter_arr, step_arr, window_pulses)
%KEEP_RECENT_PHOTONS Retain raw values for the next causal window.

n_new = numel(sig_arr);
n_keep = min(window_pulses - 1, numel(carry_sig) + n_new);
if n_keep <= 0
    carry_sig = zeros(0, 1);
    carry_backscatter = zeros(0, 1);
    carry_step = zeros(0, 1);
elseif n_new >= n_keep
    first_new = n_new - n_keep + 1;
    carry_sig = sig_arr(first_new:end);
    carry_backscatter = backscatter_arr(first_new:end);
    carry_step = step_arr(first_new:end);
else
    n_old = n_keep - n_new;
    first_old = numel(carry_sig) - n_old + 1;
    carry_sig = [carry_sig(first_old:end); sig_arr];
    carry_backscatter = [ ...
        carry_backscatter(first_old:end); backscatter_arr];
    carry_step = [carry_step(first_old:end); step_arr];
end
end

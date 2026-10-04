function [post_photons, is_complete] = append_post_detection_photons( ...
    post_photons, sig_arr, backscatter_arr, step_arr, ...
    first_post_step, last_post_step, background_photons, dark_count_photons)
%APPEND_POST_DETECTION_PHOTONS Append a continuous post-detection segment.

sig_arr = sig_arr(:);
backscatter_arr = backscatter_arr(:);
step_arr = step_arr(:);

if isempty(post_photons.pulse_steps)
    next_step = first_post_step;
else
    next_step = double(post_photons.pulse_steps(end)) + 1;
end

selected = step_arr >= next_step & step_arr <= last_post_step;
if any(selected)
    post_photons.signal_photons = [ ...
        post_photons.signal_photons; sig_arr(selected)];
    post_photons.backscatter_photons = [ ...
        post_photons.backscatter_photons; backscatter_arr(selected)];
    post_photons.background_photons = background_photons;
    post_photons.dark_count_photons = dark_count_photons;
    post_photons.pulse_steps = [ ...
        post_photons.pulse_steps; uint32(step_arr(selected))];
end

is_complete = ~isempty(post_photons.pulse_steps) && ...
    double(post_photons.pulse_steps(end)) >= last_post_step;
end

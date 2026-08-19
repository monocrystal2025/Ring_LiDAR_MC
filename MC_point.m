function [detect, first_detect, first_encounter, effective_pulses] = MC_point(fasan_D, fasan_d, path1, beam_vec, path2, f, omiga, UP, UV, init_beam_idx, R)
%MC_POINT  Compatibility wrapper for the Gaussian point-beam model.
%
% fasan_D is the full divergence angle at the exp(-2) intensity contour.
% fasan_d is ignored for point beams.

[detect, first_detect, first_encounter, effective_pulses] = ...
    MC_point_snr(fasan_D, fasan_d, ...
    path1, beam_vec, path2, f, omiga, UP, UV, init_beam_idx, R);
end

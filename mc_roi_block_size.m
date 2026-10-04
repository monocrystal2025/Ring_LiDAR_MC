function blockSize = mc_roi_block_size(fasanD, beamType)
%MC_ROI_BLOCK_SIZE Select a measured pulse block size for ROI simulations.
%   BLOCKSIZE = MC_ROI_BLOCK_SIZE(FASAND, BEAMTYPE) maps the full
%   divergence angle FASAND (radians) and BEAMTYPE to a pulse-vector block
%   size. The policy is tuned for the cone-ROI kernels with N = 50000 and a
%   six-worker thread pool. It changes performance only; window state is
%   carried across block boundaries by append_and_check_photon_window.

validateattributes(fasanD, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, mfilename, 'fasanD');
beamType = validatestring(beamType, {'ring', 'line', 'point'}, ...
    mfilename, 'beamType');

% Round only for region selection so 5 mrad colon-generated values land on
% the intended inclusive boundary despite floating-point roundoff.
fasanDMrad = round(fasanD * 1e3, 9);

if fasanDMrad <= 10
    regionIndex = 1;
elseif fasanDMrad <= 50
    regionIndex = 2;
elseif fasanDMrad <= 150
    regionIndex = 3;
elseif fasanDMrad <= 250
    regionIndex = 4;
else
    regionIndex = 5;
end

switch beamType
    case 'ring'
        regionBlockSizes = [30000, 2000, 2000, 5000, 5000];
    case 'line'
        regionBlockSizes = [30000, 2000, 2000, 5000, 5000];
    case 'point'
        regionBlockSizes = [20000, 2000, 10000, 10000, 5000];
end

blockSize = regionBlockSizes(regionIndex);
end

function SNR = MC_single_pulse_snr(beam_type, fasan_D, fasan_d, z)
%MC_SINGLE_PULSE_SNR  Single-pulse SNR used by MC hit confirmation.
%
%   SNR = MC_single_pulse_snr(beam_type, fasan_D, fasan_d, z)
%
% Inputs
%   beam_type : 'point', 'line', or 'ring'
%   fasan_D   : maximum divergence angle wD, rad
%               spot diameter / line length / annular outer diameter
%   fasan_d   : beam-width divergence wd, rad; ignored for point
%   z         : target range at the current hit event, m. Vector allowed.
%
% The physical parameters are kept identical to lidar_single_pulse_snr_demo.m.
% The MC contact criterion is still handled outside this function as a point-hit
% criterion. Once a hit occurs, this function only evaluates the flat-top
% energy-density-based SNR at the current range.

z = z(:);
SNR = nan(size(z));

%% Physical constants
h = 6.62607015e-34;               % J*s
c = 299792458;                    % m/s

%% System parameters, identical to lidar_single_pulse_snr_demo.m
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

%% Target parameters, identical to lidar_single_pulse_snr_demo.m
target_w = 0.297;                 % m
A_target = target_w * target_w;   % m^2
rho = 0.5;

%% Common quantities
E_photon = h * c / lambda;
N_tx = E_pulse / E_photon;
T2_target = exp(-2 * alpha .* z);
N_dark = dark_count_rate * tau_gate;

%% Beam-shape-dependent illuminated area and matched-FOV solid angle
switch lower(char(beam_type))
    case 'point'
        theta_D = fasan_D;
        spot_D = theta_D .* z;
        beam_area = pi .* (spot_D ./ 2).^2;
        Omega_fov = pi .* (theta_D ./ 2).^2;
        valid = beam_area > 0;

    case 'line'
        theta_D = fasan_D;
        theta_w = fasan_d;
        spot_D = theta_D .* z;
        spot_w = theta_w .* z;
        beam_area = spot_D .* spot_w;
        Omega_fov = theta_D .* theta_w;
        valid = beam_area > 0;

    case 'ring'
        theta_D = fasan_D;
        theta_w = fasan_d;
        outer_D = theta_D .* z;
        ring_w = theta_w .* z;
        inner_D = outer_D - 2 .* ring_w;

        % Allow inner_D = 0, which degenerates to a filled disk. Reject only
        % nonphysical negative inner diameter.
        valid = inner_D >= 0;
        Rout = outer_D ./ 2;
        Rin = max(inner_D, 0) ./ 2;
        beam_area = pi .* (Rout.^2 - Rin.^2);
        theta_inner = max(theta_D - 2 .* theta_w, 0);
        Omega_fov = pi/4 .* (theta_D.^2 - theta_inner.^2);
        valid = valid & beam_area > 0;

    otherwise
        error('Unknown beam_type: %s. Use point, line, or ring.', beam_type);
end

if ~any(valid)
    return;
end

%% Noise and signal photons
N_bg = calc_background_photons( ...
    L_sky_nm, delta_lambda_nm, A_rx, Omega_fov, tau_gate, eta_sys, E_photon);

N_bs = calc_backscatter_photons_vec( ...
    N_tx, eta_sys, A_rx, beta, alpha, z, tau_gate, c);

N_sig = N_tx .* (A_target ./ beam_area) .* rho .* ...
        (A_rx ./ (pi .* z.^2)) .* T2_target .* eta_sys;

N_noise = N_bs + N_bg + N_dark;
SNR(valid) = N_sig(valid) ./ sqrt(N_noise(valid));
end

function N_bs = calc_backscatter_photons_vec(N_tx, eta_sys, A_rx, beta, alpha, z, tau_gate, c)
% Vectorized version of the backscatter calculation in the SNR script.
% Integration is performed over the range gate centered at each z.

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

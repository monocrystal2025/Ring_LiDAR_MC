function [SNR, N_sig, N_noise] = MC_single_pulse_snr(beam_type, fasan_D, fasan_d, z, illumination_factor, ring_profile)
%MC_SINGLE_PULSE_SNR  Single-pulse SNR used by MC hit confirmation.
%
%   SNR = MC_single_pulse_snr(beam_type, fasan_D, fasan_d, z)
%   [SNR, N_sig, N_noise] = MC_single_pulse_snr(...)
%
% Inputs
%   beam_type : 'point', 'line', or 'ring'
%   fasan_D   : maximum divergence angle wD, rad
%               spot diameter / line length / annular outer diameter
%   fasan_d   : beam-width divergence wd, rad; ignored for point
%   z         : target range at the current hit event, m. Vector allowed.
%   illumination_factor : normalized mean target illumination, where 1 means
%                         the whole target sees the peak beam intensity.
%   ring_profile        : 'flat' keeps the legacy annular top-hat model;
%                         'gaussian' uses fasan_d as the radial 1/e^2 width.
%                         Its receiver FOV is a hard annular mask centered
%                         at fasan_D/2 with total width 6*fasan_d.
%   For point beams, fasan_D is the full angle at the exp(-2) Gaussian
%   intensity contour, and the receiver FOV full angle is 3*fasan_D.
%
% The physical parameters are kept identical to lidar_single_pulse_snr_demo.m.
% Single-pulse detectability includes signal shot noise:
%   SNR = N_sig / sqrt(N_sig + N_noise).

if nargin < 5 || isempty(illumination_factor)
    illumination_factor = 1;
end
if nargin < 6 || isempty(ring_profile)
    ring_profile = 'flat';
end

z = z(:);
if isscalar(illumination_factor)
    illumination_factor = repmat(illumination_factor, size(z));
else
    illumination_factor = illumination_factor(:);
end
if numel(illumination_factor) ~= numel(z)
    error('MC_single_pulse_snr:SizeMismatch', ...
        'illumination_factor must be scalar or have the same number of elements as z.');
end
illumination_factor = max(illumination_factor, 0);
SNR = nan(size(z));
N_sig = nan(size(z));
N_noise = nan(size(z));

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
        gaussian_half_angle = fasan_D / 2;
        fov_half_angle = min(3 * gaussian_half_angle, pi);
        Omega_beam = point_gaussian_effective_solid_angle( ...
            gaussian_half_angle);
        Omega_fov = 2 * pi * (1 - cos(fov_half_angle));
        beam_area = Omega_beam .* z.^2;
        valid = gaussian_half_angle > 0 & Omega_fov > 0 & ...
            Omega_beam > 0 & beam_area > 0;

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
        switch lower(char(ring_profile))
            case 'flat'
                outer_D = theta_D .* z;
                ring_w = theta_w .* z;
                inner_D = outer_D - 2 .* ring_w;

                % Allow inner_D = 0, which degenerates to a filled disk.
                valid = inner_D >= 0;
                Rout = outer_D ./ 2;
                Rin = max(inner_D, 0) ./ 2;
                beam_area = pi .* (Rout.^2 - Rin.^2);
                theta_inner = max(theta_D - 2 .* theta_w, 0);
                Omega_fov = pi/4 .* (theta_D.^2 - theta_inner.^2);
                valid = valid & beam_area > 0;

            case 'gaussian'
                ring_radius = theta_D / 2;
                gaussian_width = theta_w;
                Omega_beam = annular_gaussian_effective_solid_angle( ...
                    ring_radius, gaussian_width);

                % Hard annular receiver mask: transmission is one inside
                % [theta_D/2-3*theta_w, theta_D/2+3*theta_w] and zero
                % everywhere else.  Use the exact spherical solid angle;
                % only two scalar cos evaluations are needed per call.
                fov_half_width = 3 * theta_w;
                fov_center = theta_D / 2;
                theta_fov_inner = max(fov_center - fov_half_width, 0);
                theta_fov_outer = min(fov_center + fov_half_width, pi);
                Omega_fov = 2 * pi * ...
                    (cos(theta_fov_inner) - cos(theta_fov_outer));

                % The transmitted beam remains radially Gaussian.  Keep
                % its effective area separate from the hard receiver FOV.
                beam_area = Omega_beam .* z.^2;
                valid = theta_D >= 0 & gaussian_width > 0 & ...
                    theta_fov_outer > theta_fov_inner & ...
                    Omega_fov > 0 & Omega_beam > 0 & beam_area > 0;

            otherwise
                error('Unknown ring_profile: %s. Use flat or gaussian.', ...
                    ring_profile);
        end

    otherwise
        error('Unknown beam_type: %s. Use point, line, or ring.', beam_type);
end

valid = valid & isfinite(z) & z > 0 & isfinite(illumination_factor);
if ~any(valid)
    return;
end

%% Noise and signal photons
N_bg = calc_background_photons( ...
    L_sky_nm, delta_lambda_nm, A_rx, Omega_fov, tau_gate, eta_sys, E_photon);

N_bs = calc_backscatter_photons_vec( ...
    N_tx, eta_sys, A_rx, beta, alpha, z, tau_gate, c);

N_sig_all = N_tx .* (A_target .* illumination_factor ./ beam_area) .* rho .* ...
        (A_rx ./ (pi .* z.^2)) .* T2_target .* eta_sys;

N_noise_all = N_bs + N_bg + N_dark;
N_sig(valid) = N_sig_all(valid);
N_noise(valid) = N_noise_all(valid);
SNR(valid) = N_sig(valid) ./ sqrt(N_sig(valid) + N_noise(valid));
end

function Omega_eff = point_gaussian_effective_solid_angle(gaussian_half_angle)
% Exact spherical solid angle of exp(-2*(theta/theta_w)^2).

if gaussian_half_angle <= 0
    Omega_eff = nan;
    return;
end

integrand = @(theta) exp(-2 .* (theta ./ gaussian_half_angle).^2) .* sin(theta);
Omega_eff = 2 * pi * integral(integrand, 0, pi, ...
    'RelTol', 1e-10, 'AbsTol', 1e-14);
end

function Omega_eff = annular_gaussian_effective_solid_angle(ring_radius, gaussian_width)
% Integral of exp(-2*((theta-ring_radius)/gaussian_width)^2) over angle.

if gaussian_width <= 0
    Omega_eff = nan;
    return;
end

u0 = -sqrt(2) * ring_radius / gaussian_width;
radial_integral = ring_radius * gaussian_width * sqrt(pi) / (2 * sqrt(2)) .* ...
    erfc(u0) + gaussian_width.^2 / 4 .* exp(-2 * (ring_radius ./ gaussian_width).^2);
Omega_eff = 2 * pi .* radial_integral;
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

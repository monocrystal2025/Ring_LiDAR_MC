function [detect,first_detect]=MC_ring_ROI_SP(fasan_D,fasan_d,beam_vec,path2,f,omiga,UP,UV,init_beam_idx,z,a)
%MC_RING_ROI_SP  ROI-SP annular-beam MC detection with SNR-based confirmation.
%
% The geometric contact condition is kept as the original ROI point-hit test:
% the target is simplified as a point in the ROI plane, and a hit is judged by
% whether the point falls within the annular footprint at the current pulse.
%
% New detection rule:
%   SNR >= 2.80          : one hit is accepted as successful detection.
%   1.64 <= SNR < 2.80  : the hit is stored as a candidate. Detection is
%                         accepted only if another hit with SNR >= 1.64 occurs
%                         before the LiDAR pointing direction moves farther
%                         than 2*wD from the candidate hit.
%   SNR < 1.64          : the hit is ignored.
%
% first_detect is the pulse index at which the warning is actually generated.
% For a weak two-hit confirmation, this is the second-hit index.

% omiga is retained for compatibility with the original interface.

if nargin < 11
    error('MC_ring_ROI_SP requires 11 input arguments.');
end

detect=0;
first_detect=nan;
N_KERNAL=100000;

SNR_WEAK = 1.64;
SNR_STRONG = 2.80;
CONFIRM_ANGLE = 2 * fasan_D;

D = 2*z*tan(fasan_D/2);
d = 2*z*tan(fasan_d/2);
R_ring_outer = D / 2;
R_ring_inner = max(R_ring_outer - d, 0);

pending_valid = false;
pending_beam = [nan, nan];

Beam_all=[beam_vec; path2; flipud(beam_vec)];
L_cycle = size(Beam_all,1);
idx0 = mod(init_beam_idx - 1, L_cycle) + 1;

%% Check the initial pulse, corresponding to first_detect = 0
B0 = Beam_all(idx0,:);
d2_0 = sum((UP - B0).^2, 2);
if (d2_0 <= R_ring_outer.^2) && (d2_0 >= R_ring_inner.^2)
    z0 = roi_slant_range(UP, z);
    snr0 = MC_single_pulse_snr('ring', fasan_D, fasan_d, z0);
    if snr0 >= SNR_STRONG
        detect=1;
        first_detect=0;
        return
    elseif snr0 >= SNR_WEAK
        pending_valid = true;
        pending_beam = B0;
    end
end

P_now=[ones(N_KERNAL,1)*UP(1),ones(N_KERNAL,1)*UP(2)];
dt=(0:1/f:N_KERNAL*1/f)';
dt(1,:)=[];
dx=UV(1)*dt; dy=UV(2)*dt;
dp=[dx,dy];
i=1;
P_now=P_now+dp;
B_now=cut_loop_simple_anyk(beam_vec,path2,init_beam_idx+1,N_KERNAL);
beam_idx=init_beam_idx+1;

while true
    norm_P_now=vecnorm(P_now,2,2);
    d2 = sum((P_now - B_now).^2, 2);
    judge=((d2<=R_ring_outer.^2)&(d2>=R_ring_inner.^2));
    judge(norm_P_now>a)=false;

    hit_idx = find(judge);
    if ~isempty(hit_idx)
        for kk = 1:numel(hit_idx)
            this_z = roi_slant_range(P_now(hit_idx(kk),:), z);
            this_snr = MC_single_pulse_snr('ring', fasan_D, fasan_d, this_z);
            this_beam = B_now(hit_idx(kk),:);
            this_step = (i-1)*N_KERNAL + hit_idx(kk);

            if this_snr >= SNR_STRONG
                detect=1;
                first_detect=this_step;
                return
            elseif this_snr >= SNR_WEAK
                if pending_valid
                    beam_sep = roi_beam_angle(pending_beam, this_beam, z);
                    if beam_sep <= CONFIRM_ANGLE
                        detect=1;
                        first_detect=this_step;
                        return
                    else
                        % The old weak hit is outside the allowed short-time
                        % angular window. Discard it and start a new candidate.
                        pending_beam = this_beam;
                    end
                else
                    pending_valid = true;
                    pending_beam = this_beam;
                end
            end
        end
    end

    norm_P_tmp = norm_P_now;
    norm_P_tmp(norm_P_tmp>a)=nan;
    i=i+1;
    if any(isnan(norm_P_tmp))
        break
    end

    P_now(:,1)=P_now(:,1)+N_KERNAL*UV(1)/f;
    P_now(:,2)=P_now(:,2)+N_KERNAL*UV(2)/f;
    beam_idx=beam_idx+N_KERNAL;
    B_now=cut_loop_simple_anyk(beam_vec,path2,beam_idx,N_KERNAL);
end
end

function seg = cut_loop_simple_anyk(a2,a3,k,m)
% k can be any absolute index in the infinite periodic sequence.
% seg is an m-by-2 segment starting from k.

if m<=0
    seg = zeros(0,size(a2,2));
    return;
end

P = [a2; a3; flipud(a2)];
L = size(P,1);
k0 = mod(k-1, L) + 1;
idx = mod((k0-1) + (0:m-1), L) + 1;
seg = P(idx, :);
end

function zr = roi_slant_range(P, z0)
% Convert ROI-plane coordinates to instantaneous LiDAR-target slant range.
% If a strict constant-range ROI model is desired, replace this function with:
% zr = z0 * ones(size(P,1),1);
if isvector(P)
    zr = sqrt(z0.^2 + sum(P(:).'.^2, 2));
else
    zr = sqrt(z0.^2 + sum(P.^2, 2));
end
zr = z0 * ones(size(P,1),1);
end

function a = roi_beam_angle(b1, b2, z0)
% Angular separation between two ROI beam pointings as seen by the LiDAR.
v1 = [b1(:).', z0];
v2 = [b2(:).', z0];
a = safe_acos(dot(v1,v2) ./ (norm(v1)*norm(v2)));
end

function y = safe_acos(x)
y = acos(max(min(x,1),-1));
end

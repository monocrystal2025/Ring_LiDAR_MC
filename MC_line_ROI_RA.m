function [detect,first_detect]=MC_line_ROI_RA(fasan_D,fasan_d,beam_vec,f,omiga,UP,UV,init_beam_idx,z,a)
detect=0;
first_detect=nan;
% UV=[15,15,15];
N_KERNAL=100000;
D = 2*z*tan(fasan_D/2);
d = 2*z*tan(fasan_d/2);
R_ring_outer = D / 2;
R_ring_inner = R_ring_outer - d;

UP_N=UP./norm(UP);

Beam_all=[beam_vec; flipud(beam_vec)]; 

%第一发能否成功

if (abs(UP(:,1)-Beam_all(init_beam_idx,1))<=D/2)&&(abs(UP(:,2)-Beam_all(init_beam_idx,2))<=d/2)
    detect=1;
    first_detect=0;
    return
end

P_now=[ones(N_KERNAL,1)*UP(1),ones(N_KERNAL,1)*UP(2)]; %初始化位置矩阵
dt=(0:1/f:N_KERNAL*1/f)';
dt(1,:)=[];
dx=UV(1)*dt;dy=UV(2)*dt;
dp=[dx,dy];
i=1;
P_now=P_now+dp;
B_now=cut_loop_simple_anyk(beam_vec,init_beam_idx+1,N_KERNAL);% 首次循环开始前beam向量更新
beam_idx=init_beam_idx+1;


while true
    

    norm_P_now=vecnorm(P_now,2,2); %每一时刻UAV位置与原点距离的平方
    d2 = sum((P_now - B_now).^2, 2);%UAV和当前光轴的距离平方，不开方更快
    ddx=abs(P_now(:,1)-B_now(:,1));
    ddy=abs(P_now(:,2)-B_now(:,2));
    judge=((ddx<D/2)&(ddy<d/2)); %判断标准
    judge(norm_P_now>a)=0; %超出的点强制置0
    if any(judge)
        detect=1;
        % 此处后面再写，首次预警时间
        first_detect=(i-1)*N_KERNAL+find(judge,1,'first');
        
        return
    end
    norm_P_now(norm_P_now>a)=nan;
    i=i+1;
    if anynan(norm_P_now) %算完了，退出
        break
    end
    P_now(:,1)=P_now(:,1)+N_KERNAL*UV(1)/f; %预更新下次循环的位置集合
    P_now(:,2)=P_now(:,2)+N_KERNAL*UV(2)/f;
    beam_idx=beam_idx+N_KERNAL;
    B_now=cut_loop_simple_anyk(beam_vec,beam_idx,N_KERNAL);% 预更新下次循环的BEAM集合
end






% max_step_count=2*R*f/norm(UV);




end


function seg = cut_loop_simple_anyk(a2,k,m)
% k: 可以是 P 中任意索引（1..L），也可以是无限序列中的绝对索引（>L）
% 返回：seg, 尺寸 m×3

    if m<=0
        seg = zeros(0,3);
        return;
    end

    % 一个完整周期
    P = [a2; flipud(a2);];
    L = size(P,1);

    % 关键：把任意 k 映射到周期内的起点（1..L）
    k0 = mod(k-1, L) + 1;

    % 从起点 k0 开始，连续取 m 个（超出 L 自动回绕）
    idx = mod((k0-1) + (0:m-1), L) + 1;

    seg = P(idx, :);
end


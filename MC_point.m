function [detect,first_detect]=MC_point(fasan_D,fasan_d,path1,beam_vec,path2,f,omiga,UP,UV,init_beam_idx,R)
detect=0;
first_detect=nan;
% UV=[15,15,15];
N_KERNAL=100000;


UP_N=UP./norm(UP);

Beam_all=[beam_vec; path2; flipud(beam_vec); path1]; %第一发能否成功
if (acos(dot(UP_N,Beam_all(init_beam_idx,:)))<=fasan_D/2)
    detect=1;
    first_detect=0;
    return
end

P_now=[ones(N_KERNAL,1)*UP(1),ones(N_KERNAL,1)*UP(2),ones(N_KERNAL,1)*UP(3)]; %初始化位置矩阵
dt=(0:1/f:N_KERNAL*1/f)';
dt(1,:)=[];
dx=UV(1)*dt;dy=UV(2)*dt;dz=UV(3)*dt;
dp=[dx,dy,dz];
i=1;
P_now=P_now+dp;
B_now=cut_loop_simple_anyk(path1,beam_vec,path2,init_beam_idx+1,N_KERNAL);% 首次循环开始前beam向量更新
beam_idx=init_beam_idx+1;


while true
    

    norm_P_now=vecnorm(P_now,2,2); %每一时刻UAV位置的单位向量

    dtheta=acos(dot(P_now,B_now,2)./norm_P_now); %UAV位置和BEAM指向的夹角
    judge=(dtheta<=fasan_D/2)); %判断标准
    judge(norm_P_now>R)=0; %超出的点强制置0
    if any(judge)
        detect=1;
        % 此处后面再写，首次预警时间
        first_detect=(i-1)*N_KERNAL+find(judge,1,'first');
        
        return
    end
    norm_P_now(norm_P_now>R)=nan;
    i=i+1;
    if anynan(norm_P_now) || any(P_now(:,3)<0) %算完了，退出
        break
    end
    P_now(:,1)=P_now(:,1)+N_KERNAL*UV(1)/f; %预更新下次循环的位置集合
    P_now(:,2)=P_now(:,2)+N_KERNAL*UV(2)/f;
    P_now(:,3)=P_now(:,3)+N_KERNAL*UV(3)/f;
    beam_idx=beam_idx+N_KERNAL;
    B_now=cut_loop_simple_anyk(path1,beam_vec,path2,beam_idx,N_KERNAL);% 预更新下次循环的BEAM集合
end






% max_step_count=2*R*f/norm(UV);




end


function seg = cut_loop_simple_anyk(a1,a2,a3,k,m)
% k: 可以是 P 中任意索引（1..L），也可以是无限序列中的绝对索引（>L）
% 返回：seg, 尺寸 m×3

    if m<=0
        seg = zeros(0,3);
        return;
    end

    % 一个完整周期
    P = [a2; a3; flipud(a2); a1];
    L = size(P,1);

    % 关键：把任意 k 映射到周期内的起点（1..L）
    k0 = mod(k-1, L) + 1;

    % 从起点 k0 开始，连续取 m 个（超出 L 自动回绕）
    idx = mod((k0-1) + (0:m-1), L) + 1;

    seg = P(idx, :);
end


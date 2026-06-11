function [detect,first_detect]=MC_line(fasan_D,fasan_d,path1,beam_vec,path2,f,omiga,UP,UV,init_beam_idx,R)
detect=0;
first_detect=nan;
% UV=[15,15,15];
N_KERNAL=100000;


UP_N=UP./norm(UP);

Beam_all=[beam_vec; path2; flipud(beam_vec); path1]; %第一发能否成功
if hit_linebeam_fast(UP_N, Beam_all(init_beam_idx,:), fasan_D, fasan_d)
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

    %==================== 每一时刻命中判据：线光(四棱锥) ====================%
    U_dir = P_now ./ norm_P_now;  % UAV 方向单位向量
    judge = hit_linebeam_fast(U_dir, B_now, fasan_D, fasan_d);
    % judge = hit_linebeam_pythonstyle(U_dir, B_now, fasan_D, fasan_d);

    judge(norm_P_now>R)=false; % 超出R强制置0
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

function in = hit_linebeam_fast(U_dir, B_dir, fasan_D, fasan_d)
% 四棱锥(矩形角锥)命中判据（推荐版，等价于你Python构造P2/P3/P4的几何）
% U_dir: N×3（可以不是单位向量，函数内会归一化）
% B_dir: N×3（beam指向，建议已是单位向量；函数内也会归一化）
% fasan_D/fasan_d: 两个正交方向的“全角”发散角（弧度），对应长边/短边

U_dir = U_dir ./ vecnorm(U_dir,2,2);
B_dir = B_dir ./ vecnorm(B_dir,2,2);

tanD = tan(fasan_D/2);  % = D/2（你Python里 D=2*tan(fasan_D/2)）
tand = tan(fasan_d/2);  % = d/2

bx = B_dir(:,1);
by = B_dir(:,2);
bxy = hypot(bx,by);

mask = (bxy < eps);
bxy_safe = bxy;
bxy_safe(mask) = 1;

% e1：你Python里的 vector_d（与beam垂直，且z=0方向定义一致）
e1 = [-by./bxy_safe, bx./bxy_safe, zeros(size(bx))];
if any(mask)
    e1(mask,:) = repmat([1 0 0], nnz(mask), 1);  % beam ~ z轴时给一个稳定的垂直方向
end

% e2：你Python里的 vector_D（= cross(beam, e1)）
e2 = cross(B_dir, e1, 2);

% 在beam坐标系下的分量
c3 = sum(U_dir .* B_dir, 2);  % 沿beam轴
c1 = sum(U_dir .* e1, 2);     % 短边方向
c2 = sum(U_dir .* e2, 2);     % 长边方向

% 四棱锥判据：|c2/c3| <= tanD 且 |c1/c3| <= tand，且必须在beam前方(c3>0)
in = (c3 > 0) & (abs(c2./c3) <= tanD) & (abs(c1./c3) <= tand);
end

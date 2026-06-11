
function [detect,first_time]=static_MC_func_ROI(N,f,fasan_D,fasan_d,z,omi_type,beam_type,path_type,jiaodu)

% N=1000;
v_min=1;
v_max=20;
% fasan_D=150e-3;
% fasan_d=1e-3;
% beam_type='line';
% omi_type='slow';
% path_type='raster';
% R=2000;
% f=5e3;
% jiaodu=25;
% z=500;
a= z*tand(jiaodu);
D = 2*z*tan(fasan_D/2);
d = 2*z*tan(fasan_d/2);
[UAV_P,~]=static_init_UAV_ROI(a,v_min,v_max,N);
%%
%========加载光束，填充首尾===========%
load_dir='E:\beam_vector\';
if strcmpi(omi_type,'fast')
    omiga=2*pi;
elseif strcmpi(omi_type,'slow')
    omiga=pi/2;
end
step_size=omiga*z./f;
R_ring_outer = D / 2;
R_ring_inner = R_ring_outer - d;
if R_ring_inner <= 0
    error('0.5D<d');
end


% scatter3(beam_vec_all(:,1),beam_vec_all(:,2),beam_vec_all(:,3));
%%
%==================开始仿真=====================%
detect = zeros(N,1);
first_detect = zeros(N,1); %初始化结果数组



tic
if strcmpi(beam_type,'ring')&&strcmpi(path_type,'spiral')
    r0 = a + 0.5*D; %这是直接甩出去的等距螺旋ROI轨迹
    if a>D
        path1 = generate_spiral_path(r0, D, step_size);
        path2 = generate_circular_path(r0, step_size,path1(end,1),path1(end,2));
        path = [path1; path2];
        path=path1; %无需同心圆轨迹，直接赋值即可
    else
        r0=0.5*D;
        path1=[r0,0];
        path2 = generate_circular_path(r0, step_size,path1(end,1),path1(end,2));
        path=[path1; path2];
    end

    for i=1:N
        px = UAV_P(i,1); py = UAV_P(i,2);
        detect(i)=static_MC_ring_ROI(fasan_D,fasan_d,path,f,[px,py],z);
    end
end
if strcmpi(beam_type,'ring')&&strcmpi(path_type,'raster')
    path = generate_path(D, step_size, a);

    for i=1:N
        px = UAV_P(i,1); py = UAV_P(i,2);
        detect(i)=static_MC_ring_ROI(fasan_D,fasan_d,path,f,[px,py],z);
    end
end
if strcmpi(beam_type,'line')&&strcmpi(path_type,'raster')
    path = generate_path(D, step_size, a);

    for i=1:N
        px = UAV_P(i,1); py = UAV_P(i,2);
        detect(i)=static_MC_line_ROI(fasan_D,fasan_d,path,f,[px,py],z);
    end
end
if strcmpi(beam_type,'line')&&strcmpi(path_type,'spiral')
    r0 = a + 0.5*D; %这是直接甩出去的等距螺旋ROI轨迹
    if a>D
        path1 = generate_spiral_path(r0, D, step_size);
        path2 = generate_circular_path(r0, step_size,path1(end,1),path1(end,2));
        path = [path1; path2];
        path=path1; %无需同心圆轨迹，直接赋值即可
    else
        r0=0.5*D;
        path1=[r0,0];
        path2 = generate_circular_path(r0, step_size,path1(end,1),path1(end,2));
        path=[path1; path2];
    end
    for i=1:N
        px= UAV_P(i,1); py = UAV_P(i,2);
        detect(i)=static_MC_line_ROI_SP(fasan_D,fasan_d,path,f,[px,py],z);
    end
end
toc

SUCCESS=sum(detect)/length(detect);
first_time=1;
disp(SUCCESS);
end


function detect=static_MC_ring_ROI(fasan_D,fasan_d,beam_vec_all,f,UP,z)
detect=0;
N_KERNAL=100000;
D = 2*z*tan(fasan_D/2);
d = 2*z*tan(fasan_d/2);
R_ring_outer = D / 2;
R_ring_inner = R_ring_outer - d;

P_now=[ones(N_KERNAL,1)*UP(1),ones(N_KERNAL,1)*UP(2)]; %初始化位置矩阵
B_now=[ones(N_KERNAL,1)*UP(1),ones(N_KERNAL,1)*UP(2)];

for i=1:ceil(length(beam_vec_all)/N_KERNAL)
    tmp=beam_vec_all((i-1)*N_KERNAL+1:min(N_KERNAL*i,length(beam_vec_all)),:);
    B_now(1:length(tmp),:)=tmp;

    norm_P_now=vecnorm(P_now,2,2); %每一时刻UAV位置与原点距离的平方
    d2 = sum((P_now - B_now).^2, 2);%UAV和当前光轴的距离平方，不开方更快
    judge=((d2<=R_ring_outer.^2)&(d2>=R_ring_inner.^2)); %判断标准

    if i==ceil(length(beam_vec_all)/N_KERNAL) %最后一次循环，将超出的无意义部分强制置0
        judge(1+mod(length(beam_vec_all),N_KERNAL):N_KERNAL)=0;
    end

    if any(judge)
        detect=1;
        % 此处后面再写，首次预警时间
        first_detect=(i-1)*N_KERNAL+find(judge,1,'first');
        
        return
    end
end

end



function detect=static_MC_line_ROI(fasan_D,fasan_d,beam_vec_all,f,UP,z)
detect=0;
N_KERNAL=100000;
D = 2*z*tan(fasan_D/2);
d = 2*z*tan(fasan_d/2);

P_now=[ones(N_KERNAL,1)*UP(1),ones(N_KERNAL,1)*UP(2)]; %初始化位置矩阵
B_now=[ones(N_KERNAL,1)*UP(1),ones(N_KERNAL,1)*UP(2)];

for i=1:ceil(length(beam_vec_all)/N_KERNAL)
    tmp=beam_vec_all((i-1)*N_KERNAL+1:min(N_KERNAL*i,length(beam_vec_all)),:);
    B_now(1:length(tmp),:)=tmp;

    norm_P_now=vecnorm(P_now,2,2); %每一时刻UAV位置与原点距离的平方
    d2 = sum((P_now - B_now).^2, 2);%UAV和当前光轴的距离平方，不开方更快
    ddx=abs(P_now(:,1)-B_now(:,1));
    ddy=abs(P_now(:,2)-B_now(:,2));

    % judge=((d2<=R_ring_outer.^2)&(d2>=R_ring_inner.^2)); %判断标准
    judge=((ddx<D/2)&(ddy<d/2)); %判断标准

    if i==ceil(length(beam_vec_all)/N_KERNAL) %最后一次循环，将超出的无意义部分强制置0
        judge(1+mod(length(beam_vec_all),N_KERNAL):N_KERNAL)=0;
    end

    if any(judge)
        detect=1;
        % 此处后面再写，首次预警时间
        first_detect=(i-1)*N_KERNAL+find(judge,1,'first');
        
        return
    end
end

end

function detect = static_MC_line_ROI_SP(fasan_D, fasan_d, beam_vec_all, f, UP, z)
% static_MC_line_ROI_SP
% 用于静态ROI场景下，线光场采用spiral轨迹单次扫描时的命中判定。
%
% 判定规则：
%   每一发脉冲对应一个矩形线光斑；
%   其中线光斑“长边”始终与该轨迹点处的螺旋切向垂直；
%   “短边”沿螺旋切向。
%
% 输入：
%   fasan_D      线光长边发散角(rad)
%   fasan_d      线光短边发散角(rad)
%   beam_vec_all N×2，spiral轨迹上各脉冲中心点 [x, y]
%   f            脉冲频率（本函数中仅保留接口一致性，不参与运算）
%   UP           1×2，静态目标位置 [x, y]
%   z            探测距离
%
% 输出：
%   detect       0/1，单次扫描内是否命中
%
% 说明：
%   假设 beam_vec_all 来自阿基米德螺旋：
%       r = b * theta,   b = pitch / (2*pi)
%   而在你当前程序里，generate_spiral_path(r0, D, step_size) 中的 pitch 就是 D，
%   因此这里直接取 pitch = D。
%
%   阿基米德螺旋切向解析式：
%       t = dr/dtheta * e_r + r * e_theta = b*e_r + r*e_theta
%   若轨迹点为 (x,y)，r = hypot(x,y)，则其笛卡尔形式为：
%       t = (b/r) * [x, y] + [-y, x]
%
%   再将单位切向记为 e_short，则长边方向单位向量为
%       e_long = [-e_short(2), e_short(1)]
%
%   命中判据：
%       |rel · e_long | <= D/2
%       |rel · e_short| <= d/2

    %#ok<*NASGU>
    detect = 0;
    N_KERNAL = 100000;

    % 保持与原接口一致
    % f 在本静态判定函数中不显式使用
    if size(beam_vec_all,2) ~= 2
        error('beam_vec_all 必须为 N×2 矩阵。');
    end
    if numel(UP) ~= 2
        error('UP 必须为二维坐标 [x,y]。');
    end

    UP = reshape(UP, 1, 2);

    % 实际光斑尺寸
    D = 2 * z * tan(fasan_D / 2);   % 长边尺寸
    d = 2 * z * tan(fasan_d / 2);   % 短边尺寸

    % 这里直接按你 generate_spiral_path(..., D, ...) 的定义，
    % 将螺旋节距 pitch 取为 D
    pitch = D;
    b = pitch / (2*pi);

    % ------------------ 预计算每个轨迹点的局部坐标系 ------------------ %
    x = beam_vec_all(:,1);
    y = beam_vec_all(:,2);
    r = hypot(x, y);

    % 螺旋切向解析式：
    % t = (b/r)[x,y] + [-y,x]
    tx = zeros(size(x));
    ty = zeros(size(y));

    mask = (r > eps);
    tx(mask) = (b ./ r(mask)) .* x(mask) - y(mask);
    ty(mask) = (b ./ r(mask)) .* y(mask) + x(mask);

    % r -> 0 时的极限方向：阿基米德螺旋起点切向沿 +x
    tx(~mask) = 1;
    ty(~mask) = 0;

    nt = hypot(tx, ty);
    e_short = [tx ./ nt, ty ./ nt];                    % 短边方向：沿切向
    e_long  = [-e_short(:,2), e_short(:,1)];          % 长边方向：垂直切向

    % ------------------ 分块判定，避免一次性大数组 ------------------ %
    Npulse = size(beam_vec_all,1);

    for i = 1:ceil(Npulse / N_KERNAL)
        idx1 = (i-1)*N_KERNAL + 1;
        idx2 = min(i*N_KERNAL, Npulse);
        idx  = idx1:idx2;

        B_now  = beam_vec_all(idx, :);
        EL_now = e_long(idx, :);
        ES_now = e_short(idx, :);

        % 目标相对各脉冲中心的位置
        rel = UP - B_now;   % 自动扩展为 n×2

        % 投影到长边/短边坐标轴
        u = sum(rel .* EL_now, 2);   % 长边方向投影
        v = sum(rel .* ES_now, 2);   % 短边方向投影

        % 旋转矩形命中判据
        judge = (abs(u) <= D/2) & (abs(v) <= d/2);

        if any(judge)
            detect = 1;
            return
        end
    end
end




function detect=static_MC_line(fasan_D,fasan_d,UP,beam_vec_all,f)

% UAV_COUNTS=length(UP);
detect=0;
N_KERNAL=100000; % 单次仿真最大值
UP_N=UP./norm(UP);
first_detect=nan;

P_now=[ones(N_KERNAL,1)*UP(1),ones(N_KERNAL,1)*UP(2),ones(N_KERNAL,1)*UP(3)];
B_now=[ones(N_KERNAL,1)*UP(1),ones(N_KERNAL,1)*UP(2),ones(N_KERNAL,1)*UP(3)];

for i=1:ceil(length(beam_vec_all)/N_KERNAL)

    tmp=beam_vec_all((i-1)*N_KERNAL+1:min(N_KERNAL*i,length(beam_vec_all)),:);
    B_now(1:length(tmp),:)=tmp;

    norm_P_now=vecnorm(P_now,2,2); %每一时刻UAV位置的单位向量
    %==================== 每一时刻命中判据：线光(四棱锥) ====================%
    U_dir = P_now ./ norm_P_now;  % UAV 方向单位向量
    judge = hit_linebeam_fast(U_dir, B_now, fasan_D, fasan_d);


    if i==ceil(length(beam_vec_all)/N_KERNAL) %最后一次循环，将超出的无意义部分强制置0
        judge(1+mod(length(beam_vec_all),N_KERNAL):N_KERNAL)=0;
    end

    if any(judge)
        detect=1;
        % 此处后面再写，首次预警时间
        first_detect=(i-1)*N_KERNAL+find(judge,1,'first');

        return
    end
end

end



%%
%=========轨迹生成函数===========%
function path = generate_spiral_path(r_end, pitch, step_size)
a = pitch / (2 * pi);
theta_max = r_end / a;
L = 0.5 * a * (theta_max * sqrt(theta_max^2 + 1) + log(theta_max + sqrt(theta_max^2 + 1)));
n_points = ceil(L / step_size) + 1;
path = zeros(n_points, 2);
path(1, :) = [0, 0];
point_count = 1;
theta = 0;
while true
    ds_dtheta = a * sqrt(theta^2 + 1);
    dtheta = step_size / ds_dtheta;
    theta_new = theta + dtheta;
    r_new = a * theta_new;
    if r_new >= r_end
        theta_end = r_end / a;
        x_end = r_end * cos(theta_end);
        y_end = r_end * sin(theta_end);
        point_count = point_count + 1;
        path(point_count, :) = [x_end, y_end];
        break;
    end
    point_count = point_count + 1;
    path(point_count, :) = [r_new * cos(theta_new), r_new * sin(theta_new)];
    theta = theta_new;
end

path = path(1:point_count, :);
end

function path = generate_circular_path(radius, step_size,x_start,y_start)
arc_step = step_size / radius;
[theta_start,~]=cart2pol(x_start,y_start);
theta = theta_start:arc_step:theta_start+2*pi;
x = radius * cos(theta);
y = radius * sin(theta);
path = [x(:), y(:)];
end

% function path=generate_path(D,d,a)
% % 尖刺RASTER轨迹
% if D>=2*a
% 
%     y=-a-d:d:a+d;
%     x=0.*ones(length(y),1);
%     path=[x(:),y(:)];
% else
%     step_count_single=ceil(2*a/d);
%     path=[-a+D/2,-a];
%     for i=1:ceil(2*a/D)
%         x0=-(a-D/2)+((2*a-D)/(ceil(2*a/D)-1))*(i-1);
%         x0=-a+D/2+D*(i-1); %无重叠时每个条带的起始中心点x坐标
%         % if mod(i,2)==1 %判断是否奇偶，横平竖直时奇偶反向
%         %     y0=-a;
%         %     %             y=y0:d:-y0;
%         %     y=linspace(y0,-y0,ceil(-2*y0/d));
%         % else
%         %     y0=a;
%         %     y=y0:-d:-y0;
%         %     y=linspace(y0,-y0,ceil(2*y0/d));
%         % end
%         y0=-a;
%         y=linspace(y0,-y0,ceil(-2*y0/d));
%         x=x0.*ones(step_count_single,1);
%         path=[path;x(:),y(:)];
%         if i~=ceil(2*a/D)
%             % 以下3行为原横平竖直链接线代码
%             % x=x0:d:x0+((2*a-D)/(ceil(2*a/D)-1));
%             % y=-y0.*ones(length(x),1);
%             % path=[path;x(:),y(:)];
%             % 当前扫描线终点 -> 下一条扫描线起点
%             p_start=[x0,a];
%             p_end=[x0+D,-a];
%             % 连接线总长度
%             L = hypot(p_end(1)-p_start(1), p_end(2)-p_start(2));
% 
%             % 按接近 d 的间距进行等间隔采样
%             n_conn = max(1, round(L/d));      % 分成 n_conn 段
%             t = linspace(0, 1, n_conn+1).';   % 参数 0~1
% 
%             x_conn = p_start(1) + (p_end(1)-p_start(1)) * t;
%             y_conn = p_start(2) + (p_end(2)-p_start(2)) * t;
% 
%             % 去掉起点，避免与当前扫描线终点重复
%             % 若下一条扫描线本身包含起点 [x1,-a]，则终点也去掉
%             path = [path; x_conn(2:end-1), y_conn(2:end-1)];
% 
%         end
%     end
% end
% end


function path=generate_path(D,d,a)
if D>=2*a

    y=-a-d:d:a+d;
    x=0.*ones(length(y),1);
    path=[x(:),y(:)];
else
    step_count_single=ceil(2*a/d);
    path=[-a+D/2,-a];
    for i=1:ceil(2*a/D)
        % x0=-(a-D/2)+((2*a-D)/(ceil(2*a/D)-1))*(i-1);
        x0=-a+D/2+D*(i-1); %无重叠时每个条带的起始中心点x坐标
        if mod(i,2)==1%����
            y0=-a;
            %             y=y0:d:-y0;
            y=linspace(y0,-y0,ceil(-2*y0/d));
        else
            y0=a;
            y=y0:-d:-y0;
            y=linspace(y0,-y0,ceil(2*y0/d));
        end
        x=x0.*ones(step_count_single,1);
        path=[path;x(:),y(:)];
        if i~=ceil(2*a/D)
            x=x0:d:x0+((2*a-D)/(ceil(2*a/D)-1));
            y=-y0.*ones(length(x),1);
            path=[path;x(:),y(:)];
        end
    end
end
end
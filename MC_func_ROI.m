function [detect,first_time]=MC_func_ROI(z,N,f,jiaodu,beam_type,fasan_D,fasan_d,v_min,v_max,omi_type,path_type)
%%
%=========参数与初始化===========%
% N=5000;
% f=5e3;
% omi_type='fast';
% beam_type='ring';
% path_type='raster';
% fasan_point=15e-3;
% fasan_D=150e-3;
% fasan_d=1e-3;
% v_max=40;
% v_min=0.5;
% z=500;
a= z*tand(jiaodu);
D = 2*z*tan(fasan_D/2);
d = 2*z*tan(fasan_d/2);
[UAV_P,UAV_V]=init_UAV_ROI(a,v_min,v_max,N); %无人机初始化

%%
%========加载光束，填充首尾===========%

if omi_type=='fast'
    omiga=2*pi;
elseif omi_type=='slow'
    omiga=pi/2;
end
step_size=omiga*z./f;
R_ring_outer = D / 2;
R_ring_inner = R_ring_outer - d;
if R_ring_inner <= 0
    error('0.5D<d');
end

% ===== 并行池：只创建一次（MC_func多次调用时复用），避免反复启动开销 =====
useParallel = license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
if useParallel
    p = gcp('nocreate');
    if isempty(p)
        % 优先线程池（启动快、共享内存；对你这种大只读数据广播更友好）
        try
            parpool('threads');
        catch
            % 旧版本或不支持 threads 时回退到默认 local 进程池
            parpool('local');
        end
    end
end
detect = zeros(N,1);
first_detect = zeros(N,1);


tic
if useParallel
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
        beam_idx_rand = randi([1,length(path)], N, 1);

        parfor i=1:N
            px = UAV_P(i,1); py = UAV_P(i,2);
            vx = UAV_V(i,1); vy = UAV_V(i,2);
            init_beam_idx=beam_idx_rand(i);
            % [detect(i), first_detect(i)]=MC_ring_ROI_SP(fasan_D,fasan_d,path1,path2,f,omiga,[px,py],[vx,vy],init_beam_idx,z,a);
            [detect(i), first_detect(i)]=MC_ring_ROI_SP_snr(fasan_D,fasan_d,path1,[],f,omiga,[px,py],[vx,vy],init_beam_idx,z,a);
        end
    end
    if strcmpi(beam_type,'ring')&&strcmpi(path_type,'raster')
        path = generate_path(D, step_size, a);
        beam_idx_rand = randi([1,length(path)], N, 1);

        parfor i=1:N
            px = UAV_P(i,1); py = UAV_P(i,2);
            vx = UAV_V(i,1); vy = UAV_V(i,2);
            init_beam_idx=beam_idx_rand(i);
            [detect(i), first_detect(i)]=MC_ring_ROI_RA_snr(fasan_D,fasan_d,path,f,omiga,[px,py],[vx,vy],init_beam_idx,z,a);
        end
    end
    if strcmpi(beam_type,'line')&&strcmpi(path_type,'raster')
        path = generate_path(D, step_size, a);
        beam_idx_rand = randi([1,length(path)], N, 1);

        parfor i=1:N
            px = UAV_P(i,1); py = UAV_P(i,2);
            vx = UAV_V(i,1); vy = UAV_V(i,2);
            init_beam_idx=beam_idx_rand(i);
            [detect(i), first_detect(i)]=MC_line_ROI_RA(fasan_D,fasan_d,path,f,omiga,[px,py],[vx,vy],init_beam_idx,z,a);
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
        beam_idx_rand = randi([1,length(path)], N, 1);

        parfor i=1:N
            px = UAV_P(i,1); py = UAV_P(i,2);
            vx = UAV_V(i,1); vy = UAV_V(i,2);
            init_beam_idx=beam_idx_rand(i);
            % [detect(i), first_detect(i)]=MC_ring_ROI_SP(fasan_D,fasan_d,path1,path2,f,omiga,[px,py],[vx,vy],init_beam_idx,z,a);
            [detect(i), first_detect(i)]=MC_line_ROI_SP(fasan_D,fasan_d,path1,[],f,omiga,[px,py],[vx,vy],init_beam_idx,z,a);
        end
    end
else
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
        beam_idx_rand = randi([1,length(path)], N, 1);

        for i=1:N
            px = UAV_P(i,1); py = UAV_P(i,2);
            vx = UAV_V(i,1); vy = UAV_V(i,2);
            init_beam_idx=beam_idx_rand(i);
            % [detect(i), first_detect(i)]=MC_ring_ROI_SP(fasan_D,fasan_d,path1,path2,f,omiga,[px,py],[vx,vy],init_beam_idx,z,a);
            [detect(i), first_detect(i)]=MC_ring_ROI_SP(fasan_D,fasan_d,path1,[],f,omiga,[px,py],[vx,vy],init_beam_idx,z,a);
        end
    end
    if strcmpi(beam_type,'ring')&&strcmpi(path_type,'raster')
        path = generate_path(D, step_size, a);
        beam_idx_rand = randi([1,length(path)], N, 1);

        for i=1:N
            px = UAV_P(i,1); py = UAV_P(i,2);
            vx = UAV_V(i,1); vy = UAV_V(i,2);
            init_beam_idx=beam_idx_rand(i);
            [detect(i), first_detect(i)]=MC_ring_ROI_RA(fasan_D,fasan_d,path,f,omiga,[px,py],[vx,vy],init_beam_idx,z,a);
        end
    end
    if strcmpi(beam_type,'line')&&strcmpi(path_type,'raster')
        path = generate_path(D, step_size, a);
        beam_idx_rand = randi([1,length(path)], N, 1);

        for i=1:N
            px = UAV_P(i,1); py = UAV_P(i,2);
            vx = UAV_V(i,1); vy = UAV_V(i,2);
            init_beam_idx=beam_idx_rand(i);
            [detect(i), first_detect(i)]=MC_line_ROI_RA(fasan_D,fasan_d,path,f,omiga,[px,py],[vx,vy],init_beam_idx,z,a);
        end
    end
end
toc
SUCCESS=sum(detect)/length(detect);
avetime=mean(first_detect,'omitnan')/f;
first_time=first_detect./f;
disp([omi_type,num2str(fasan_D),'_',num2str(fasan_d),'_',num2str(SUCCESS),'_',num2str(avetime)]);

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

% 以下为带重叠的横平竖直型raster轨迹
% function path=generate_path(D,d,a)
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
%         if mod(i,2)==1%����
%             y0=-a;
%             %             y=y0:d:-y0;
%             y=linspace(y0,-y0,ceil(-2*y0/d));
%         else
%             y0=a;
%             y=y0:-d:-y0;
%             y=linspace(y0,-y0,ceil(2*y0/d));
%         end
%         x=x0.*ones(step_count_single,1);
%         path=[path;x(:),y(:)];
%         if i~=ceil(2*a/D)
%             x=x0:d:x0+((2*a-D)/(ceil(2*a/D)-1));
%             y=-y0.*ones(length(x),1);
%             path=[path;x(:),y(:)];
%         end
%     end
% end
% end

% function path=generate_path(D,d,a)
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
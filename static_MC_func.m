
function [detect,first_time]=static_MC_func(N,f,fasan_D,fasan_d,R,omi_type,beam_type)
% clc;clear all;close all;
% N=1000;
v_min=1;
v_max=20;
% fasan_D=150e-3;
% fasan_d=1e-3;
% beam_type='line';
% omi_type='fast';
% R=2000;
% f=5e3;
[UAV_P,~]=init_UAV(R,v_min,v_max,N);
%%
%========加载光束，填充首尾===========%
load_dir='D:\pythonandconda\MonteCarlo\path_pi2_2pi\';
if strcmpi(omi_type,'fast')
    omiga=2*pi;
elseif strcmpi(omi_type,'slow')
    omiga=pi/2;
end
beam_vec=readNPY([load_dir,'beam_vector_',omi_type,'_',num2str(f),'Hz_h3000m_D',num2str(fasan_D*1e6),'u.npy']);
x0=beam_vec(1,1);y0=beam_vec(1,2);z0=beam_vec(1,3);
theta0=atan2(y0,x0);r0=sqrt(x0^2+y0^2);
theta=(theta0:omiga/f/r0:theta0+2*pi)';
x=r0*cos(theta);y=r0*sin(theta);z=ones(length(x),1)*z0;
path1=[x,y,z];
x0=beam_vec(end,1);y0=beam_vec(end,2);z0=beam_vec(end,3);
theta0=atan2(y0,x0);r0=sqrt(x0^2+y0^2);
theta=(theta0:omiga/f/r0:theta0+2*pi)';
x=r0*cos(theta);y=r0*sin(theta);z=ones(length(x),1)*z0;
path2=[x,y,z];
beam_vec_all=[path1;beam_vec;path2];

scatter3(beam_vec_all(:,1),beam_vec_all(:,2),beam_vec_all(:,3));
%%
%==================开始仿真=====================%
detect = zeros(N,1); %初始化结果数组

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


tic
if useParallel
    parfor i=1:N
        px = UAV_P(i,1); py = UAV_P(i,2); pz = UAV_P(i,3);
        if strcmpi(beam_type,'ring')
            detect(i)=static_MC_ring(fasan_D,fasan_d,[px,py,pz],beam_vec_all,f);
        end
        if strcmpi(beam_type,'line')
            detect(i)=static_MC_line(fasan_D,fasan_d,[px,py,pz],beam_vec_all,f);
        end
    end
else
    for i=1:N
        px = UAV_P(i,1); py = UAV_P(i,2); pz = UAV_P(i,3);
        if strcmpi(beam_type,'ring')
            detect(i)=static_MC_ring(fasan_D,fasan_d,[px,py,pz],beam_vec_all,f);
        end
        if strcmpi(beam_type,'line')
            detect(i)=static_MC_line(fasan_D,fasan_d,[px,py,pz],beam_vec_all,f);
        end
    end
end
toc


% if strcmpi(beam_type,'ring')
%     for i=1:N
%         px = UAV_P(i,1); py = UAV_P(i,2); pz = UAV_P(i,3);
%         detect(i)=static_MC_ring(fasan_D,fasan_d,[px,py,pz],beam_vec_all,f);
%     end
% end
% if strcmpi(beam_type,'line')
%     for i=1:N
%         px = UAV_P(i,1); py = UAV_P(i,2); pz = UAV_P(i,3);
%         detect(i)=static_MC_line(fasan_D,fasan_d,[px,py,pz],beam_vec_all,f);
%     end
% end








SUCCESS=sum(detect)/length(detect);
first_time=1;
disp(SUCCESS);
end

function detect=static_MC_ring(fasan_D,fasan_d,UP,beam_vec_all,f)

% UAV_COUNTS=length(UP);
detect=0;
N_KERNAL=100000; % 单次仿真最大值
first_detect=nan;
UAV_COUNTS=0;

P_now=[ones(N_KERNAL,1)*UP(1),ones(N_KERNAL,1)*UP(2),ones(N_KERNAL,1)*UP(3)];
B_now=[ones(N_KERNAL,1)*UP(1),ones(N_KERNAL,1)*UP(2),ones(N_KERNAL,1)*UP(3)];

for i=1:ceil(length(beam_vec_all)/N_KERNAL)

    tmp=beam_vec_all((i-1)*N_KERNAL+1:min(N_KERNAL*i,length(beam_vec_all)),:);
    B_now(1:length(tmp),:)=tmp;


    norm_P_now=vecnorm(P_now,2,2); %每一时刻UAV位置的单位向量
    dtheta=acos(dot(P_now,B_now,2)./norm_P_now); %UAV位置和BEAM指向的夹角
    judge=(((fasan_D/2-fasan_d)<=dtheta)&(dtheta<=fasan_D/2)); %判断标准


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
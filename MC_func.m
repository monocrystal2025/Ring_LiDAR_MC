function [detect,first_time]=MC_func(R,N,f,beam_type,fasan_D,fasan_d,v_min,v_max,omi_type)
%%
%=========参数与初始化===========%
% N=5000;
% f=20e3;
% omi_type='slow';
% beam_type='ring';
% fasan_point=15e-3;
% fasan_D=150e-3;
% fasan_d=0.5*fasan_D-(sqrt(fasan_D^2-fasan_point^2))/2;
% v_max=40;
% v_min=0.5;
% R=2000;
[UAV_P,UAV_V]=init_UAV(R,v_min,v_max,N); %无人机初始化



%%
%========加载光束，填充首尾===========%
load_dir='G:\BeamVEC_NEW\';
if strcmpi(omi_type, 'fast')
    omiga=2*pi;
elseif strcmpi(omi_type, 'slow')
    omiga=pi/2;
else
    error('MC_func:UnknownOmiType', ...
        'Unknown omi_type: %s. Use fast or slow.', omi_type);
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
%==========MC===========%
% beam_idx_rand=randi([1,length(beam_vec_all)],N,1);
% direction_rand=randi([0,1],N,1)*2-1;
% detect=zeros(N,1);
% first_detect=zeros(N,1);
% 
% 
% % % ---- 进度条初始化 ----
% % hwb = waitbar(0, '0%', 'Name', 'MC Progress', 'CreateCancelBtn', '');
% % tStart = tic;
% % updateStep = max(1, round(N/100));   % 每 1% 更新一次（N很小时至少每次更新）
% tic
% for i=1:N
%     px=UAV_P(i,1);py=UAV_P(i,2);pz=UAV_P(i,3);
%     vx=UAV_V(i,1);vy=UAV_V(i,2);vz=UAV_V(i,3);
%     beam_idx=beam_idx_rand(i);
%     direction=direction_rand(i);
%     if strcmpi(beam_type,'ring')
%         [detect(i),first_detect(i)]=MC_ring(fasan_D,fasan_d, ...
%             path1,beam_vec,path2, ...
%             f, omiga, [px,py,pz], [vx,vy,vz], beam_idx, R);
%     end
%     if strcmpi(beam_type,'point')
%         [detect(i),first_detect(i)]=MC_point(fasan_D,fasan_d, ...
%             path1,beam_vec,path2, ...
%             f, omiga, [px,py,pz], [vx,vy,vz], beam_idx, R);
%     end
%     if strcmpi(beam_type,'line')
%         [detect(i),first_detect(i)]=MC_line(fasan_D,fasan_d, ...
%             path1,beam_vec,path2, ...
%             f, omiga, [px,py,pz], [vx,vy,vz], beam_idx, R);
%     end
% 
%     % % ---- 进度条更新（控制刷新频率，减少开销）----
%     % if mod(i, updateStep) == 0 || i == N
%     %     elapsed = toc(tStart);
%     %     eta = elapsed * (N - i) / max(i,1);
%     %     msg = sprintf('%.1f%%  (%d/%d)  Elapsed: %.1fs  ETA: %.1fs', ...
%     %         100*i/N, i, N, elapsed, eta);
%     %     if ishandle(hwb)
%     %         waitbar(i/N, hwb, msg);
%     %     end
%     % end
% 
% 
% 
% end
% toc



%%
%==========MC===========%
beam_idx_rand = randi([1,length(beam_vec_all)], N, 1);
direction_rand = randi([0,1], N, 1)*2 - 1;
detect = zeros(N,1);
first_detect = zeros(N,1);

isRing = strcmpi(beam_type,'ring');
isPoint = strcmpi(beam_type,'point');
isLine = strcmpi(beam_type,'line');
if ~(isRing || isPoint || isLine)
    error('MC_func:UnknownBeamType', ...
        'Unknown beam_type: %s. Use ring, point, or line.', beam_type);
end

UAV_P_x = UAV_P(:,1); UAV_P_y = UAV_P(:,2); UAV_P_z = UAV_P(:,3);
UAV_V_x = UAV_V(:,1); UAV_V_y = UAV_V(:,2); UAV_V_z = UAV_V(:,3);

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
    parfor i = 1:N
        px = UAV_P_x(i); py = UAV_P_y(i); pz = UAV_P_z(i);
        vx = UAV_V_x(i); vy = UAV_V_y(i); vz = UAV_V_z(i);
        beam_idx = beam_idx_rand(i);
        direction = direction_rand(i); %#ok<NASGU>  % 你目前没用到direction，先保留不改

        if isRing
            [detect(i), first_detect(i)] = MC_ring_snr(fasan_D,fasan_d, ...
                path1,beam_vec,path2, ...
                f, omiga, [px,py,pz], [vx,vy,vz], beam_idx, R);
        elseif isPoint
            [detect(i), first_detect(i)] = MC_point_snr(fasan_D,fasan_d, ...
                path1,beam_vec,path2, ...
                f, omiga, [px,py,pz], [vx,vy,vz], beam_idx, R);
        elseif isLine
            [detect(i), first_detect(i)] = MC_line_snr(fasan_D,fasan_d, ...
                path1,beam_vec,path2, ...
                f, omiga, [px,py,pz], [vx,vy,vz], beam_idx, R);
        end
    end
else
    for i = 1:N
        px = UAV_P_x(i); py = UAV_P_y(i); pz = UAV_P_z(i);
        vx = UAV_V_x(i); vy = UAV_V_y(i); vz = UAV_V_z(i);
        beam_idx = beam_idx_rand(i);
        direction = direction_rand(i); %#ok<NASGU>

        if isRing
            [detect(i), first_detect(i)] = MC_ring_snr(fasan_D,fasan_d, ...
                path1,beam_vec,path2, ...
                f, omiga, [px,py,pz], [vx,vy,vz], beam_idx, R);
        elseif isPoint
            [detect(i), first_detect(i)] = MC_point_snr(fasan_D,fasan_d, ...
                path1,beam_vec,path2, ...
                f, omiga, [px,py,pz], [vx,vy,vz], beam_idx, R);
        elseif isLine
            [detect(i), first_detect(i)] = MC_line_snr(fasan_D,fasan_d, ...
                path1,beam_vec,path2, ...
                f, omiga, [px,py,pz], [vx,vy,vz], beam_idx, R);
        end
    end
end
toc









SUCCESS=sum(detect)/length(detect);
avetime=mean(first_detect,'omitnan')/f;
disp([omi_type,num2str(fasan_D),'_',num2str(fasan_d),'_',num2str(SUCCESS),'_',num2str(avetime)]);
first_time=first_detect./f;
end

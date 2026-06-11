clc;clear;close all;
%=========参数与初始化===========%
N=50000;
f=5e3;
omi_type='slow';
beam_type='ring';
fasan_point=15e-3;
fasan_D_LIST=5e-3:5e-3:200e-3;
v_max=[20,40];
v_min=0.5;
R=2000;
detect_R=zeros(N,2);detect_L=zeros(N,2);
first_time_R=zeros(N,2);first_time_L=zeros(N,2);
% [a,b]=MC_func(R,N,f,'ring',150e-3,1.5e-3,v_min,v_max(1),'fast');
path='D:\matlab\STATIC_MC_RESULTS\';

for i=length(fasan_D_LIST):-1:5
    tic
    fasan_D=fasan_D_LIST(i);
    fasan_d=10e-3;
    fasan_d=fasan_D/2;
    [a,b]=static_MC_func(N,f,fasan_D,fasan_d,R,'fast','ring');
    detect_R(:,1)=a;first_time_R(:,1)=b;
    % [a,b]=static_MC_func(N,f,fasan_D,fasan_d,R,'fast','ring');
    detect_R(:,2)=a;first_time_R(:,2)=b;
    save([path,'STATIC_MC_1par_EA_RING_D',num2str(fasan_D*1000),'d',num2str(fasan_d*1000),'mrad.mat'],"detect_R","first_time_R");

    % [a,b]=static_MC_func(N,f,fasan_D,fasan_d,R,'fast','line');
    % detect_L(:,1)=a;first_time_L(:,1)=b;
    % % [a,b]=static_MC_func(N,f,fasan_D,fasan_d,R,'fast','line');
    % detect_L(:,2)=a;first_time_L(:,2)=b;
    % save([path,'STATIC_MC_1par_EA_LINE_D',num2str(fasan_D*1000),'d',num2str(fasan_d*1000),'mrad.mat'],"detect_L","first_time_L");
    toc



end
% fasan_D=5e-3;
% fasan_d=2.5e-3;
% [a,b]=static_MC_func(N,f,fasan_D,fasan_d,R,'fast','ring');
% detect_R(:,1)=a;first_time_R(:,1)=b;
% % [a,b]=static_MC_func(N,f,fasan_D,fasan_d,R,'fast','ring');
% detect_R(:,2)=a;first_time_R(:,2)=b;
% save([path,'STATIC_MC_1par_EA_RING_D',num2str(fasan_D*1000),'d',num2str(fasan_d*1000),'mrad.mat'],"detect_R","first_time_R");
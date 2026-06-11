clc;clear;close all;
%=========参数与初始化===========%
N=50000;
f=5e3;
omi_type='slow';
beam_type='ring';
fasan_point=15e-3;
fasan_D_LIST=5e-3:5e-3:180e-3;
v_max=[20,40];
v_min=1;
R=2000;
detect_R=zeros(N,4);detect_L=zeros(N,4);
first_time_R=zeros(N,4);first_time_L=zeros(N,4);
% [a,b]=MC_func(R,N,f,'line',150e-3,1e-3,v_min,v_max(1),'fast');
path='D:\matlab\MC_RESULTS_SNR\';

for i=length(fasan_D_LIST):-1:1
    fasan_D=fasan_D_LIST(i);
    % fasan_d=0.5*fasan_D-(sqrt(fasan_D^2-fasan_point^2))/2;
    fasan_d=0.8e-3;
    % fasan_d=fasan_D/2;
    [a,b]=MC_func(R,N,f,'ring',fasan_D,fasan_d,v_min,30,'fast');
    detect_R(:,1)=a;first_time_R(:,1)=b;
    % [a,b]=MC_func(R,N,f,'ring',fasan_D,fasan_d,v_min,v_max(1),'fast');
    detect_R(:,2)=a;first_time_R(:,2)=b;
    % [a,b]=MC_func(R,N,f,'ring',fasan_D,fasan_d,v_min,v_max(2),'slow');
    detect_R(:,3)=a;first_time_R(:,3)=b;
    % [a,b]=MC_func(R,N,f,'ring',fasan_D,fasan_d,v_min,v_max(2),'fast');
    detect_R(:,4)=a;first_time_R(:,4)=b;
    save([path,'MC_1par_EA_RING_D',num2str(fasan_D*1000),'d',num2str(fasan_d*1000),'mrad.mat'],"detect_R","first_time_R");

    % fasan_d=(0.25*pi*fasan_point^2)/fasan_D;
    % fasan_d=0.8e-3;
    % [a,b]=MC_func(R,N,f,'line',fasan_D,fasan_d,v_min,30,'fast');
    % detect_L(:,1)=a;first_time_L(:,1)=b;
    % % [a,b]=MC_func(R,N,f,'line',fasan_D,fasan_d,v_min,v_max(1),'fast');
    % detect_L(:,2)=a;first_time_L(:,2)=b;
    % % [a,b]=MC_func(R,N,f,'line',fasan_D,fasan_d,v_min,v_max(2),'slow');
    % detect_L(:,3)=a;first_time_L(:,3)=b;
    % % [a,b]=MC_func(R,N,f,'line',fasan_D,fasan_d,v_min,v_max(2),'fast');
    % detect_L(:,4)=a;first_time_L(:,4)=b;
    % save([path,'MC_1par_EA_LINE_D',num2str(fasan_D*1000),'d',num2str(fasan_d*1000),'mrad.mat'],"detect_L","first_time_L");




end
clc;clear;close all;
%=========参数与初始化===========%
N=3000;
f=5e3;
omi_type='slow';
beam_type='ring';
fasan_point=15e-3;
fasan_D_LIST=25e-3:10e-3:390e-3;
v_max=[20,40];
v_min=1;
R=2000;
detect_R=zeros(N,4);detect_L=zeros(N,4);detect_P=zeros(N,4);
first_time_R=zeros(N,4);first_time_L=zeros(N,4);first_time_P=zeros(N,4);
first_encounter_time_R=zeros(N,4);first_encounter_time_L=zeros(N,4);first_encounter_time_P=zeros(N,4);
effective_pulses_R=zeros(N,4);effective_pulses_L=zeros(N,4);effective_pulses_P=zeros(N,4);
% [a,b,c,d]=MC_func(R,N,f,'line',150e-3,1e-3,v_min,v_max(1),'fast');
path='D:\lzx\MC_RESULTS\';

for i=length(fasan_D_LIST):-1:1
    fasan_D=fasan_D_LIST(i);
    % fasan_d=0.5*fasan_D-(sqrt(fasan_D^2-fasan_point^2))/2;
    fasan_d=1e-3;
    % fasan_d=fasan_D/2;
    [a,b,c,d]=MC_func(R,N,f,'ring',fasan_D,fasan_d,v_min,30,'fast');
    detect_R(:,1)=a;first_time_R(:,1)=b;first_encounter_time_R(:,1)=c;effective_pulses_R(:,1)=d;
    % % [a,b,c,d]=MC_func(R,N,f,'ring',fasan_D,fasan_d,v_min,v_max(1),'fast');
    % detect_R(:,2)=a;first_time_R(:,2)=b;first_encounter_time_R(:,2)=c;effective_pulses_R(:,2)=d;
    % % [a,b,c,d]=MC_func(R,N,f,'ring',fasan_D,fasan_d,v_min,v_max(2),'slow');
    % detect_R(:,3)=a;first_time_R(:,3)=b;first_encounter_time_R(:,3)=c;effective_pulses_R(:,3)=d;
    % % [a,b,c,d]=MC_func(R,N,f,'ring',fasan_D,fasan_d,v_min,v_max(2),'fast');
    % detect_R(:,4)=a;first_time_R(:,4)=b;first_encounter_time_R(:,4)=c;effective_pulses_R(:,4)=d;
    % save([path,'MC_1par_EA_RING_D',num2str(fasan_D*1000),'d',num2str(fasan_d*1000),'mrad.mat'],"detect_R","first_time_R","first_encounter_time_R","effective_pulses_R");

    % fasan_d=(0.25*pi*fasan_point^2)/fasan_D;
    % fasan_d=0.8e-3;
    [a,b,c,d]=MC_func(R,N,f,'line',fasan_D,fasan_d,v_min,30,'fast');
    detect_L(:,1)=a;first_time_L(:,1)=b;first_encounter_time_L(:,1)=c;effective_pulses_L(:,1)=d;
    % [a,b,c,d]=MC_func(R,N,f,'line',fasan_D,fasan_d,v_min,v_max(1),'fast');
    % detect_L(:,2)=a;first_time_L(:,2)=b;first_encounter_time_L(:,2)=c;effective_pulses_L(:,2)=d;
    % [a,b,c,d]=MC_func(R,N,f,'line',fasan_D,fasan_d,v_min,v_max(2),'slow');
    % detect_L(:,3)=a;first_time_L(:,3)=b;first_encounter_time_L(:,3)=c;effective_pulses_L(:,3)=d;
    % [a,b,c,d]=MC_func(R,N,f,'line',fasan_D,fasan_d,v_min,v_max(2),'fast');
    % detect_L(:,4)=a;first_time_L(:,4)=b;first_encounter_time_L(:,4)=c;effective_pulses_L(:,4)=d;
    % save([path,'MC_1par_EA_LINE_D',num2str(fasan_D*1000),'d',num2str(fasan_d*1000),'mrad.mat'],"detect_L","first_time_L","first_encounter_time_L","effective_pulses_L");

    % Point Gaussian beam:
    % fasan_D is the exp(-2) full divergence angle. fasan_d is ignored.
    % Receiver FOV full angle is fixed inside MC_point_snr as 3*fasan_D.
    % [a,b,c,d]=MC_func(R,N,f,'point',fasan_D,[],v_min,30,'fast');
    % detect_P(:,1)=a;first_time_P(:,1)=b;first_encounter_time_P(:,1)=c;effective_pulses_P(:,1)=d;
    % save([path,'MC_1par_EA_POINT_D',num2str(fasan_D*1000),'mrad.mat'],"detect_P","first_time_P","first_encounter_time_P","effective_pulses_P");




end

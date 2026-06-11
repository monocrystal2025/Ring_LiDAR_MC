clc;clear;close all;
%=========参数与初始化===========%
N=50000;
f=5e3;
z=500;
JIAODU=25;%以度计量的半角
fasan_point=15e-3;
fasan_D_LIST=5e-3:5e-3:200e-3;
v_max=30;
v_min=1;

detect_R_RA=zeros(N,4);detect_L_RA=zeros(N,4);detect_R_SP=zeros(N,4);detect_L_SP=zeros(N,4);
first_time_R_RA=zeros(N,4);first_time_L_RA=zeros(N,4);first_time_R_SP=zeros(N,4);first_time_L_SP=zeros(N,4);
% [a,b]=MC_func(R,N,f,'ring',150e-3,1.5e-3,v_min,v_max(1),'fast');
path='D:\matlab\MC_ROI_RESULTS_SNR\';

for i=length(fasan_D_LIST):-1:5
    fasan_D=fasan_D_LIST(i);
    fasan_d=0.8e-3;
    % fasan_d=fasan_D/2;
    [a,b]=MC_func_ROI(z,N,f,JIAODU,'ring',fasan_D,fasan_d,v_min,v_max,'fast','raster');
    detect_R_RA(:,1)=a;first_time_R_RA(:,1)=b;
    % [a,b]=MC_func_ROI(z,N,f,JIAODU,'ring',fasan_D,fasan_d,v_min,v_max(1),'fast','raster');
    detect_R_RA(:,2)=a;first_time_R_RA(:,2)=b;
    % [a,b]=MC_func_ROI(z,N,f,JIAODU,'ring',fasan_D,fasan_d,v_min,v_max(2),'slow','raster');
    detect_R_RA(:,3)=a;first_time_R_RA(:,3)=b;
    % [a,b]=MC_func_ROI(z,N,f,JIAODU,'ring',fasan_D,fasan_d,v_min,v_max(2),'fast','raster');
    detect_R_RA(:,4)=a;first_time_R_RA(:,4)=b;    
    save([path,'MC_1par_ROI_RING_RA_D',num2str(fasan_D*1000),'d',num2str(fasan_d*1000),'mrad.mat'],"detect_R_RA","first_time_R_RA");

    [a,b]=MC_func_ROI(z,N,f,JIAODU,'ring',fasan_D,fasan_d,v_min,v_max,'fast','spiral');
    detect_R_SP(:,1)=a;first_time_R_SP(:,1)=b;
    % [a,b]=MC_func_ROI(z,N,f,JIAODU,'ring',fasan_D,fasan_d,v_min,v_max(1),'fast','spiral');
    detect_R_SP(:,2)=a;first_time_R_SP(:,2)=b;
    % [a,b]=MC_func_ROI(z,N,f,JIAODU,'ring',fasan_D,fasan_d,v_min,v_max(2),'slow','spiral');
    detect_R_SP(:,3)=a;first_time_R_SP(:,3)=b;
    % [a,b]=MC_func_ROI(z,N,f,JIAODU,'ring',fasan_D,fasan_d,v_min,v_max(2),'fast','spiral');
    detect_R_SP(:,4)=a;first_time_R_SP(:,4)=b;
    save([path,'MC_1par_ROI_RING_SP_D',num2str(fasan_D*1000),'d',num2str(fasan_d*1000),'mrad.mat'],"detect_R_SP","first_time_R_SP");
    % 
    % [a,b]=MC_func_ROI(z,N,f,JIAODU,'line',fasan_D,fasan_d,v_min,v_max,'fast','raster');
    % detect_L_RA(:,1)=a;first_time_L_RA(:,1)=b;
    % % [a,b]=MC_func_ROI(z,N,f,JIAODU,'line',fasan_D,fasan_d,v_min,v_max(1),'fast','raster');
    % detect_L_RA(:,2)=a;first_time_L_RA(:,2)=b;
    % % [a,b]=MC_func_ROI(z,N,f,JIAODU,'line',fasan_D,fasan_d,v_min,v_max(2),'slow','raster');
    % detect_L_RA(:,3)=a;first_time_L_RA(:,3)=b;
    % % [a,b]=MC_func_ROI(z,N,f,JIAODU,'line',fasan_D,fasan_d,v_min,v_max(2),'fast','raster');
    % detect_L_RA(:,4)=a;first_time_L_RA(:,4)=b;
    % save([path,'MC_1par_ROI_LINE_RA_D',num2str(fasan_D*1000),'d',num2str(fasan_d*1000),'mrad.mat'],"detect_L_RA","first_time_L_RA");

    % [a,b]=MC_func_ROI(z,N,f,JIAODU,'line',fasan_D,fasan_d,v_min,v_max,'fast','spiral');
    % detect_L_SP(:,1)=a;first_time_L_SP(:,1)=b;
    % % [a,b]=MC_func_ROI(z,N,f,JIAODU,'line',fasan_D,fasan_d,v_min,v_max(1),'fast','raster');
    % detect_L_SP(:,2)=a;first_time_L_SP(:,2)=b;
    % % [a,b]=MC_func_ROI(z,N,f,JIAODU,'line',fasan_D,fasan_d,v_min,v_max(2),'slow','raster');
    % detect_L_SP(:,3)=a;first_time_L_SP(:,3)=b;
    % % [a,b]=MC_func_ROI(z,N,f,JIAODU,'line',fasan_D,fasan_d,v_min,v_max(2),'fast','raster');
    % detect_L_SP(:,4)=a;first_time_L_SP(:,4)=b;
    % save([path,'MC_1par_ROI_LINE_SP_D',num2str(fasan_D*1000),'d',num2str(fasan_d*1000),'mrad.mat'],"detect_L_SP","first_time_L_SP");


end
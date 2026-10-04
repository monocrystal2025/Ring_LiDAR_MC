clc;clear;close all;
%=========参数与初始化===========%
N=50000;
f=5e3;
fasan_D_LIST=5e-3:5e-3:400e-3;
v_min=1;
R=2000;

save_dir='D:\lzx\MC_RESULTS';
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

for i=length(fasan_D_LIST):-1:1
    fasan_D=fasan_D_LIST(i);
    % Use the measured piecewise optimum only for pulse-vector block size.
    ring_block_size=mc_ea_block_size(fasan_D,'ring');
    line_block_size=mc_ea_block_size(fasan_D,'line');
    point_block_size=mc_ea_block_size(fasan_D,'point');
    fasan_d=1e-3;

    % fprintf('Starting ring beam, D=%g mrad.\n',fasan_D*1000);
    % [detect_R,first_time_R,first_encounter_time_R,effective_pulses_R]= ...
    %     MC_func(R,N,f,'ring',fasan_D,fasan_d, ...
    %     v_min,30,'fast',ring_block_size,'fixed');
    % detect_R=logical(detect_R);
    % output_file=fullfile(save_dir, ...
    %     ['MC_1par_EA_RING_D',num2str(fasan_D*1000), ...
    %     'd',num2str(fasan_d*1000),'mrad.mat']);
    % save(output_file,'detect_R','first_time_R', ...
    %     'first_encounter_time_R','effective_pulses_R');
    % fprintf('Saved %s\n',output_file);
    % clear detect_R first_time_R first_encounter_time_R effective_pulses_R

    fprintf('Starting line beam, D=%g mrad.\n',fasan_D*1000);
    [detect_L,first_time_L,first_encounter_time_L,effective_pulses_L]= ...
        MC_func(R,N,f,'line',fasan_D,fasan_d, ...
        v_min,30,'fast',line_block_size,'fixed');
    detect_L=logical(detect_L);
    output_file=fullfile(save_dir, ...
        ['MC_1par_EA_LINE_D',num2str(fasan_D*1000), ...
        'd',num2str(fasan_d*1000),'mrad.mat']);
    save(output_file,'detect_L','first_time_L', ...
        'first_encounter_time_L','effective_pulses_L');
    fprintf('Saved %s\n',output_file);
    clear detect_L first_time_L first_encounter_time_L effective_pulses_L

    % fprintf('Starting point beam, D=%g mrad.\n',fasan_D*1000);
    % [detect_P,first_time_P,first_encounter_time_P,effective_pulses_P]= ...
    %     MC_func(R,N,f,'point',fasan_D,[], ...
    %     v_min,30,'fast',point_block_size,'fixed');
    % detect_P=logical(detect_P);
    % output_file=fullfile(save_dir, ...
    %     ['MC_1par_EA_POINT_D',num2str(fasan_D*1000),'mrad.mat']);
    % save(output_file,'detect_P','first_time_P', ...
    %     'first_encounter_time_P','effective_pulses_P');
    % fprintf('Saved %s\n',output_file);
    % clear detect_P first_time_P first_encounter_time_P effective_pulses_P
end

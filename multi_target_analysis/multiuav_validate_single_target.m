function validation = multiuav_validate_single_target(N, seed)
%MULTIUAV_VALIDATE_SINGLE_TARGET Compare M=1 with current single-target code.

arguments
    N (1,1) double {mustBeInteger, mustBePositive} = 300
    seed (1,1) double {mustBeInteger, mustBeNonnegative} = 20260819
end

moduleDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(moduleDir);
addpath(moduleDir, projectRoot);
cfg = multiuav_default_config();
cfg.runN = N;
cfg.system.maxObservation_s = 70;
cfg.scene.entrySpan_s = 0;
cfg.scene.turnRateRms_deg_s = 0;
cfg.scene.turnRateLimit_deg_s = 0;
bank = multiuav_build_scenario_bank(N, 1, seed);

Trajectory = ["RA"; "SP"];
Beam = ["line"; "ring"];
wD_mrad = [125; 125];
PathType = ["raster"; "spiral"];
NewP = zeros(2, 1);
CurrentSingleP = zeros(2, 1);
NewMeanTime_s = zeros(2, 1);
CurrentSingleMeanTime_s = zeros(2, 1);

for index = 1:2
    beam = struct('name', Beam(index), 'label', Beam(index), ...
        'wD_rad', wD_mrad(index) * 1e-3, 'wd_rad', 0.8e-3);
    result = multiuav_simulate_case(cfg, beam, "ROI", ...
        Trajectory(index), "sync_dispersed", bank);
    NewP(index) = mean(result.detect);
    NewMeanTime_s(index) = mean(result.firstTime_s, 'omitnan');

    rng(seed, 'twister');
    [detected, firstTime] = MC_func_ROI_NEW(2000, N, 5000, 30, ...
        char(Beam(index)), beam.wD_rad, beam.wd_rad, 30, 30, ...
        'fast', char(PathType(index)));
    CurrentSingleP(index) = mean(detected);
    CurrentSingleMeanTime_s(index) = mean(firstTime, 'omitnan');
end

validation = table(Trajectory, Beam, wD_mrad, NewP, CurrentSingleP, ...
    100 .* (NewP - CurrentSingleP), NewMeanTime_s, ...
    CurrentSingleMeanTime_s, NewMeanTime_s - CurrentSingleMeanTime_s, ...
    'VariableNames', {'Trajectory', 'Beam', 'wD_mrad', 'NewP', ...
    'CurrentSingleP', 'DeltaP_pp', 'NewMeanTime_s', ...
    'CurrentSingleMeanTime_s', 'DeltaMeanTime_s'});
end

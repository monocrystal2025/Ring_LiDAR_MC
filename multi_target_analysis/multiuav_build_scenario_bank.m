function bank = multiuav_build_scenario_bank(N, maxTargets, seed)
%MULTIUAV_BUILD_SCENARIO_BANK Generate paired multi-target random draws.

arguments
    N (1,1) double {mustBeInteger, mustBePositive}
    maxTargets (1,1) double {mustBeInteger, mustBePositive}
    seed (1,1) double {mustBeInteger, mustBeNonnegative}
end

stream = RandStream('mt19937ar', 'Seed', seed);
drawSize = [N, maxTargets];

bank = struct();
bank.N = N;
bank.maxTargets = maxTargets;
bank.seed = seed;
bank.scenarioId = (1:N).';
bank.targetId = 1:maxTargets;
bank.entryAzimuth_u = rand(stream, drawSize);
bank.entryRange_u = rand(stream, drawSize);
bank.velocityAzimuth_u = rand(stream, drawSize);
bank.velocityCosPolar_u = rand(stream, drawSize);
bank.entryTime_u = rand(stream, drawSize);
bank.scanPhase_u = rand(stream, N, 1);
bank.formationRotation_u = rand(stream, N, 1);
end

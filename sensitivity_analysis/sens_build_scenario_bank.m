function bank = sens_build_scenario_bank(N, seed)
%SENS_BUILD_SCENARIO_BANK Generate normalized paired Monte Carlo scenarios.
%
% The bank is independent of range, ROI angle, target speed, beam shape,
% trajectory, and path length. sens_simulate_roi_case maps these normalized
% draws to each physical case, ensuring paired comparisons.

arguments
    N (1,1) double {mustBeInteger, mustBePositive}
    seed (1,1) double {mustBeInteger, mustBeNonnegative}
end

stream = RandStream('mt19937ar', 'Seed', seed);

bank = struct();
bank.N = N;
bank.seed = seed;
bank.entryAzimuth_u = rand(stream, N, 1);
bank.entryRange_u = rand(stream, N, 1);
bank.velocityAzimuth_u = rand(stream, N, 1);
bank.velocityCosPolar_u = rand(stream, N, 1);
bank.scanPhase_u = rand(stream, N, 1);
bank.scenarioId = (1:N).';
end

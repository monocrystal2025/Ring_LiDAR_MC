function bank = ea_sens_build_scenario_bank(N, seed)
%EA_SENS_BUILD_SCENARIO_BANK Generate paired ROE Monte Carlo scenarios.
%
% Normalized random draws are mapped to each physical case by the simulator.
% Reusing one bank across beams, divergence angles, and factor levels applies
% common random numbers and substantially reduces comparison noise.

arguments
    N (1,1) double {mustBeInteger, mustBePositive}
    seed (1,1) double {mustBeInteger, mustBeNonnegative}
end

stream = RandStream("mt19937ar", "Seed", seed);

bank = struct();
bank.N = N;
bank.seed = seed;
bank.entryAzimuth_u = rand(stream, N, 1);
bank.entryCosColatitude_u = rand(stream, N, 1);
bank.velocityAzimuth_u = rand(stream, N, 1);
bank.velocityCosPolar_u = rand(stream, N, 1);
bank.scanPhase_u = rand(stream, N, 1);
bank.scenarioId = (1:N).';
end

classdef TestMultiUavAnalysis < matlab.unittest.TestCase
    %TESTMULTIUAVANALYSIS Unit tests for public multi-UAV interfaces.

    methods (TestClassSetup)
        function addModulePaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            moduleDir = fileparts(fileparts(mfilename('fullpath')));
            projectRoot = fileparts(moduleDir);
            testCase.applyFixture(PathFixture(moduleDir));
            testCase.applyFixture(PathFixture(projectRoot));
        end
    end

    methods (Test)
        function scenarioBankIsReproducible(testCase)
            first = multiuav_build_scenario_bank(8, 4, 1234);
            second = multiuav_build_scenario_bank(8, 4, 1234);

            testCase.verifyEqual(first, second);
            testCase.verifySize(first.scanPhase_u, [8, 1]);
            testCase.verifySize(first.entryAzimuth_u, [8, 4]);
        end

        function rangeCellCollisionMarksEveryConflictingTarget(testCase)
            ranges = [100, 104, 120; 200, 204, 220];
            candidates = logical([1, 1, 0; 1, 0, 1]);

            [resolved, ambiguous, rangeBin] = ...
                multiuav_apply_range_resolution(ranges, candidates, 7.5);

            testCase.verifyEqual(ambiguous, logical([1, 1, 0; 0, 0, 0]));
            testCase.verifyEqual(resolved, logical([0, 0, 0; 1, 0, 1]));
            testCase.verifyEqual(rangeBin(1, 1), rangeBin(1, 2));
        end

        function scanPhaseIsSharedAndSimulationIsReproducible(testCase)
            cfg = TestMultiUavAnalysis.testConfig();
            bank = multiuav_build_scenario_bank(5, 3, 9001);
            beam = TestMultiUavAnalysis.testBeam("ring", 125);

            first = multiuav_simulate_case(cfg, beam, "ROI", "RA", ...
                "sync_dispersed", bank);
            second = multiuav_simulate_case(cfg, beam, "ROI", "RA", ...
                "sync_dispersed", bank);

            testCase.verifyEqual(first.initialScanPhase_u, bank.scanPhase_u);
            testCase.verifySize(first.initialScanPhase_u, [5, 1]);
            testCase.verifyEqual(first.detect, second.detect);
            testCase.verifyEqual(first.firstTime_s, second.firstTime_s);
            testCase.verifyLessThanOrEqual( ...
                max(first.maxInterceptedEnergyFraction), 1 + 1e-12);
        end

        function zeroTurnRateCollapsesToStraightLine(testCase)
            cfg = TestMultiUavAnalysis.testConfig();
            cfg.scene.turnRateRms_deg_s = 0;
            cfg.scene.turnRateLimit_deg_s = 0;
            cfg.scene.entrySpan_s = 0;
            cfg.system.snrThreshold = Inf;
            bank = multiuav_build_scenario_bank(4, 1, 7001);
            beam = TestMultiUavAnalysis.testBeam("point", 125);

            maneuver = multiuav_simulate_case(cfg, beam, "ROI", "SP", ...
                "async_maneuver", bank);
            straight = multiuav_simulate_case(cfg, beam, "ROI", "SP", ...
                "sync_dispersed", bank);

            testCase.verifyEqual(maneuver.finalPosition_m, ...
                straight.finalPosition_m, 'AbsTol', 1e-12);
            testCase.verifyEqual(maneuver.finalVelocity_m_s, ...
                straight.finalVelocity_m_s, 'AbsTol', 1e-12);
        end

        function summaryContainsRequiredPaperMetrics(testCase)
            cfg = TestMultiUavAnalysis.testConfig();
            cfg.scene.timelyThresholds_s = 0.1;
            bank = multiuav_build_scenario_bank(4, 2, 8001);
            beam = TestMultiUavAnalysis.testBeam("line", 125);
            result = multiuav_simulate_case(cfg, beam, "ROI", "RA", ...
                "sync_formation", bank);

            [summary, comparisons] = multiuav_summarize({result}, cfg);

            testCase.verifyEqual(height(summary), 1);
            testCase.verifyTrue(all(ismember({'PAll', 'MeanRecall', ...
                'MeanCompleteTime_s', 'MeanRepeatedEncounters', ...
                'MeanMaxTargetsPerPulse', 'RangeAmbiguityRate'}, ...
                summary.Properties.VariableNames)));
            testCase.verifyEmpty(comparisons);
        end

        function formationOptimizationFindsIndependentBeamOptima(testCase)
            cfg = TestMultiUavAnalysis.testConfig();
            cfg.optimization.threshold_s = 0.5;
            cfg.optimization.cvFolds = 5;
            results = TestMultiUavAnalysis.syntheticOptimizationResults();

            [optima, comparisons] = ...
                multiuav_analyze_formation_optima(results, cfg);

            line = optima.Beam == "line";
            ring = optima.Beam == "ring";
            testCase.verifyEqual(optima.BestWD_mrad(line), 200);
            testCase.verifyEqual(optima.BestWD_mrad(ring), 100);
            testCase.verifyEqual(comparisons.Winner95, "Ring");
            testCase.verifyGreaterThan(comparisons.DeltaCVPAllCILow_pp, 0);
            testCase.verifyEqual(cfg.manuscriptWD_mrad, 5:5:500);
        end

        function narrowBeamFormationConservesPulseEnergy(testCase)
            cfg = TestMultiUavAnalysis.testConfig();
            cfg.scene.formationSpacing_m = 5;
            bank = multiuav_build_scenario_bank(5, 8, 9101);
            beam = TestMultiUavAnalysis.testBeam("ring", 5);

            result = multiuav_simulate_case(cfg, beam, "ROI", "SP", ...
                "sync_formation", bank);

            testCase.verifyLessThanOrEqual( ...
                max(result.maxInterceptedEnergyFraction), 1 + 1e-12);
        end
    end

    methods (Static, Access = private)
        function cfg = testConfig()
            cfg = multiuav_default_config();
            cfg.system.maxRange_m = 220;
            cfg.system.minimumEntryRange_m = 100;
            cfg.system.prf_Hz = 50;
            cfg.system.scanRate_rad_s = 2 * pi;
            cfg.system.maxObservation_s = 0.1;
            cfg.system.noiseLutStep_m = 2;
            cfg.system.targetSampleGrid = 3;
            cfg.execution.motionUpdate_s = 0.02;
            cfg.scene.entrySpan_s = 0;
            cfg.scene.timelyThresholds_s = 0.1;
        end

        function beam = testBeam(name, width_mrad)
            beam = struct('name', name, 'label', name, ...
                'wD_rad', width_mrad * 1e-3, 'wd_rad', 1e-3);
            if name == "point"
                beam.wd_rad = NaN;
            end
        end

        function results = syntheticOptimizationResults()
            N = 20;
            M = 2;
            definitions = {
                "line", 100, false(N, 1);
                "line", 200, [true(5, 1); false(15, 1)];
                "ring", 100, true(N, 1);
                "ring", 200, [true(10, 1); false(10, 1)]};
            results = cell(1, size(definitions, 1));
            for index = 1:size(definitions, 1)
                captured = definitions{index, 3};
                detect = repmat(captured, 1, M);
                firstTime = nan(N, M);
                firstTime(detect) = 0.1;
                results{index} = struct( ...
                    'region', "ROI", 'trajectory', "RA", ...
                    'scenario', "sync_formation", ...
                    'targetCount', M, 'beam', definitions{index, 1}, ...
                    'wD_mrad', definitions{index, 2}, ...
                    'scenarioId', (1:N).', 'detect', detect, ...
                    'firstTime_s', firstTime, ...
                    'ambiguousHitCount', zeros(N, M), ...
                    'resolvedHitCount', ones(N, M));
            end
        end
    end
end

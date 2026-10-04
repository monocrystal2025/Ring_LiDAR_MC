classdef eaSensitivityTest < matlab.unittest.TestCase
    %EASENSITIVITYTEST Unit tests for the standalone ROE sensitivity module.

    methods (TestClassSetup)
        function addModuleToPath(testCase)
            moduleDir = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture( ...
                matlab.unittest.fixtures.PathFixture(moduleDir));
        end
    end

    methods (Test)
        function testDefaultBatchSizeAndFactors(testCase)
            cfg = ea_sens_default_config();

            testCase.verifyEqual(cfg.screenN, 3000);
            testCase.verifyEqual(cfg.refineN, 3000);
            testCase.verifyEqual(numel(cfg.factors), 9);
            testCase.verifyEqual(cfg.system.beamWidth_rad, 0.8e-3, ...
                AbsTol=1e-15);
            testCase.verifyLessThanOrEqual(max([cfg.screenN, cfg.refineN]), ...
                3000);
        end

        function testScenarioBankIsReproducible(testCase)
            first = ea_sens_build_scenario_bank(12, 42);
            second = ea_sens_build_scenario_bank(12, 42);

            testCase.verifyEqual(first.entryAzimuth_u, ...
                second.entryAzimuth_u, AbsTol=0);
            testCase.verifyEqual(first.velocityCosPolar_u, ...
                second.velocityCosPolar_u, AbsTol=0);
            testCase.verifyEqual(first.scanPhase_u, ...
                second.scanPhase_u, AbsTol=0);
            testCase.verifyGreaterThanOrEqual(first.entryCosColatitude_u, 0);
            testCase.verifyLessThanOrEqual(first.entryCosColatitude_u, 1);
        end

        function testSyntheticIsoPerformance(testCase)
            cfg = ea_sens_default_config();
            cfg.timelyThresholds_s = 10;
            cfg.targetProbabilities = 0.5;
            cfg.primaryTimelyThreshold_s = 10;
            cfg.primaryProbability = 0.5;
            cfg.bootstrapN = 0;
            raw = eaSensitivityTest.syntheticRaw();

            [summary, timely, iso, effects] = ...
                ea_sens_analyze_results(raw, cfg);

            testCase.verifyEqual(height(summary), 4);
            testCase.verifyEqual(height(timely), 4);
            testCase.verifyEqual(height(iso), 1);
            testCase.verifyEqual(iso.RingMinWD_mrad, 100, AbsTol=1e-10);
            testCase.verifyEqual(iso.LineMinWD_mrad, 150, AbsTol=1e-10);
            testCase.verifyEqual(iso.KappaRingLine, 2/3, AbsTol=1e-10);
            testCase.verifyTrue(iso.AnnularAdvantage);
            testCase.verifyEmpty(effects);
        end

        function testSmallParameterizedSimulation(testCase)
            cfg = ea_sens_default_config();
            cfg.execution.useParallel = false;
            cfg.execution.blockSize = 100;
            cfg.system.maxRange_m = 50;
            cfg.system.prf_Hz = 50;
            cfg.system.scanRate_rad_s = 2 * pi;
            cfg.system.targetSpeed_m_s = 20;
            cfg.system.noiseLutStep_m = 2;
            bank = ea_sens_build_scenario_bank(3, 7);
            beam = struct("name", "ring", "label", "Annular", ...
                "wD_rad", 0.2, "wd_rad", 0.8e-3);

            result = ea_sens_simulate_case(cfg, beam, bank);

            testCase.verifySize(result.detect, [3, 1]);
            testCase.verifyEqual(result.scenarioId, bank.scenarioId);
            testCase.verifyGreaterThan(result.cycleLength, 0);
            testCase.verifyGreaterThan(result.windowPulses, 0);
            testCase.verifyTrue(all( ...
                result.firstTime_s(isfinite(result.firstTime_s)) >= 0));
        end

        function testUnknownBeamErrors(testCase)
            cfg = ea_sens_default_config();
            bank = ea_sens_build_scenario_bank(1, 1);
            beam = struct("name", "invalid", "label", "Invalid", ...
                "wD_rad", 0.1, "wd_rad", 0.8e-3);

            testCase.verifyError( ...
                @() ea_sens_simulate_case(cfg, beam, bank), ...
                "ea_sens_simulate_case:UnknownBeam");
        end
    end

    methods (Test, TestTags={'Integration'})
        function testSmokePipeline(testCase)
            outputRoot = string(tempname);
            mkdir(outputRoot);
            testCase.addTeardown(@() rmdir(outputRoot, "s"));

            outputs = run_ea_sensitivity_analysis( ...
                "Mode", "smoke", ...
                "OutputRoot", outputRoot, ...
                "UseParallel", false, ...
                "Force", true, ...
                "BootstrapN", 0);

            testCase.verifyEqual(outputs.config.screenN, 6);
            testCase.verifyEqual(outputs.config.refineN, 6);
            testCase.verifyNotEmpty(outputs.screenSummary);
            testCase.verifyNotEmpty(outputs.refineSummary);
            testCase.verifyNotEmpty(outputs.interactionSummary);
            testCase.verifyTrue(isfile(fullfile(outputRoot, "results", ...
                "ea_sensitivity_screen_summary.csv")));
            testCase.verifyTrue(isfile(fullfile(outputRoot, "results", ...
                "ea_sensitivity_refined_iso_performance.csv")));
            testCase.verifyTrue(isfile(fullfile(outputRoot, "results", ...
                "ea_sensitivity_interaction_iso_performance.csv")));
            testCase.verifyTrue(isfile(fullfile(outputRoot, "figures", ...
                "ea_sensitivity_main.png")));
        end
    end

    methods (Static, Access=private)
        function raw = syntheticRaw()
            scenarioId = (1:10).';
            raw = repmat(eaSensitivityTest.rawTemplate(), 1, 4);
            raw(1) = eaSensitivityTest.makeRaw( ...
                "ring", 0.1, 6, scenarioId);
            raw(2) = eaSensitivityTest.makeRaw( ...
                "ring", 0.2, 8, scenarioId);
            raw(3) = eaSensitivityTest.makeRaw( ...
                "line", 0.1, 2, scenarioId);
            raw(4) = eaSensitivityTest.makeRaw( ...
                "line", 0.2, 8, scenarioId);
        end

        function item = makeRaw(beam, wD, timelyCount, scenarioId)
            item = eaSensitivityTest.rawTemplate();
            item.beam = beam;
            item.wD_rad = wD;
            item.scenarioId = scenarioId;
            item.detect = false(numel(scenarioId), 1);
            item.detect(1:timelyCount) = true;
            item.firstTime_s = nan(numel(scenarioId), 1);
            item.firstTime_s(1:timelyCount) = 1;
        end

        function item = rawTemplate()
            item = struct( ...
                "caseId", "baseline", ...
                "factor", "baseline", ...
                "factorLabel", "Baseline", ...
                "level", "Baseline", ...
                "levelValue", 0, ...
                "unit", "", ...
                "beam", "", ...
                "wD_rad", nan, ...
                "wd_rad", 0.8e-3, ...
                "scenarioId", zeros(0, 1), ...
                "detect", false(0, 1), ...
                "firstTime_s", zeros(0, 1));
        end
    end
end

classdef test_sensitivity_analysis < matlab.unittest.TestCase
    %TEST_SENSITIVITY_ANALYSIS Public-interface tests for the new module.

    methods (TestClassSetup)
        function addModuleToPath(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            moduleFolder = fileparts(testFolder);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(moduleFolder));
        end
    end

    methods (Test)
        function defaultConfigurationMatchesStudy(testCase)
            cfg = sens_default_config();

            testCase.verifyEqual(cfg.system.maxRange_m, 2000, AbsTol=0);
            testCase.verifyEqual(cfg.system.prf_Hz, 5000, AbsTol=0);
            testCase.verifyEqual(cfg.system.scanRate_rad_s, 2*pi, AbsTol=1e-14);
            testCase.verifyEqual([cfg.beams.wD_rad], [20e-3, 180e-3, 150e-3], ...
                AbsTol=1e-14);
            testCase.verifyEqual(numel(cfg.factors), 11);
        end

        function scenarioBankIsDeterministic(testCase)
            first = sens_build_scenario_bank(16, 42);
            second = sens_build_scenario_bank(16, 42);

            testCase.verifyEqual(first.entryAzimuth_u, second.entryAzimuth_u, AbsTol=0);
            testCase.verifyEqual(first.velocityCosPolar_u, second.velocityCosPolar_u, AbsTol=0);
            testCase.verifyEqual(first.scanPhase_u, second.scanPhase_u, AbsTol=0);
            testCase.verifyEqual(first.scenarioId, second.scenarioId);
        end

        function scenarioSeedsProduceDifferentBanks(testCase)
            first = sens_build_scenario_bank(16, 42);
            second = sens_build_scenario_bank(16, 43);

            testCase.verifyNotEqual(first.entryRange_u, second.entryRange_u);
        end

        function simulationIsReproducible(testCase)
            cfg = test_sensitivity_analysis.smallConfig();
            bank = sens_build_scenario_bank(4, 101);
            beam = cfg.beams(1);

            first = sens_simulate_roi_case(cfg, beam, "RA", bank);
            second = sens_simulate_roi_case(cfg, beam, "RA", bank);

            testCase.verifyEqual(first.detect, second.detect);
            testCase.verifyEqual(first.firstTime_s, second.firstTime_s, AbsTol=0);
            testCase.verifyEqual(first.scenarioId, bank.scenarioId);
            testCase.verifySize(first.detect, [4, 1]);
        end

        function simulationRejectsUnknownBeam(testCase)
            cfg = test_sensitivity_analysis.smallConfig();
            bank = sens_build_scenario_bank(2, 4);
            beam = cfg.beams(1);
            beam.name = "unknown";

            testCase.verifyError( ...
                @() sens_simulate_roi_case(cfg, beam, "RA", bank), ...
                'sens_simulate_roi_case:UnknownBeam');
        end

        function analysisComputesPairedEffects(testCase)
            cfg = sens_default_config();
            cfg.bootstrapN = 20;
            raw = test_sensitivity_analysis.syntheticRawResults();

            [summary, comparisons, effects] = sens_analyze_results(raw, cfg);

            testCase.verifyEqual(height(summary), 6);
            testCase.verifyEqual(height(comparisons), 2);
            testCase.verifyEqual(height(effects), 1);
            testCase.verifyEqual(effects.ProbabilityEffect_pp, 25, AbsTol=1e-12);
            testCase.verifyLessThanOrEqual( ...
                effects.ProbabilityEffect_CI_Low_pp, effects.ProbabilityEffect_pp);
            testCase.verifyGreaterThanOrEqual( ...
                effects.ProbabilityEffect_CI_High_pp, effects.ProbabilityEffect_pp);
            testCase.verifyEqual(summary.PDetect(summary.Beam == "ring" & ...
                summary.Level == "High"), 0.75, AbsTol=1e-12);
        end
    end

    methods (Test, TestTags = {'Integration'})
        function smokePipelineWritesOutputs(testCase)
            outputFolder = string(tempname);
            mkdir(outputFolder);
            testCase.addTeardown(@() rmdir(outputFolder, 's'));

            outputs = run_sensitivity_analysis('Mode', 'smoke', ...
                'OutputRoot', outputFolder, 'ScreenN', 4, 'RefineN', 5, ...
                'BootstrapN', 5, 'UseParallel', false, 'Force', true);

            testCase.verifyTrue(isfile(fullfile(outputFolder, 'results', ...
                'sensitivity_screen_summary.csv')));
            testCase.verifyTrue(isfile(fullfile(outputFolder, 'results', ...
                'sensitivity_refined_summary.csv')));
            testCase.verifyTrue(isfile(fullfile(outputFolder, 'figures', ...
                'sensitivity_analysis_main.png')));
            testCase.verifyTrue(isfile(fullfile(outputFolder, 'results', ...
                'effect_direction_validation.csv')));
            testCase.verifyEqual(height(outputs.selectedFactors), 1);
        end
    end

    methods (Static, Access=private)
        function cfg = smallConfig()
            cfg = sens_default_config();
            cfg.execution.useParallel = false;
            cfg.execution.blockSize = 500;
            cfg.system.maxRange_m = 150;
            cfg.system.minimumEntryRange_m = 100;
            cfg.system.prf_Hz = 100;
            cfg.system.noiseLutStep_m = 2;
        end

        function raw = syntheticRawResults()
            ids = (1:4).';
            template = struct('caseId', "", 'factor', "factorX", ...
                'factorLabel', "Synthetic factor", 'level', "", ...
                'levelValue', nan, 'unit', "1", 'beam', "", ...
                'trajectory', "RA", 'scenarioId', ids, 'detect', false(4,1), ...
                'firstTime_s', nan(4,1));
            raw = repmat(template, 1, 6);

            raw(1) = test_sensitivity_analysis.makeSynthetic(template, ...
                "x_low", "Low", 0, "point", [1;0;0;0], [3;nan;nan;nan]);
            raw(2) = test_sensitivity_analysis.makeSynthetic(template, ...
                "x_low", "Low", 0, "line", [1;0;0;0], [2;nan;nan;nan]);
            raw(3) = test_sensitivity_analysis.makeSynthetic(template, ...
                "x_low", "Low", 0, "ring", [1;1;0;0], [1;2;nan;nan]);
            raw(4) = test_sensitivity_analysis.makeSynthetic(template, ...
                "x_high", "High", 1, "point", [1;0;0;0], [3;nan;nan;nan]);
            raw(5) = test_sensitivity_analysis.makeSynthetic(template, ...
                "x_high", "High", 1, "line", [1;0;0;0], [2;nan;nan;nan]);
            raw(6) = test_sensitivity_analysis.makeSynthetic(template, ...
                "x_high", "High", 1, "ring", [1;1;1;0], [1;2;3;nan]);
        end

        function item = makeSynthetic(template, caseId, level, value, beam, detect, time)
            item = template;
            item.caseId = caseId;
            item.level = level;
            item.levelValue = value;
            item.beam = beam;
            item.detect = logical(detect);
            item.firstTime_s = time;
        end
    end
end

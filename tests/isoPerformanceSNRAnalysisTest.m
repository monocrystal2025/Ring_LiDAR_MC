classdef isoPerformanceSNRAnalysisTest < matlab.unittest.TestCase
    %ISOPERFORMANCESNRANALYSISTEST Tests the compact paper analysis.

    properties (Constant)
        ProjectFolder = 'D:\lzx\MatlabCode\codexagent_MC'
        ResultFolder = 'D:\lzx\MC_RESULTS'
    end

    methods (TestClassSetup)
        function addProjectFolder(testCase)
            import matlab.unittest.fixtures.PathFixture
            testCase.applyFixture(PathFixture(testCase.ProjectFolder));
        end
    end

    methods (Test)
        function testEqualPerformanceMonteCarloValues(testCase)
            result = iso_performance_snr_analysis('', ...
                testCase.ResultFolder, false);

            testCase.verifyEqual( ...
                result.mc.annular.detectionProbability, 0.9151, ...
                'AbsTol', 1e-12);
            testCase.verifyEqual( ...
                result.mc.line.detectionProbability, 0.91296, ...
                'AbsTol', 1e-12);
            testCase.verifyLessThan(abs( ...
                result.mc.annular.detectionProbability - ...
                result.mc.line.detectionProbability), 0.005);
        end

        function testEffectivePulseStatistics(testCase)
            result = iso_performance_snr_analysis('', ...
                testCase.ResultFolder, false);

            testCase.verifyEqual( ...
                result.mc.annular.meanEffectivePulses, ...
                1.73827996940225, 'AbsTol', 1e-12);
            testCase.verifyEqual( ...
                result.mc.line.meanEffectivePulses, ...
                1.00383368384157, 'AbsTol', 1e-12);
            testCase.verifyGreaterThan( ...
                result.metrics.probabilityAtLeastTwoPulsesAnnular, 0.26);
            testCase.verifyLessThan( ...
                result.metrics.probabilityAtLeastTwoPulsesLine, 0.003);
        end

        function testPhotonBudgetIsPhysical(testCase)
            result = iso_performance_snr_analysis('', ...
                testCase.ResultFolder, false);

            testCase.verifyGreaterThan(result.annular.signalPhotons, 0);
            testCase.verifyGreaterThan(result.line.signalPhotons, 0);
            testCase.verifyGreaterThan(result.annular.noisePhotons, 0);
            testCase.verifyGreaterThan(result.line.noisePhotons, 0);
            testCase.verifyGreaterThanOrEqual( ...
                result.annular.illuminationFactor, 0);
            testCase.verifyLessThanOrEqual( ...
                result.annular.illuminationFactor, 1);
            testCase.verifyGreaterThanOrEqual( ...
                result.line.illuminationFactor, 0);
            testCase.verifyLessThanOrEqual( ...
                result.line.illuminationFactor, 1);
        end

        function testSNRDecreasesWithRange(testCase)
            result = iso_performance_snr_analysis('', ...
                testCase.ResultFolder, false);

            testCase.verifyLessThan(diff(result.annular.snr), 0);
            testCase.verifyLessThan(diff(result.line.snr), 0);
            testCase.verifyLessThan(result.annular.snr, result.line.snr);
        end

        function testExpectedFieldReduction(testCase)
            result = iso_performance_snr_analysis('', ...
                testCase.ResultFolder, false);

            testCase.verifyEqual(result.metrics.wDReductionFraction, ...
                1 - 120/195, 'AbsTol', 1e-12);
            testCase.verifyGreaterThan( ...
                result.metrics.outerFieldReductionFraction, 0.35);
            testCase.verifyLessThan( ...
                result.metrics.outerFieldReductionFraction, 0.40);
        end

        function testExportsExpectedArtifacts(testCase)
            outputFolder = tempname;
            mkdir(outputFolder);
            testCase.addTeardown(@() rmdir(outputFolder, 's'));

            iso_performance_snr_analysis(outputFolder, ...
                testCase.ResultFolder, true);

            testCase.verifyTrue(isfile(fullfile(outputFolder, ...
                'iso_performance_snr.png')));
            testCase.verifyTrue(isfile(fullfile(outputFolder, ...
                'iso_performance_snr.fig')));
            testCase.verifyTrue(isfile(fullfile(outputFolder, ...
                'iso_performance_snr_photon_budget.csv')));
            testCase.verifyTrue(isfile(fullfile(outputFolder, ...
                'iso_performance_snr_summary.csv')));
            testCase.verifyTrue(isfile(fullfile(outputFolder, ...
                'iso_performance_snr_effective_pulses.csv')));
        end
    end
end

classdef sensitivityMinimumWidthTest < matlab.unittest.TestCase
    % Analytic checks for first-crossing estimation without shape correction.
    methods (TestClassSetup)
        function addAnalysisPath(testCase)
            folder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(folder));
        end
    end
    methods (Test)
        function respectsDipImmediatelyBeforeCrossing(testCase)
            widths = [10; 20; 30; 40];
            probability = [0.4; 0.3; 0.6; 0.2];
            [minimum, gridMinimum] = sensitivity_min_width(widths, probability, 0.5);
            testCase.verifyEqual(minimum, 20 + 20/3, AbsTol=1e-12);
            testCase.verifyEqual(gridMinimum, 30);
        end
        function doesNotInterpolateAnUnattainedRequirement(testCase)
            widths = [10; 20; 30];
            probability = [0.2 0.1; 0.4 0.3; 0.3 0.5];
            [minimum, gridMinimum] = sensitivity_min_width(widths, probability, 0.6);
            testCase.verifyTrue(all(isnan(minimum)));
            testCase.verifyTrue(all(isnan(gridMinimum)));
        end
        function retainsDifferentFeasibilityForPairedBeams(testCase)
            widths = [10; 20; 30];
            probability = [0.2 0.1; 0.6 0.3; 0.8 0.4];
            [minimum, gridMinimum] = sensitivity_min_width(widths, probability, 0.5);
            testCase.verifyEqual(minimum(1), 17.5, AbsTol=1e-12);
            testCase.verifyTrue(isnan(minimum(2)));
            testCase.verifyEqual(gridMinimum(1), 20);
        end
        function reportsLowerSearchBoundaryWithoutExtrapolation(testCase)
            widths = [10; 20; 30];
            probability = [0.7; 0.8; 0.6];
            minimum = sensitivity_min_width(widths, probability, 0.6);
            testCase.verifyEqual(minimum, 10);
        end
        function rejectsUnorderedDivergenceGrid(testCase)
            testCase.verifyError(@() sensitivity_min_width([10; 30; 20], ...
                [0.1; 0.2; 0.3], 0.5), 'Sensitivity:Grid');
        end
    end
end

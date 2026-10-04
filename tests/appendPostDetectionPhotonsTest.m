classdef appendPostDetectionPhotonsTest < matlab.unittest.TestCase
    %APPENDPOSTDETECTIONPHOTONSTEST Tests continuous post-detection capture.

    methods (TestClassSetup)
        function addProjectFolder(testCase)
            projectFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                projectFolder));
        end
    end

    methods (Test)
        function testStartsImmediatelyAfterDetection(testCase)
            photons = empty_detection_window_photons();

            [actual, complete] = append_post_detection_photons( ...
                photons, (1:6)', (11:16)', (8:13)', 10, 13, 0.2, 0.01);

            testCase.verifyTrue(complete);
            testCase.verifyEqual(actual.signal_photons, (3:6)');
            testCase.verifyEqual(actual.backscatter_photons, (13:16)');
            testCase.verifyEqual(actual.pulse_steps, uint32((10:13)'));
            testCase.verifyEqual(actual.background_photons, 0.2, AbsTol=0);
            testCase.verifyEqual(actual.dark_count_photons, 0.01, AbsTol=0);
        end

        function testCaptureContinuesAcrossBlocks(testCase)
            photons = empty_detection_window_photons();
            [firstPart, firstComplete] = append_post_detection_photons( ...
                photons, [1; 0], [0.1; 0.2], [21; 22], ...
                21, 24, 0.3, 0.02);

            [actual, complete] = append_post_detection_photons( ...
                firstPart, [0; 4; 9], [0.3; 0.4; 0.5], [23; 24; 25], ...
                21, 24, 0.3, 0.02);

            testCase.verifyFalse(firstComplete);
            testCase.verifyTrue(complete);
            testCase.verifyEqual(actual.signal_photons, [1; 0; 0; 4]);
            testCase.verifyEqual(actual.backscatter_photons, ...
                [0.1; 0.2; 0.3; 0.4], AbsTol=0);
            testCase.verifyEqual(actual.pulse_steps, uint32((21:24)'));
        end

        function testDoesNotDuplicateAlreadyCapturedPulses(testCase)
            photons = empty_detection_window_photons();
            [firstPart, ~] = append_post_detection_photons( ...
                photons, [2; 3], [0.1; 0.2], [31; 32], ...
                31, 34, 0.3, 0.02);

            [actual, complete] = append_post_detection_photons( ...
                firstPart, [3; 4; 5], [0.2; 0.3; 0.4], [32; 33; 34], ...
                31, 34, 0.3, 0.02);

            testCase.verifyTrue(complete);
            testCase.verifyEqual(actual.signal_photons, [2; 3; 4; 5]);
            testCase.verifyEqual(actual.pulse_steps, uint32((31:34)'));
        end
    end
end

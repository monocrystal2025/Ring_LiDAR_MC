classdef eaPostDetectionIntegrationTest < matlab.unittest.TestCase
    %EAPOSTDETECTIONINTEGRATIONTEST Tests EA five-output capture interfaces.

    methods (TestClassSetup)
        function addProjectFolder(testCase)
            projectFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                projectFolder));
        end
    end

    methods (Test)
        function testRingCaptureStartsAfterDetection(testCase)
            regionRange = 101;
            detectionRange = 100;
            f = 5000;
            fasanD = 0.15;
            fasanDsmall = 0.001;
            omiga = 2*pi;
            targetAngle = fasanD / 2;
            targetPosition = detectionRange * ...
                [sin(targetAngle), 0, cos(targetAngle)];
            beam = [0, 0, 1];

            [detected, firstDetect, ~, ~, postPhotons] = MC_ring_snr( ...
                fasanD, fasanDsmall, beam, beam, beam, f, omiga, ...
                targetPosition, [0, 0, 0], 1, regionRange, 1000);
            expectedCount = post_detection_pulse_count( ...
                f, fasanD + 3*fasanDsmall, omiga, detectionRange, 0);

            testCase.verifyEqual(detected, 1);
            testCase.verifyNumElements(postPhotons.pulse_steps, expectedCount);
            testCase.verifyEqual(double(postPhotons.pulse_steps(1)), ...
                firstDetect + 1);
        end

        function testLineCaptureStartsAfterDetection(testCase)
            regionRange = 101;
            detectionRange = 100;
            f = 5000;
            fasanD = 0.15;
            fasanDsmall = 0.001;
            omiga = 2*pi;
            beam = [0, 0, 1];

            [detected, firstDetect, ~, ~, postPhotons] = MC_line_snr( ...
                fasanD, fasanDsmall, beam, beam, beam, f, omiga, ...
                [0, 0, detectionRange], [0, 0, 0], 1, regionRange, 1000);
            expectedCount = post_detection_pulse_count( ...
                f, fasanD + 3*fasanDsmall, omiga, detectionRange, 0);

            testCase.verifyEqual(detected, 1);
            testCase.verifyNumElements(postPhotons.pulse_steps, expectedCount);
            testCase.verifyEqual(double(postPhotons.pulse_steps(1)), ...
                firstDetect + 1);
        end

        function testPointCaptureStartsAfterDetection(testCase)
            regionRange = 101;
            detectionRange = 100;
            f = 5000;
            fasanD = 0.05;
            omiga = 2*pi;
            beam = [0, 0, 1];

            [detected, firstDetect, ~, ~, postPhotons] = MC_point_snr( ...
                fasanD, [], beam, beam, beam, f, omiga, ...
                [0, 0, detectionRange], [0, 0, 0], 1, regionRange);
            expectedCount = post_detection_pulse_count( ...
                f, fasanD + 3*0, omiga, detectionRange, 0);

            testCase.verifyEqual(detected, 1);
            testCase.verifyNumElements(postPhotons.pulse_steps, expectedCount);
            testCase.verifyEqual(double(postPhotons.pulse_steps(1)), ...
                firstDetect + 1);
        end

        function testDetectionInsideBlockUsesFollowingPulse(testCase)
            regionRange = 101;
            detectionRange = 100;
            f = 5000;
            fasanD = 0.15;
            fasanDsmall = 0.001;
            omiga = 2*pi;
            targetPosition = [0, 0, detectionRange];
            beamOff = [0, 0, 1];
            targetAngle = fasanD / 2;
            beamHit = [sin(targetAngle), 0, cos(targetAngle)];

            [detected, firstDetect, ~, ~, postPhotons] = MC_ring_snr( ...
                fasanD, fasanDsmall, beamHit, beamOff, beamHit, f, ...
                omiga, targetPosition, [0, 0, 0], 1, regionRange, 1000);
            expectedCount = post_detection_pulse_count( ...
                f, fasanD + 3*fasanDsmall, omiga, detectionRange, 0);

            testCase.verifyEqual(detected, 1);
            testCase.verifyGreaterThan(firstDetect, 0);
            testCase.verifyNumElements(postPhotons.pulse_steps, expectedCount);
            testCase.verifyEqual(double(postPhotons.pulse_steps(1)), ...
                firstDetect + 1);
        end
    end
end

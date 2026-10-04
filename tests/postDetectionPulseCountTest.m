classdef postDetectionPulseCountTest < matlab.unittest.TestCase
    %POSTDETECTIONPULSECOUNTTEST Tests the dynamic capture-length formula.

    methods (TestClassSetup)
        function addProjectFolder(testCase)
            projectFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                projectFolder));
        end
    end

    methods (Test)
        function testUsesDetectionRangeAboveFloor(testCase)
            f = 5000;
            baseAngle = 0.15 + 3e-3;
            omiga = 2*pi;
            detectionRange = 2000;
            targetSpeed = 30;
            targetWidth = 0.297;
            targetAngle = atan((sqrt(2)*targetWidth/2) / detectionRange);
            expected = ceil(f*(baseAngle + 2*targetAngle) / ...
                (omiga - targetSpeed/detectionRange));

            actual = post_detection_pulse_count( ...
                f, baseAngle, omiga, detectionRange, targetSpeed);

            testCase.verifyEqual(actual, expected);
        end

        function testClampsRangesBelowOneHundredMeters(testCase)
            belowFloor = post_detection_pulse_count( ...
                5000, 0.053, 2*pi, 25, 30);
            atFloor = post_detection_pulse_count( ...
                5000, 0.053, 2*pi, 100, 30);

            testCase.verifyEqual(belowFloor, atFloor);
        end

        function testSameDirectionSpeedIncreasesCaptureLength(testCase)
            stationary = post_detection_pulse_count( ...
                5000, 0.153, pi/2, 100, 0);
            moving = post_detection_pulse_count( ...
                5000, 0.153, pi/2, 100, 30);

            testCase.verifyGreaterThan(moving, stationary);
        end
    end
end

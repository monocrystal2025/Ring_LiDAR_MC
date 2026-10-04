classdef mcEaBlockSizeTest < matlab.unittest.TestCase
    %MCEABLOCKSIZETEST Verify the measured piecewise block-size policy.

    properties (TestParameter)
        blockCase = struct( ...
            'fiveMrad', {[5e-3, 200000]}, ...
            'aboveFiveMrad', {[5.001e-3, 100000]}, ...
            'twentyFiveMrad', {[25e-3, 100000]}, ...
            'aboveTwentyFiveMrad', {[25.001e-3, 50000]}, ...
            'fortyMrad', {[40e-3, 50000]}, ...
            'aboveFortyMrad', {[40.001e-3, 20000]}, ...
            'oneHundredFiftyMrad', {[150e-3, 20000]});
    end

    methods (TestClassSetup)
        function addProjectFolder(testCase)
            projectFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                projectFolder));
        end
    end

    methods (Test)
        function testMeasuredPiecewiseMapping(testCase, blockCase)
            testCase.verifyEqual(mc_ea_block_size(blockCase(1)), ...
                blockCase(2));
        end

        function testRejectsNonpositiveAngle(testCase)
            testCase.verifyError(@() mc_ea_block_size(0), ...
                'MATLAB:mc_ea_block_size:expectedPositive');
        end
    end
end

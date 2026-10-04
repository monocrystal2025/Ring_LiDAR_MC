classdef appendAndCheckPhotonWindowTest < matlab.unittest.TestCase
    %APPENDANDCHECKPHOTONWINDOWTEST Tests the shared detection-window logic.

    properties (TestParameter)
        eaSnrFunction = struct( ...
            'ring', 'MC_ring_snr', ...
            'line', 'MC_line_snr', ...
            'point', 'MC_point_snr');
        fourOutputSnrFunction = struct( ...
            'eaPointCompatibility', 'MC_point', ...
            'roiRasterRing', 'MC_ring_RA_snr_NEW', ...
            'roiSpiralRing', 'MC_ring_SP_snr_NEW', ...
            'roiRasterLine', 'MC_line_RA_snr_NEW', ...
            'roiSpiralLine', 'MC_line_SP_snr_NEW', ...
            'roiRasterPoint', 'MC_point_RA_snr_NEW', ...
            'roiSpiralPoint', 'MC_point_SP_snr_NEW');
    end

    methods (TestClassSetup)
        function addProjectFolder(testCase)
            projectFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                projectFolder));
        end
    end

    methods (Test)
        function testEarliestTailAndShortestSuccessfulWindow(testCase)
            empty = zeros(0, 1);
            signal = [1; 1; 3; 10];
            backscatter = zeros(4, 1);
            steps = (1:4)';

            [detected, firstDetect, photons] = ...
                append_and_check_photon_window( ...
                empty, empty, empty, signal, backscatter, steps, ...
                0, 0, 4, 2);

            testCase.verifyEqual(detected, 1);
            testCase.verifyEqual(firstDetect, 3);
            testCase.verifyEqual(photons.signal_photons, [1; 3], ...
                AbsTol=0);
            testCase.verifyEqual(photons.backscatter_photons, [0; 0], ...
                AbsTol=0);
            testCase.verifyEqual(photons.pulse_steps, uint32([2; 3]));
        end

        function testFixedModeUsesAvailableStartupPulses(testCase)
            empty = zeros(0, 1);

            [detected, firstDetect, photons] = ...
                append_and_check_photon_window( ...
                empty, empty, empty, 5, 0, 0, ...
                0, 0, 4, 2, 'fixed');

            testCase.verifyEqual(detected, 1);
            testCase.verifyEqual(firstDetect, 0);
            testCase.verifyEqual(photons.signal_photons, 5, AbsTol=0);
            testCase.verifyEqual(photons.pulse_steps, uint32(0));
        end

        function testFixedModeStoresFullWindowAtLimit(testCase)
            empty = zeros(0, 1);
            signal = [1; 1; 1; 2];
            steps = (1:4)';

            [detected, firstDetect, photons] = ...
                append_and_check_photon_window( ...
                empty, empty, empty, signal, zeros(4, 1), steps, ...
                0, 0, 4, 2, 'fixed');

            testCase.verifyEqual(detected, 1);
            testCase.verifyEqual(firstDetect, 4);
            testCase.verifyEqual(photons.signal_photons, signal, AbsTol=0);
            testCase.verifyEqual(photons.pulse_steps, uint32(steps));
        end

        function testFixedModeUsesExactlyMWhenCrossingBlocks(testCase)
            empty = zeros(0, 1);
            [~, ~, ~, carrySig, carryBackscatter, carryStep] = ...
                append_and_check_photon_window( ...
                empty, empty, empty, [1; 1; 1; 0], zeros(4, 1), ...
                (1:4)', 0, 0, 4, 2, 'fixed');

            [detected, firstDetect, photons] = ...
                append_and_check_photon_window( ...
                carrySig, carryBackscatter, carryStep, 2, 0, 5, ...
                0, 0, 4, 2, 'fixed');

            testCase.verifyEqual(detected, 1);
            testCase.verifyEqual(firstDetect, 5);
            testCase.verifyEqual(photons.signal_photons, [1; 1; 0; 2], ...
                AbsTol=0);
            testCase.verifyEqual(photons.pulse_steps, uint32((2:5)'));
        end

        function testPhotonComponentsAndTypes(testCase)
            empty = zeros(0, 1);
            signal = [2; 4];
            backscatter = [0.5; 1.5];
            steps = [20; 21];
            background = 0.25;
            darkCount = 0.01;

            [detected, firstDetect, photons] = ...
                append_and_check_photon_window( ...
                empty, empty, empty, signal, backscatter, steps, ...
                background, darkCount, 2, 1.5);

            testCase.verifyEqual(detected, 1);
            testCase.verifyEqual(firstDetect, 21);
            testCase.verifyEqual(photons.signal_photons, 4, AbsTol=0);
            testCase.verifyEqual(photons.backscatter_photons, 1.5, ...
                AbsTol=0);
            testCase.verifyEqual(photons.background_photons, background, ...
                AbsTol=0);
            testCase.verifyEqual(photons.dark_count_photons, darkCount, ...
                AbsTol=0);
            testCase.verifyClass(photons.signal_photons, 'double');
            testCase.verifyClass(photons.backscatter_photons, 'double');
            testCase.verifyClass(photons.background_photons, 'double');
            testCase.verifyClass(photons.dark_count_photons, 'double');
            testCase.verifyClass(photons.pulse_steps, 'uint32');
        end

        function testSuccessfulWindowCanCrossBlocks(testCase)
            empty = zeros(0, 1);
            [~, ~, ~, carrySig, carryBackscatter, carryStep] = ...
                append_and_check_photon_window( ...
                empty, empty, empty, [1; 1], [0; 0], [1; 2], ...
                0, 0, 3, 2);

            [detected, firstDetect, photons] = ...
                append_and_check_photon_window( ...
                carrySig, carryBackscatter, carryStep, 3, 0, 3, ...
                0, 0, 3, 2);

            testCase.verifyEqual(detected, 1);
            testCase.verifyEqual(firstDetect, 3);
            testCase.verifyEqual(photons.signal_photons, [1; 3], ...
                AbsTol=0);
            testCase.verifyEqual(photons.pulse_steps, uint32([2; 3]));
        end

        function testZeroSignalPulsesInsideWindowAreRetained(testCase)
            empty = zeros(0, 1);
            signal = [1; 0; 1; 2];
            backscatter = [0.1; 0.2; 0.3; 0.4];

            [detected, firstDetect, photons] = ...
                append_and_check_photon_window( ...
                empty, empty, empty, signal, backscatter, (4:7)', ...
                0, 0, 4, 1.8);

            testCase.verifyEqual(detected, 1);
            testCase.verifyEqual(firstDetect, 7);
            testCase.verifyEqual(photons.signal_photons, signal, AbsTol=0);
            testCase.verifyEqual(photons.backscatter_photons, ...
                backscatter, AbsTol=0);
            testCase.verifyEqual(photons.pulse_steps, uint32((4:7)'));
        end

        function testNoDetectionReturnsCanonicalEmptyValue(testCase)
            empty = zeros(0, 1);

            [detected, firstDetect, photons] = ...
                append_and_check_photon_window( ...
                empty, empty, empty, [0.1; 0], [1; 1], [7; 8], ...
                0.2, 0.01, 2, 2);

            testCase.verifyEqual(detected, 0);
            testCase.verifyTrue(isnan(firstDetect));
            testCase.verifyEqual(photons, ...
                empty_detection_window_photons());
        end

        function testEaSNRInterfacesHaveFiveOutputs(testCase, eaSnrFunction)
            testCase.verifyEqual(nargout(eaSnrFunction), 5);
        end

        function testOtherSNRInterfacesKeepFourOutputs( ...
                testCase, fourOutputSnrFunction)
            testCase.verifyEqual(nargout(fourOutputSnrFunction), 4);
        end

        function testDispatcherInterfacesKeepFourOutputs(testCase)
            testCase.verifyEqual(nargout('MC_func'), 5);
            testCase.verifyEqual(nargout('MC_func_ROI_NEW'), 4);
        end

        function testEaDescriptiveMatFileRoundTrip(testCase)
            matFile = [tempname, '.mat'];
            testCase.addTeardown(@() delete(matFile));
            detect_R = true;
            first_time_R = 2;
            first_encounter_time_R = 3;
            effective_pulses_R = empty_detection_window_photons();
            effective_pulses_R.signal_photons = [1.25; 2.5];
            effective_pulses_R.backscatter_photons = [0.1; 0.2];
            effective_pulses_R.background_photons = 0.01;
            effective_pulses_R.dark_count_photons = 2e-5;
            effective_pulses_R.pulse_steps = uint32([8; 9]);

            save(matFile, 'detect_R', 'first_time_R', ...
                'first_encounter_time_R', 'effective_pulses_R');
            loaded = load(matFile);

            testCase.verifyEqual(sort(fieldnames(loaded)), ...
                {'detect_R'; 'effective_pulses_R'; ...
                'first_encounter_time_R'; 'first_time_R'});
            testCase.verifyClass(loaded.detect_R, 'logical');
            testCase.verifyEqual(loaded.effective_pulses_R, ...
                effective_pulses_R);
        end

        function testRoiDescriptiveMatFileRoundTrip(testCase)
            matFile = [tempname, '.mat'];
            testCase.addTeardown(@() delete(matFile));
            detect_P_RA = false;
            first_time_P_RA = nan;
            first_encounter_time_P_RA = 1.5;
            effective_pulses_P_RA = empty_detection_window_photons();

            save(matFile, 'detect_P_RA', 'first_time_P_RA', ...
                'first_encounter_time_P_RA', 'effective_pulses_P_RA');
            loaded = load(matFile);

            testCase.verifyEqual(sort(fieldnames(loaded)), ...
                {'detect_P_RA'; 'effective_pulses_P_RA'; ...
                'first_encounter_time_P_RA'; 'first_time_P_RA'});
            testCase.verifyClass(loaded.detect_P_RA, 'logical');
            testCase.verifyEqual(loaded.effective_pulses_P_RA, ...
                effective_pulses_P_RA);
        end
    end
end

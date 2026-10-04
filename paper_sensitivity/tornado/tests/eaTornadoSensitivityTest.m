classdef eaTornadoSensitivityTest < matlab.unittest.TestCase
    %EATORNADOSENSITIVITYTEST Tests for the public tornado workflow.

    methods (TestClassSetup)
        function addProjectPaths(testCase)
            testFolder = fileparts(mfilename("fullpath"));
            moduleFolder = fileparts(testFolder);
            paperFolder = fileparts(moduleFolder);
            projectRoot = fileparts(paperFolder);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                projectRoot));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                paperFolder));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                moduleFolder));
        end
    end

    methods (Test)
        function interpolatesFirstRawCrossing(testCase)
            widths = [5; 10; 15];
            probability = [0.30; 0.70; 0.40];

            [minimum, grid, bracket] = ea_tornado_min_width( ...
                widths, probability, 0.60);

            testCase.verifyEqual(minimum, 8.75, AbsTol=1e-12);
            testCase.verifyEqual(grid, 10, AbsTol=1e-12);
            testCase.verifyEqual(bracket, [5, 10], AbsTol=1e-12);
        end

        function doesNotExtrapolateUnattainedRequirement(testCase)
            widths = [5; 10; 15];
            probability = [0.10; 0.45; 0.59];

            [minimum, grid, bracket] = ea_tornado_min_width( ...
                widths, probability, 0.60);

            testCase.verifyTrue(isnan(minimum));
            testCase.verifyTrue(isnan(grid));
            testCase.verifyTrue(all(isnan(bracket)));
        end

        function reportsLowerBoundaryWithoutExtrapolation(testCase)
            widths = [5; 10; 15];
            probability = [0.65; 0.70; 0.80];

            [minimum, grid, bracket] = ea_tornado_min_width( ...
                widths, probability, 0.60);

            testCase.verifyEqual(minimum, 5, AbsTol=1e-12);
            testCase.verifyEqual(grid, 5, AbsTol=1e-12);
            testCase.verifyEqual(bracket, [5, 5], AbsTol=1e-12);
        end

        function rejectsUnorderedWidthGrid(testCase)
            testCase.verifyError(@() ea_tornado_min_width( ...
                [5; 15; 10], [0.2; 0.5; 0.7], 0.6), ...
                "EATornado:WidthGrid");
        end

        function ranksEffectsAndFlagsInteriorExcursion(testCase)
            levels = eaTornadoSensitivityTest.syntheticLevels();

            effects = ea_tornado_build_effects(levels);

            testCase.verifyEqual(effects.Factor(1), "energy");
            testCase.verifyEqual(effects.Rank, [1; 2]);
            testCase.verifyTrue( ...
                effects.IntermediateOutsideEndpointRange(1));
            testCase.verifyFalse( ...
                effects.IntermediateOutsideEndpointRange(2));
        end

        function createsSeventeenCaseSharedBaselineDesign(testCase)
            outputDirectory = string(tempname);
            mkdir(outputDirectory);
            testCase.addTeardown(@() rmdir(outputDirectory, "s"));

            outputs = run_ea_tornado_sensitivity( ...
                Mode="simulate", N=3, BootstrapReplicates=5, ...
                UseParallel=false, TaskLimit=0, ...
                InitialWidthsMrad=[50, 100, 150], MaxWidthMrad=150, ...
                OutputDirectory=outputDirectory);

            testCase.verifyEqual(numel(outputs.design.cases), 17);
            testCase.verifyEqual(nnz([outputs.design.cases.factorIndex] == 0), 1);
            testCase.verifyEqual(size(outputs.design.position_m), [3, 3]);
            testCase.verifyEqual(nnz(outputs.design.required), 51);
        end

        function rejectsChangedCachedConfiguration(testCase)
            outputDirectory = string(tempname);
            mkdir(outputDirectory);
            testCase.addTeardown(@() rmdir(outputDirectory, "s"));
            run_ea_tornado_sensitivity(Mode="simulate", N=3, ...
                BootstrapReplicates=5, UseParallel=false, TaskLimit=0, ...
                InitialWidthsMrad=[50, 100, 150], MaxWidthMrad=150, ...
                OutputDirectory=outputDirectory, Seed=10);

            call = @() run_ea_tornado_sensitivity(Mode="simulate", N=3, ...
                BootstrapReplicates=5, UseParallel=false, TaskLimit=0, ...
                InitialWidthsMrad=[50, 100, 150], MaxWidthMrad=150, ...
                OutputDirectory=outputDirectory, Seed=11);

            testCase.verifyError(call, "EATornado:CacheSignature");
        end
    end

    methods (Test, TestTags={'Integration', 'Slow'})
        function completesSmokeWorkflow(testCase)
            outputDirectory = string(tempname);
            mkdir(outputDirectory);
            testCase.addTeardown(@() rmdir(outputDirectory, "s"));

            outputs = run_ea_tornado_sensitivity(Mode="smoke", ...
                N=2, BootstrapReplicates=10, UseParallel=false, ...
                OutputDirectory=outputDirectory);
            repeated = ea_tornado_analyze(outputDirectory, 10, ...
                outputs.config.bootstrapSeed);
            lowFeasibility = repeated.levels.BothFeasibleBootstrap < ...
                outputs.config.feasibilityThreshold;

            testCase.verifyTrue(outputs.verification.Passed);
            testCase.verifyTrue(isfile(outputs.figures.png));
            testCase.verifyTrue(isfile(fullfile(outputDirectory, ...
                "ea_tornado_effects.csv")));
            testCase.verifyEqual(height(outputs.analysis.levels), 20);
            testCase.verifyEqual(repeated.levels.KappaCI_Low, ...
                outputs.analysis.levels.KappaCI_Low, AbsTol=0);
            testCase.verifyTrue(all(isnan( ...
                repeated.levels.KappaCI_Low(lowFeasibility))));
        end
    end

    methods (Static, Access=private)
        function levels = syntheticLevels()
            Factor = repelem(["energy"; "speed"], 5);
            FactorLabel = repelem(["Pulse energy"; "Target speed"], 5);
            Scale = repmat([0.8; 0.9; 1.0; 1.1; 1.2], 2, 1);
            Value = [120; 135; 150; 165; 180; 24; 27; 30; 33; 36];
            Unit = repelem(["uJ"; "m/s"], 5);
            Kappa = [0.80; 1.02; 0.90; 0.94; 0.95; ...
                0.84; 0.86; 0.87; 0.88; 0.90];
            KappaCI_Low = Kappa - 0.01;
            KappaCI_High = Kappa + 0.01;
            baseline = [repmat(0.90, 5, 1); repmat(0.87, 5, 1)];
            DeltaKappa = Kappa - baseline;
            DeltaKappaCI_Low = DeltaKappa - 0.015;
            DeltaKappaCI_High = DeltaKappa + 0.015;
            BothFeasibleBootstrap = ones(10, 1);
            levels = table(Factor, FactorLabel, Scale, Value, Unit, ...
                Kappa, KappaCI_Low, KappaCI_High, DeltaKappa, ...
                DeltaKappaCI_Low, DeltaKappaCI_High, ...
                BothFeasibleBootstrap);
        end
    end
end

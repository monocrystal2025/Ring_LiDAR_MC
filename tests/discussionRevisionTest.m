classdef discussionRevisionTest < matlab.unittest.TestCase
    %DISCUSSIONREVISIONTEST Tests of physical definitions, sparse windows and parity.
    methods (TestClassSetup)
        function projectPath(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end
    methods (Test)
        function gapPreservesWindowDuration(testCase)
            [cfg,model] = unitPhotonModel();
            events = [0,1,2;5,1,2];
            actual = discussion2.scoreEvents(events,model,cfg,5);
            testCase.verifyFalse(isfinite(actual.firstTime));
            testCase.verifyEqual(actual.maxQ,2,'AbsTol',1e-12);
        end
        function twoSeparatedQualifiedClusters(testCase)
            [cfg,model] = unitPhotonModel();
            events = [0,1,1.5;1,1,0.5;4,1,2];
            actual = discussion2.scoreEvents(events,model,cfg,5);
            testCase.verifyEqual(actual.effectivePulses,3,'AbsTol',1e-12);
            testCase.verifyEqual(actual.clusters,2,'AbsTol',1e-12);
            testCase.verifyEqual(actual.qualifiedClusters,2,'AbsTol',1e-12);
            testCase.verifyEqual(actual.clusterSpan_ms,0.75,'AbsTol',1e-12);
        end
        function tinyTailDoesNotQualify(testCase)
            [cfg,model] = unitPhotonModel();
            actual = discussion2.scoreEvents([0,1,0.001;4,1,4],model,cfg,5);
            testCase.verifyEqual(actual.clusters,2,'AbsTol',1e-12);
            testCase.verifyEqual(actual.qualifiedClusters,1,'AbsTol',1e-12);
            testCase.verifyTrue(isnan(actual.clusterSpan_ms));
        end
        function firstDetectionExcludesFuturePulses(testCase)
            [cfg,model] = unitPhotonModel();
            actual = discussion2.scoreEvents([0,1,2;1,1,2;2,1,100],model,cfg,5);
            testCase.verifyEqual(actual.firstTime,1/5000,'AbsTol',1e-12);
            testCase.verifyEqual(actual.effectivePulses,2,'AbsTol',1e-12);
        end
        function absentSignalStaysMissing(testCase)
            [cfg,model] = unitPhotonModel();
            actual = discussion2.scoreEvents(zeros(0,3),model,cfg,5);
            testCase.verifyFalse(isfinite(actual.firstTime));
            testCase.verifyTrue(isnan(actual.effectivePulses));
        end
        function staticUsesVolumeSampling(testCase)
            cfg = discussion2.config("smoke");
            bank = struct('N',1,'u',0.5*ones(1,8));
            bank.u(8) = 1/8;
            [p,v,t] = discussion2.positions(cfg,bank,"static");
            testCase.verifyEqual(norm(p),cfg.system.maxRange_m/2,'AbsTol',1e-10);
            testCase.verifyEqual(v,zeros(1,3),'AbsTol',1e-12);
            testCase.verifyTrue(isinf(t));
        end
        function dynamicHorizonIsActualExit(testCase)
            cfg = discussion2.config("smoke");
            bank = discussion.bank(16,814);
            [p,v,t] = discussion2.positions(cfg,bank,"dynamic");
            finish = p+v.*t;
            boundary = min(abs(vecnorm(finish,2,2)-cfg.system.maxRange_m),abs(finish(:,3)));
            testCase.verifyLessThan(max(boundary),1e-8);
            testCase.verifyGreaterThanOrEqual(min(finish(:,3)),-1e-8);
        end
        function minimumKeepsMissesAndNonmonotonicity(testCase)
            curve = struct('widths_mrad',[30,60,120,195], ...
                'probability',[0,0;0.8,0.7;0.6,0.85;0.9,0.5],'N',100);
            actual = discussion2.minima(curve,0.8);
            testCase.verifyEqual(actual.widths_mrad,[60,120],'AbsTol',1e-12);
            testCase.verifyEqual(actual.kappa,2,'AbsTol',1e-12);
        end
        function fullRoeLineParity(testCase)
            [actual,expected] = compareOriginal("line","dynamic");
            testCase.verifyEqual(actual,expected,'AbsTol',1e-9);
        end
        function fullRoeRingParity(testCase)
            [actual,expected] = compareOriginal("ring","dynamic");
            testCase.verifyEqual(actual,expected,'AbsTol',1e-9);
        end
        function staticLineOneCycleParity(testCase)
            [actual,expected] = compareOriginal("line","static");
            testCase.verifyEqual(actual,expected,'AbsTol',1e-9);
        end
        function staticRingOneCycleParity(testCase)
            [actual,expected] = compareOriginal("ring","static");
            testCase.verifyEqual(actual,expected,'AbsTol',1e-9);
        end
    end
end

function [cfg,model] = unitPhotonModel()
cfg = discussion2.config("smoke");
model = discussion.photonModel(cfg,discussion.beam("line",120));
model.signalConstant = 1;
model.extinction = 0;
model.background = 0;
model.dark = 0;
model.backscatter(:) = 0;
end

function [actual,expected] = compareOriginal(name,kind)
cfg = discussion2.config("smoke");
bank = discussion.bank(12,90631);
beam = discussion.beam(name,120);
record = discussion2.geometry(cfg,beam,bank,kind);
out = discussion2.evaluate(cfg,beam,bank,kind);
actual = [out.metrics.firstTime,out.metrics.effectivePulses,out.metrics.clusters];
scan = discussion.scan(cfg,beam.wD_rad);
expected = nan(bank.N,3);
legacy = str2func("MC_"+name+"_snr");
for j = 1:bank.N
    initial = floor(bank.u(j,5)*scan.initialSupport)+1;
    if kind == "static"
        initial = 1;
    end
    maxStep = floor(record.horizon_s(j)*cfg.system.prf_Hz);
    if kind == "static"
        maxStep = scan.length-1;
    end
    [det,time,~,photons] = legacy(beam.wD_rad,beam.wd_rad, ...
        scan.path1,scan.spiral,scan.path2,5000,2*pi, ...
        record.position_m(j,:),record.velocity_m_s(j,:),initial,2000,25000,maxStep,'fixed');
    if det
        active = photons.signal_photons>0;
        expected(j,:) = [time/5000,nnz(active),nnz(active & [true;~active(1:end-1)])];
    else
        expected(j,1) = Inf;
    end
end
end

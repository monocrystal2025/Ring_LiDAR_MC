classdef discussionSimulationTest < matlab.unittest.TestCase
    %DISCUSSIONSIMULATIONTEST Physical parity and statistical edge cases.
    methods (TestClassSetup)
        function addProjectPath(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end
    methods (Test)
        function originalLineParity(testCase)
            [newTimes,oldTimes] = legacyComparison("line");
            testCase.verifyEqual(newTimes,oldTimes,'AbsTol',1e-10);
        end
        function originalRingParity(testCase)
            [newTimes,oldTimes] = legacyComparison("ring");
            testCase.verifyEqual(newTimes,oldTimes,'AbsTol',1e-10);
        end
        function annularRotationInvariance(testCase)
            cfg = discussion.config("smoke");
            b = discussion.beam("ring",120);
            bank = discussion.bank(32,7351);
            a = discussion.localCrossing(cfg,b,bank,0);
            c = discussion.localCrossing(cfg,b,bank,67);
            testCase.verifyEqual(a.metrics.maxQ,c.metrics.maxQ,'AbsTol',1e-9);
            testCase.verifyEqual(a.metrics.firstTime,c.metrics.firstTime,'AbsTol',1e-10);
        end
        function accumulationDominatesSingle(testCase)
            cfg = discussion.config("smoke");
            raw = discussion.simulate(cfg,discussion.beam("ring",195),discussion.bank(32,7301));
            testCase.verifyLessThanOrEqual(raw.metrics.firstTime,raw.metrics.singleTime);
            testCase.verifyGreaterThanOrEqual(raw.metrics.maxQ,raw.metrics.maxSingleQ-1e-10);
        end
        function nonmonotonicFrontier(testCase)
            times = [zeros(5,4);inf(5,4)];
            times(6:9,2) = 0;
            times(5,3) = inf;
            curve = struct('times',cat(3,times,times),'widths_mrad',[30,60,120,180],'N',10);
            f = discussion.frontier(curve,15,0.8);
            testCase.verifyEqual(f.lineWidth_mrad,60,'AbsTol',1e-12);
            testCase.verifyEqual(f.probability(:,1),[0.5;0.9;0.4;0.5],'AbsTol',1e-12);
        end
        function unreachableRemainsMissing(testCase)
            curve = struct('times',inf(10,4,2),'widths_mrad',[30,60,120,180],'N',10);
            f = discussion.frontier(curve,15,0.8);
            testCase.verifyTrue(isnan(f.kappa));
            testCase.verifyEqual(f.lineStatus,"unattained-on-grid");
        end
        function zeroDiscordanceStillHasUncertainty(testCase)
            cfg = discussion.config("smoke");
            stat = discussion.paired(zeros(100,1),zeros(100,1),cfg);
            testCase.verifyLessThan(stat.low,0);
            testCase.verifyGreaterThan(stat.high,0);
        end
        function photonSolidAngleNormalization(testCase)
            cfg = discussion.config("smoke");
            ring = discussion.photonModel(cfg,discussion.beam("ring",120));
            line = discussion.photonModel(cfg,discussion.beam("line",120));
            testCase.verifyEqual(ring.omegaBeam/line.omegaBeam,pi,'AbsTol',1e-10);
            testCase.verifyGreaterThan(ring.omegaMatched,line.omegaMatched);
        end
    end
end

function [newTimes,oldTimes] = legacyComparison(name)
cfg = discussion.config("smoke");
bank = discussion.bank(16,90210);
b = discussion.beam(name,120);
scan = discussion.scan(cfg,b.wD_rad);
new = discussion.simulate(cfg,b,bank);
newTimes = new.metrics.firstTime;
u = bank.u;
az = 2*pi*u(:,1);
cz = u(:,2);
p = 2000*[sqrt(1-cz.^2).*cos(az),sqrt(1-cz.^2).*sin(az),cz];
vaz = 2*pi*u(:,3);
vc = 2*u(:,4)-1;
v = 30*[sqrt(1-vc.^2).*cos(vaz),sqrt(1-vc.^2).*sin(vaz),vc];
flip = sum(p.*v,2)>0;
v(flip,:) = -v(flip,:);
oldTimes = inf(bank.N,1);
legacy = str2func("MC_"+name+"_snr");
for j = 1:bank.N
    initial = floor(u(j,5)*scan.initialSupport)+1;
    [det,time] = legacy(.12,.001,scan.path1,scan.spiral,scan.path2,5000,2*pi, ...
        p(j,:),v(j,:),initial,2000,25000,75000,'fixed');
    if det
        oldTimes(j) = time/5000;
    end
end
end

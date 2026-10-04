function result = tmp_1(mode,options)
%TMP_1 Manuscript-aligned sensitivity of angular width, blind volume and pulses.
%   tmp_1                    simulates and plots the revised preview.
%   tmp_1("paper")           uses 3000 moving and 10000 static targets.
%   tmp_1("preview",PlotOnly=true) redraws the saved revision.
% PNGs are written directly to this workspace root. See discussion_revision.md.
arguments
    mode (1,1) string {mustBeMember(mode,["smoke","preview","paper"])} = "preview"
    options.N (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.StaticN (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.PlotOnly (1,1) logical = false
end
cfg = discussion2.config(mode);
if options.N>0
    cfg.N = options.N;
end
if options.StaticN>0
    cfg.staticN = options.StaticN;
end
if options.PlotOnly
    saved = load(fullfile(cfg.outputRoot,"tmp_1.mat"),'result');
    result = saved.result;
else
    result = discussion2.computeStudy(cfg);
end
result = discussion2.addCoverage(result);
result = discussion2.addRequirements(result);
figures = discussion2.plotStudy(result);
discussion2.exportFigures(result.config,"tmp_1",figures);
end

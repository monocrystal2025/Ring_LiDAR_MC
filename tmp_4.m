function result = tmp_4(mode,options)
%TMP_4 Background-dependent Pd(wD) and receiver acceptance controls.
%   tmp_4 runs the revised background/width simulations and saves root PNGs.
%   tmp_4("preview",PlotOnly=true) redraws saved simulation results.
% Figure 1: four background levels, shape-matched receivers.
% Figure 2: receiver controls. Figure 3: archived 15 s engineering checks.
arguments
    mode (1,1) string {mustBeMember(mode,["smoke","preview","paper"])} = "preview"
    options.N (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.PlotOnly (1,1) logical = false
end
cfg = discussion2.config(mode);
if options.N>0
    cfg.N = options.N;
end
if options.PlotOnly
    saved = load(fullfile(cfg.outputRoot,"tmp_4.mat"),'result');
    result = saved.result;
else
    result = discussion2.computeBackground(cfg);
end
figures = discussion2.plotBackground(result);
discussion2.exportFigures(result.config,"tmp_4",figures);
end

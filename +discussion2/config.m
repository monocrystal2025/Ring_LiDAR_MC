function cfg = config(mode)
%CONFIG Settings for manuscript-aligned Discussion revisions.
arguments
    mode (1,1) string {mustBeMember(mode,["smoke","preview","paper"])} = "preview"
end
cfg = discussion.config(mode);
cfg.version = "discussion-revision-2.0";
cfg.outputRoot = fullfile(cfg.root,"discussion_results","revision",mode);
cfg.cacheRoot = fullfile(cfg.outputRoot,"cache");
cfg.N = 256;
cfg.staticN = 512;
cfg.widthGrid_mrad = unique([30:15:240,270,300]);
cfg.maxTime_s = Inf; % Overall detection before sphere/ground exit, as MC_EA.
cfg.targetProbability = 0.8;
cfg.bootstrapN = 400;
cfg.blockSize = 25000;
cfg.clusterQThreshold = 1; % Qualified cluster: accumulated squared SNR >= 1.
if mode == "smoke"
    cfg.N = 12;
    cfg.staticN = 24;
    cfg.widthGrid_mrad = [60,120,195];
elseif mode == "paper"
    cfg.N = 3000;
    cfg.staticN = 10000;
    cfg.widthGrid_mrad = 15:5:325;
    cfg.bootstrapN = 2000;
end
end

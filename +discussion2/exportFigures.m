function exportFigures(cfg,name,figures)
%EXPORTFIGURES Root-level PNGs plus archived vector PDF and editable FIG.
arguments
    cfg (1,1) struct
    name (1,1) string
    figures
end
if ~isfolder(cfg.outputRoot)
    mkdir(cfg.outputRoot);
end
for j = 1:numel(figures)
    stem = name+"_fig"+j;
    drawnow;
    exportgraphics(figures(j),fullfile(cfg.root,stem+".png"),'Resolution',600);
    exportgraphics(figures(j),fullfile(cfg.outputRoot,stem+".pdf"),'ContentType','vector');
    savefig(figures(j),fullfile(cfg.outputRoot,stem+".fig"));
    fprintf('PNG: %s\n',fullfile(cfg.root,stem+".png"));
end
end

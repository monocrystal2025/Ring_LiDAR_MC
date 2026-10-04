function finish(cfg,name,result,figures,summary)
%FINISH Save reproducible data, summary CSV, vector PDF, PNG, and FIG.
arguments
    cfg (1,1) struct
    name (1,1) string
    result (1,1) struct
    figures
    summary table
end
if ~isfolder(cfg.outputRoot)
    mkdir(cfg.outputRoot);
end
result.config = cfg;
result.implementationSHA256 = discussion.fingerprint();
result.createdAt = string(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
result.interpretation = "Simulated expected-SNR feasibility; not calibrated stochastic Pd/Pfa.";
result.statisticalNote = "Preview is exploratory. Paper conclusions require declared grid, adequate N and independent validation.";
save(fullfile(cfg.outputRoot,name+".mat"),'result','-v7.3');
writetable(summary,fullfile(cfg.outputRoot,name+"_summary.csv"));
for j = 1:numel(figures)
    stem = fullfile(cfg.outputRoot,name+"_fig"+j);
    drawnow;
    exportgraphics(figures(j),stem+".pdf",'ContentType','vector');
    exportgraphics(figures(j),stem+".png",'Resolution',600);
    savefig(figures(j),stem+".fig");
end
fprintf('%s finished (%s, N=%d). Outputs: %s\n',name,cfg.mode,cfg.N,cfg.outputRoot);
end

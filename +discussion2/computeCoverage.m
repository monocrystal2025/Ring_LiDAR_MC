function result = computeCoverage(cfg)
%COMPUTECOVERAGE Separate the paper's geometric blind zone from SNR blindness.
arguments
    cfg (1,1) struct
end
cases = discussion2.cases(cfg);
indices = [1,4,5,6,7];
curves = cell(numel(indices),1);
bank = discussion.bank(cfg.staticN,cfg.seed+1701);
for j = 1:numel(indices)
    fprintf('Coverage parameter case %d\n',indices(j));
    curves{j} = discussion2.coverageCurve(cases(indices(j)).config,bank);
end
result = struct('caseIndices',indices,'curves',{curves},'config',cfg);
save(fullfile(cfg.outputRoot,"tmp_1_geometric_coverage.mat"),'result','-v7.3');
end

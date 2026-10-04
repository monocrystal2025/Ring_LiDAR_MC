function finish_unified_robustness()
%FINISH_UNIFIED_ROBUSTNESS Analyze, export, and verify a completed simulation.
folder=fileparts(mfilename('fullpath'));addpath(folder,fileparts(folder));
analyze_unified_robustness;
plot_unified_robustness;
write_unified_manuscript;
verify_unified_robustness;
tests=runtests(fullfile(fileparts(folder),'tests'));
assertSuccess(tests);writetable(table(tests),fullfile(folder,'results','threshold_tests.csv'));
fprintf('UNIFIED_DELIVERABLES_COMPLETE\n');
end

%RUN_EA_SENSITIVITY_FORMAL Run the publication-scale N=3000 workflow.
%
% Results are isolated from the quick preview. Compatible checkpoints are
% resumed automatically, so this script can be rerun after interruption.

moduleDir = fileparts(mfilename("fullpath"));
addpath(moduleDir);
formalRoot = fullfile(moduleDir, "formal_output");
if ~isfolder(formalRoot)
    mkdir(formalRoot);
end

logFile = fullfile(formalRoot, "formal_run.log");
diary(logFile);
diary on;
diaryCleanup = onCleanup(@() diary("off")); %#ok<NASGU>

fprintf("Formal ROE sensitivity run started: %s\n", ...
    string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")));
fprintf("Output root: %s\n", formalRoot);
fprintf("Screen N = 3000; refine N = 3000; paired bootstrap = 200\n");

outputs = run_ea_sensitivity_analysis( ...
    "Mode", "full", ...
    "OutputRoot", formalRoot, ...
    "ScreenN", 3000, ...
    "RefineN", 3000, ...
    "BootstrapN", 200, ...
    "UseParallel", true, ...
    "Force", false);

completionFile = fullfile(formalRoot, "formal_run_complete.mat");
save(completionFile, "outputs", "-v7.3");
fprintf("Formal ROE sensitivity run completed: %s\n", ...
    string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")));
fprintf("Figure: %s\n", ...
    fullfile(formalRoot, "figures", "ea_sensitivity_main.png"));

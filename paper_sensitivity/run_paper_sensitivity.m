function results = run_paper_sensitivity(options)
%RUN_PAPER_SENSITIVITY Reproduce the redesigned Results-section analysis.
%   run_paper_sensitivity                 % resume N=512 joint preview, analyze
%   run_paper_sensitivity(JointN=3000)     % larger joint validation run
%   plot_paper_sensitivity                % redraw saved results, no simulation
arguments
    options.JointN (1,1) double {mustBeInteger,mustBePositive} = 512
    options.JointScenarios (1,1) double {mustBeInteger,mustBePositive} = 24
    options.JointStaticN (1,1) double {mustBeInteger,mustBePositive} = 2048
    options.BootstrapReplicates (1,1) double {mustBeInteger,mustBePositive} = 1000
    options.UseParallel (1,1) logical = true
end
folder=fileparts(mfilename('fullpath'));root=fileparts(folder);
addpath(root,folder);
out=fullfile(folder,'results',sprintf('joint_N%d',options.JointN));
run_joint_sensitivity(N=options.JointN,ScenarioCount=options.JointScenarios, ...
    UseParallel=options.UseParallel,OutputDirectory=out);
jointFile=fullfile(out,'joint_results.mat');
if options.JointStaticN>options.JointN
    if options.JointN==512
        staticOut=fullfile(folder,'results',sprintf('joint_static_N%d',options.JointStaticN));
    else
        staticOut=fullfile(folder,'results',sprintf('joint_N%d_static_N%d',options.JointN,options.JointStaticN));
    end
    extend_joint_static(InputFile=jointFile,OutputDirectory=staticOut, ...
        StaticN=options.JointStaticN,UseParallel=options.UseParallel);
    jointFile=fullfile(staticOut,'joint_results.mat');
end
results=analyze_paper_sensitivity(JointFile=jointFile, ...
    BootstrapReplicates=options.BootstrapReplicates);
plot_paper_sensitivity();
write_sensitivity_manuscript();
verify_paper_sensitivity();
end

function write_sensitivity_manuscript(options)
%WRITE_SENSITIVITY_MANUSCRIPT Write data-linked Chinese/English Results text.
arguments
    options.ResultsFile (1,1) string = ""
    options.OutputFile (1,1) string = ""
end
folder=fileparts(mfilename('fullpath'));
if options.ResultsFile==""
    options.ResultsFile=fullfile(folder,'results','paper_sensitivity.mat');
end
if options.OutputFile==""
    options.OutputFile=fullfile(folder,'manuscript_final.md');
end
saved=load(options.ResultsFile,'results');
r=saved.results;
t=r.joint(r.joint.Scenario>0,:);
base=r.oat(r.oat.Case=="energy_s03",:);
pg=100*[min(t.PdGain),max(t.PdGain)];
bg=100*[min(t.BlindReduction),max(t.BlindReduction)];
pdSupported=nnz(t.PdSimLow>0);
blindSupported=nnz(t.BlindSimLow>0);
fid=fopen(options.OutputFile,'w','n','UTF-8');
assert(fid>=0,'Sensitivity:Write','Cannot open the manuscript output.');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Results 最后一节建议稿\n\n');
fprintf(fid,'建议题名：**Parameter Sensitivity and Joint-Perturbation Robustness**。插入原稿 Pulse Utilization and Multi-Cluster Returns 之后、Discussion 之前。以下正文对应 figures_final/paper_parameter_sensitivity，图号可顺延为 Fig. 10。\n\n');
fprintf(fid,'## 本次运行的定量结果\n\n');
fprintf(fid,'- 单因素：每配置 N=%d，完整保留原有参数范围及不可达端点。\n',r.manifest.N);
fprintf(fid,'- 联合扰动：%d组四参数拉丁超立方组合，另加一个基准；每组合动态 N=%d、静态 N=%d。\n', ...
    r.jointMetadata.ScenarioCount,r.jointMetadata.N,r.jointMetadata.StaticN);
fprintf(fid,'- 基准最小角宽：线形 %.2f mrad、环形 %.2f mrad，κ=%.3f，配对95%%区间 %.3f–%.3f。\n', ...
    base.LineWidth,base.RingWidth,base.Kappa,base.KappaLow,base.KappaHigh);
fprintf(fid,'- 联合组合中的限时检出比例差：%.2f–%.2f 个百分点；静态盲区差：%.2f–%.2f 个百分点，点估计均偏向环形。\n',pg(1),pg(2),bg(1),bg(2));
fprintf(fid,'- 对25配置×2指标作联合同时区间控制后，%d/%d组的检出比例差区间下限大于0，%d/%d组的盲区差区间下限大于0。序号0基准未计入这些组数。\n\n', ...
    pdSupported,height(t),blindSupported,height(t));
fprintf(fid,'## 中文正文\n\n');
fprintf(fid,'为考察环形扫描优势对工作条件的依赖性，进一步开展单因素参数扫描及多参数联合扰动分析。单因素分析分别改变单脉冲能量、目标反射率、目标速度及大气消光系数，其余参数保持基准值。每个配置使用%d个目标样本，并以15 s内检出比例不低于60%%作为统一要求，比较两种光束所需的最小特征发散角。最小角宽由离散曲线的首次达标位置在相邻角宽点之间插值估计，未达标的配置明确标注为NA。大气参数变化时保持消光与后向散射系数之比为50 sr。\n\n',r.manifest.N);
fprintf(fid,'如图10(a)–(d)所示，提高脉冲能量或目标反射率降低了两种方案的角宽需求；在本研究的扫描范围内，提高目标速度或增加大气消光与后向散射则提高了所需角宽。基准条件下，线形和环形光束的最小角宽分别为%.2f和%.2f mrad，对应κ=%.3f，即约%.1f%%的角宽降低。速度和大气参数扫描中，环形光束在全部已测试水平上保持较小的角宽需求；提高回波强度后两种方案逐渐接近，说明能量充分时相对角宽收益减弱。弱回波端点仍可能无法满足规定要求，因此结果并不意味着环形光束在任意工作条件下都更优。\n\n', ...
    base.LineWidth,base.RingWidth,base.Kappa,100*(1-base.Kappa));
fprintf(fid,'为避免仅依据逐个改变参数的结果判断稳健性，采用拉丁超立方抽样在上述四个参数各为基准值±20%%的区间内构造%d组联合扰动组合，并额外计算基准配置。在100 mrad的共同角宽下，分别使用%d个动态目标和%d个体积均匀分布的静态目标，比较15 s内检出比例及一个完整扫描周期的静态盲区比例。图10(e)–(f)给出按目标配对重采样得到的差值及近似95%%同时区间。全部%d组组合的限时检出比例差和盲区降低量的点估计均为正，分别为%.2f–%.2f和%.2f–%.2f个百分点；同时区间下限大于零的组合分别为%d组和%d组。\n\n', ...
    height(t),r.jointMetadata.N,r.jointMetadata.StaticN,height(t),pg(1),pg(2),bg(1),bg(2),pdSupported,blindSupported);
fprintf(fid,'这些结果支持环形扫描的角宽、限时检测与覆盖收益并非仅出现在单一基准配置，而是在所测试的参数水平和联合扰动组合中具有一定持续性。该现象与环形几何提供重复相遇机会、同时受到局部能量稀释约束的机制一致。本节结论限定于ROE模式、所测试的参数与目标分布及期望光子数加权SNR判据；24组组合不代表对连续四维参数区间的穷举，亦不构成Sobol全局敏感度分解。\n\n');
fprintf(fid,'## English text\n\n');
fprintf(fid,'To examine the dependence of the annular-beam advantage on operating conditions, one-factor parameter sweeps were complemented by a joint-perturbation experiment. Pulse energy, target reflectivity, target speed, and atmospheric extinction were varied individually about the baseline. Each configuration contained %d target realizations. The minimum characteristic divergence was estimated for a common requirement of a detection fraction of at least 60%% within 15 s, using interpolation between the adjacent divergence samples at the first threshold crossing. Unattained requirements are explicitly marked as NA. The extinction-to-backscatter ratio was held at 50 sr throughout the atmospheric sweep.\n\n',r.manifest.N);
fprintf(fid,'As shown in Fig. 10(a)–(d), increasing pulse energy or target reflectivity reduced the angular-width demand of both geometries. Within the tested range, increasing target speed or atmospheric extinction and backscatter increased the required angular width. At the baseline, the line and annular beams required %.2f and %.2f mrad, respectively, giving κ=%.3f and an estimated angular-width reduction of %.1f%%. The annular beam retained a smaller required angular width at every tested speed and atmospheric level. With sufficiently strong returns, the two geometries approached similar angular requirements, whereas weak-return configurations could fail to meet the prescribed requirement. Thus, the relative benefit is conditional on the encounter–SNR balance.\n\n', ...
    base.LineWidth,base.RingWidth,base.Kappa,100*(1-base.Kappa));
fprintf(fid,'For the joint experiment, %d Latin-hypercube parameter combinations were sampled within ±20%% of the baseline values of all four factors, with the baseline evaluated separately. At a common divergence of 100 mrad, %d moving-target and %d stationary-target realizations were used per configuration. Figures 10(e)–(f) report the paired differences in the 15 s detection fraction and the one-cycle static blind-zone fraction, with approximate 95%% simultaneous bootstrap intervals. All %d perturbed configurations had positive point estimates for both the detection-fraction gain and blind-zone reduction, ranging from %.2f to %.2f and from %.2f to %.2f percentage points, respectively. The simultaneous lower limits exceeded zero in %d configurations for the detection-fraction gain and %d configurations for the blind-zone reduction.\n\n', ...
    height(t),r.jointMetadata.N,r.jointMetadata.StaticN,height(t),pg(1),pg(2),bg(1),bg(2),pdSupported,blindSupported);
fprintf(fid,'These results support benefits extending beyond the single baseline configuration, while identifying weak-return conditions where the advantage can diminish or disappear. The findings are consistent with repeated interrogation opportunities competing with energy dilution. They are restricted to the ROE model, the tested parameter combinations and target distributions, and the adopted expected-photon SNR criterion. The joint sample does not exhaust the continuous parameter box or constitute a variance-based global sensitivity decomposition.\n\n');
fprintf(fid,'## Figure caption\n\n');
fprintf(fid,'**Fig. 10. Parameter sensitivity and joint-perturbation robustness of annular-beam scanning in ROE surveillance.** (a)–(d) Minimum characteristic divergence required for a detection fraction of at least 60%% within 15 s, under variations in pulse energy, reflectivity, speed, and atmospheric extinction, respectively. Blue circles and orange triangles denote line and annular beams; vertical dotted lines mark the baseline values. Energy is shown on a logarithmic axis. NA denotes a requirement not attained over the simulated 5–400 mrad range and is not a numerical angular-width value. (e) Detection-fraction gain P_D,R−P_D,L and (f) static blind-zone reduction B_L−B_R at 100 mrad for %d joint parameter perturbations within ±20%% of the baseline; scenario 0 is the separate baseline. Positive contrasts favor annular illumination. Whiskers in (a)–(d) are pointwise 95%% paired-bootstrap intervals; those in (e)–(f) are approximate simultaneous 95%% intervals across all 25 configurations and both response contrasts. OAT results use N=%d targets per configuration; joint results use N_moving=%d and N_static=%d. All intervals quantify Monte Carlo sampling uncertainty, not angular-grid or model uncertainty.\n\n', ...
    height(t),r.manifest.N,r.jointMetadata.N,r.jointMetadata.StaticN);
fprintf(fid,'## 使用边界与投稿建议\n\n');
fprintf(fid,'- 这是已实际运行的联合扰动小批量结果。若论文需要对更小的差值作强结论，可统一增加联合动态样本；不能只追加不显著的组合直到显著。\n');
fprintf(fid,'- 24组点估计均有利，不等于连续参数区间中所有组合都更优，也不等于24组效应均显著。文中的区间支持组数必须保留或改成更谨慎的定性表述。\n');
fprintf(fid,'- 角宽插值不能消除5 mrad网格误差；高能量下差值小于网格间距时应表述为两种方案接近。\n');
fprintf(fid,'- 不应把角宽降低百分比直接写成接收立体角或SNR的同比例变化。\n');
fprintf(fid,'- 不把本节推论扩展为ROI模式、不同扫描硬件参数或实测虚警性能均已验证。\n');
fprintf(fid,'- 单因素分析的局限及联合设计的动机可引用 Saltelli and Annoni (2010), DOI 10.1016/j.envsoft.2010.04.012；详见README中的[JRC原始记录](https://publications.jrc.ec.europa.eu/repository/handle/JRC49564)。\n');
end

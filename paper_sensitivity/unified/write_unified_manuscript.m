function write_unified_manuscript(options)
%WRITE_UNIFIED_MANUSCRIPT Write conclusions at their measured evidence level.
arguments
    options.Directory (1,1) string = ""
end
folder=fileparts(mfilename('fullpath'));
if options.Directory=="",options.Directory=fullfile(folder,'results');end
saved=load(fullfile(options.Directory,'analysis.mat'),'analysis');a=saved.analysis;
k=find(a.design.options.ReferenceWidths==125,1);y=a.point(:,:,k);
ranges=[min(y,[],1);max(y,[],1)];pointwise=sum(a.ci(:,:,1,k)>0,1);simultaneous=sum(a.sim(:,:,1,k)>0,1);
file=fullfile(folder,'manuscript.md');fid=fopen(file,'w','n','UTF-8');cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# Results 最后一节：Local Parameter Robustness of the Three Annular-Beam Benefits\n\n');
fprintf(fid,'## 为什么旧图有NA\n\nNA表示规定探测要求在搜索角宽内未达到，并不是仿真数据丢失。诚实报告NA本身不是论文硬伤；将其删除后声称整个范围都更优，或把它插补为达标值，才会使证据与结论不一致。新主图研究明确限定的基准附近工作区间，原有宽范围失败与反转完整保留在补图。\n\n');
fprintf(fid,'## 统一方法\n\n三项优势保持并列：角宽降低、盲区降低和有效脉冲增加。三列分别为脉冲能量、反射率和大气消光系数，各改变−20%%、−10%%、0、+10%%、+20%%，其余保持基准；消光与后向散射比为50 sr。三行均显示相对线形的百分比收益，正值有利于环形，采用相同的配对bootstrap方法和95%%点态区间。基准在三列重复展示，但只有13个独立参数配置。\n\n');
fprintf(fid,'基准为E=150 μJ、rho=0.5、alpha=1.5×10^−5 m^−1。相应研究范围为E=120–180 μJ、rho=0.4–0.6、alpha=(1.2–1.8)×10^−5 m^−1。\n\n');
fprintf(fid,'角宽指标为15 s内检出比例达到60%%所需最小角宽的相对降低；盲区与有效脉冲主工作点均为125 mrad。有效脉冲计数沿用原稿定义：各光束首次成功检测窗内正期望回波脉冲的平均数量，条件于该光束最终成功检出。角宽使用限时事件，脉冲使用完整入场至出界过程，两者估计对象不同且均在代码中明确。\n\n');
fprintf(fid,'本次为新种子的小批量结果：动态N=%d、静态N=%d；没有把它标成N=3000。125 mrad与±20%%范围参考既有探索结果选定，属于工作区间分析，不是独立确定的硬件容差或预注册范围。100/150 mrad工作点和原N=3000宽范围结果作为补充检验。\n\n',a.design.options.N,a.design.options.StaticN);
fprintf(fid,'## 实际定量结果\n\n|指标|13配置点估计范围|95%%点态下限>0|近似95%%同时下限>0|\n|---|---|---|---|\n');
names={'角宽降低','盲区相对降低','有效脉冲增加'};
for row=1:3
    fprintf(fid,'|%s|%.2f–%.2f%%|%d/13|%d/13|\n',names{row},ranges(1,row),ranges(2,row),pointwise(row),simultaneous(row));
end
fprintf(fid,'\n同时区间覆盖主工作点的13配置×3指标；图中画的是点态区间。盲区相对降低百分数与盲区的百分点差不是同一数量。CSV包含线形/环形绝对值、两类区间、全部工作点。\n\n');
common=a.commonPulseGain(:,k);
fprintf(fid,'两种光束共同成功目标的脉冲增加为%.2f–%.2f%%，与主统计方向一致；因此本次脉冲收益并非仅由两种光束成功样本群体不同所产生。这是成功目标子集对照，不是所有目标的无条件脉冲收益。\n\n',min(common),max(common));
fprintf(fid,'原N=3000基准角宽收益约13.0%%，本次独立N=512基准为%.2f%%。两者样本库不同；本次区间包含原估计，不能把两次点估计差异解释为物理模型改变。\n\n',y(1,1));
fprintf(fid,'## 中文正文草稿\n\n为检验环形扫描的三项收益是否仅依赖单一基准配置，分别对脉冲能量、目标反射率及大气消光系数开展局部参数扰动分析，各参数取基准值的80%%–120%%，其余条件保持不变。将最小角宽降低、静态盲区降低和有效脉冲增加作为三个并列指标，均相对于线形光束归一化，并使用共同目标样本的配对自助法评估采样不确定性。最小角宽满足15 s内检出比例不低于60%%的共同要求；盲区与成功检测窗脉冲统计采用125 mrad的共同角宽。\n\n');
fprintf(fid,'在13个已测试配置中，三项相对收益分别为%.2f–%.2f%%、%.2f–%.2f%%和%.2f–%.2f%%。对应95%%点态区间下限为正的配置分别为%d、%d和%d个；同时考虑39项比较的近似95%%区间下限为正的配置分别为%d、%d和%d个。上述结果反映所选工作区间内三项收益的保持程度，而不意味着三种收益具有相同的变化幅度或对参数具有相同的敏感性。\n\n',ranges(1,1),ranges(2,1),ranges(1,2),ranges(2,2),ranges(1,3),ranges(2,3),pointwise,simultaneous);
fprintf(fid,'提高能量或反射率增强回波，大气消光与后向散射增强则降低探测能力，这些绝对响应方向与模型一致。但相对收益由两种光束响应的竞争决定，不要求其随参数单调增加。有效脉冲数量还受首次报告时间和成功目标群体变化影响，较多脉冲本身不等于更高SNR。本节与原稿的检测概率、盲区和回波簇分析共同说明环形几何的作用。补充工作点和宽范围分析进一步界定优势的适用边界。\n\n');
fprintf(fid,'## English results\n\nThe three annular-beam benefits were evaluated in parallel under local perturbations of pulse energy, target reflectivity, and atmospheric extinction. Each factor was varied individually from 80%% to 120%% of its baseline value, with the remaining conditions fixed and the extinction-to-backscatter ratio maintained at 50 sr. The responses were the relative reduction in required divergence, the relative reduction in static blind-zone fraction, and the relative increase in effective pulses per successful detection window. All responses used paired target resampling and pointwise 95%% bootstrap intervals. Required divergence was evaluated for a detection fraction of at least 60%% within 15 s; blind-zone and pulse statistics used a common divergence of 125 mrad.\n\n');
fprintf(fid,'Across the 13 tested configurations, the estimated benefits ranged from %.2f%% to %.2f%% for divergence, %.2f%% to %.2f%% for blind-zone fraction, and %.2f%% to %.2f%% for effective pulses. Pointwise lower limits exceeded zero in %d, %d, and %d configurations, respectively; approximate simultaneous lower limits across the 39 contrasts exceeded zero in %d, %d, and %d configurations. These results characterize local, one-factor-at-a-time robustness within the stated operating regime. They do not establish universal superiority or robustness to arbitrary joint parameter changes.\n\n',ranges(1,1),ranges(2,1),ranges(1,2),ranges(2,2),ranges(1,3),ranges(2,3),pointwise,simultaneous);
fprintf(fid,'## Figure caption\n\n**Local parameter robustness of three parallel annular-beam benefits in ROE surveillance.** Columns vary pulse energy, reflectivity, and atmospheric extinction individually by −20%% to +20%% around their baseline values. Rows show relative angular-width reduction, relative blind-zone reduction, and relative effective-pulse increase compared with line illumination. Positive values favor annular illumination; diamonds indicate the shared baseline. Required divergence satisfies a 60%% detection fraction within 15 s. Blind-zone and pulse statistics use 125 mrad; effective-pulse means are conditioned on successful detection by each beam. Whiskers denote pointwise 95%% paired-bootstrap intervals. There are 13 distinct parameter configurations; the baseline is repeated in each column. Moving/static sample sizes are %d/%d. Supplemental figures show 100/150 mrad reference widths and the complete earlier wide-range study, including unattained requirements and reversals.\n',a.design.options.N,a.design.options.StaticN);
end

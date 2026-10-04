# 论文 Results 参数敏感度与联合扰动分析

本目录提供一套新的六子图方案。核心问题是：在参数发生变化时，环形光束能否仍以较小角宽满足同一探测要求，并在多个参数同时变化时保持探测和覆盖收益？

## 方法

1. **同要求的角宽响应，子图 a–d。**复用 `para_sensative_N3000_range2` 的完整逐目标原始结果（每配置 N=3000），比较 15 s 内检出比例达到 60% 所需的最小角宽。对 E、ρ、v、α 分别使用原有全部参数水平，包括未达标的弱回波端点。直接寻找原始经验曲线第一次达到门限的位置，仅在相邻两个角宽点之间插值，不做平滑、累计最大包络或人为单调修正。未达标显示 NA；CSV 同时给出第一个实际模拟达标角宽，供判断插值精度。
2. **多参数联合扰动，子图 e–f。**在 E、ρ、v、α 各为基准值 ±20% 的四维区间内生成 24 组拉丁超立方组合，另加一个基准。所有组合都计算，没有按输出筛选组合。固定 w_D=100 mrad，比较 15 s 内检出比例差 ΔP=P_R−P_L，以及一个完整扫描周期的静态盲区差 ΔB=B_L−B_R。两者为正均表示环形更优。每组合默认 N=512 个动态目标；静态目标先计算512个，再对所有组合统一追加1536个，总计2048个。追加计划在静态扩充计算前固定，不按某个组合是否显著决定停止。
3. **配对与不确定性。**同一实验内跨光束、参数水平和角宽保留共同目标与相位。OAT 对目标索引作 1000 次配对 bootstrap，给出点态 95%区间和每种光束的重采样可达比例。联合试验对动态和静态两个独立样本库分别重采样，库内保持光束与场景配对；使用自助法最大标准化偏差，构造同时覆盖 25 个配置、2 个差值指标的近似 95%同时区间。联合图的序号0表示额外基准，不计入24组拉丁超立方组合。区间只反映目标采样误差，不含模型、插值或输入分布不确定性。
4. **稳健性补查。**同时导出 15 s 下 50%、60%、70%、80% 四种要求的全部结果；不能因某项要求未达到就删去该结果。导出原始绝对检出概率和静态盲区的物理方向检查结果。

选择 15 s/60% 是沿用已有分析中的一组任务要求；100 mrad 与论文已有盲区例子一致。±20% 是明确指定的工程扰动区间，不是从实测统计得到的参数分布。24组组合不是整个连续参数区间的穷举，也不能用其有利比例估计现实中的优势概率。

本方案补充了共同变化的参数组合，但不是 Sobol 全局敏感度分解；不估计参数方差贡献，也不把单因素曲线斜率当作可跨参数直接比较的敏感度排名。联合试验没有重新估计每个组合的最小角宽，所以其作用是验证固定角宽下的检测与覆盖收益，不能称为联合最小角宽优势证明。

## 物理模型与论文对应

- 来源：`D:/lzx/小论文/李宗禧 20260915 面向低空小目标探测的低角域需求环形光束激光雷达-单栏.docx`。
- 原稿 SHA256：`E031C236DA57EBD47DE80B6AFCA5649A981D41262B52CB31E3D6A3519AE29015`。
- 分析对应 ROE 模式，作为 Results 的最后一个小节，建议题名 **Parameter Sensitivity and Joint-Perturbation Robustness**。
- 直接调用现有 `MC_line_snr`、`MC_ring_snr` 与 `append_and_check_photon_window`，保持原有扫描轨迹、入场方向和相位约定、有限目标照明与加权 SNR=2 判据。没有使用旧独立模块中 0.8 mrad 的参数；本方案小角宽为 1 mrad，与当前论文及 N3000 数据一致。
- 动态目标从 R=2000 m 的半球边界进入；静态点按半球体积均匀采样。动态计算在15 s或出界时结束，有限时长与所定义的限时检出事件完全一致。静态计算覆盖一个完整扫描周期。
- 大气始终满足 α/β=50 sr，因此 α 扫描表示固定激光雷达比下消光与后向散射的联合变化。
- 核心是模型中的期望光子数与最优加权 SNR 判据；没有新增显式光子噪声抽样或受控虚警率实验。
- 本节只验证 ROE。不能直接将结果推广为 ROI-RA、ROI-SP 或所有接收机/扫描硬件参数的稳健性验证。

期望物理趋势：固定角宽时 E、ρ 增加提高回波；在本参数范围内 α、β 增加使检测变差。相同探测要求下的最小角宽相应下降/上升。在当前速度扫描中驻留机会减少导致所需角宽增加，但不把速度对所有指标的单调性当作普遍定理。环形额外相遇与能量稀释存在竞争，弱回波端点可以出现优势消失或反转，不能修改数据来消除这些现象。

具体地，固定几何与距离时，当前代码有 `N_s ∝ E rho exp(-2 alpha r)`，单脉冲对累计 SNR 平方的贡献为 `N_s^2/(N_s+N_bs+N_bg+N_dark)`。改变能量时，信号与大气散射项分别随 E 线性增长，其贡献可写成 `E^2 a^2/(E b+c)`，对 E 严格递增；反射率增加也使该贡献递增。固定 `beta=alpha/50` 时，散射积分核包含 `alpha exp(-2 alpha r)`；本扫描中 `2 alpha r < 1`，该噪声项随 alpha 增加，同时目标回波下降。因此共同目标下 E、rho、alpha 的上述方向有模型内的直接依据。速度会同时改变相遇与距离，不能使用同一种单调性证明。

当两个光束已接近几何相遇所能达到的表现时，继续增加能量对最小角宽的改善趋于饱和，二者曲线接近是合理现象。不能要求每一个相对优势指标随参数单调变化，也不能仅因一条曲线不够平滑就修改统计结果。

参考：Saltelli and Annoni, *How to avoid a perfunctory sensitivity analysis*, Environmental Modelling & Software 25 (2010), 1508–1517, DOI 10.1016/j.envsoft.2010.04.012。该文讨论单因素分析的局限，是本方案增加联合扰动检验的依据，而不是采用了该文的某一幅图或完成了全局敏感度分解。[JRC 原始记录](https://publications.jrc.ec.europa.eu/repository/handle/JRC49564)

## 复现

MATLAB R2025a；原项目需在路径中，`readNPY` 可用，扫描轨迹位于原始 manifest 指定的 `G:/BeamVEC_NEW`。默认使用4个线程；无并行工具箱时显式设置 `UseParallel=false`。运行已有结果重绘不需要轨迹盘或原始仿真数据。

```matlab
addpath('D:/lzx/MatlabCode/codexagent_MC')
addpath('D:/lzx/MatlabCode/codexagent_MC/paper_sensitivity')

% 复现当前小批量联合验证；完成的配置会从文件中恢复
results = run_paper_sensitivity;

% 只重画已保存的图，不重新仿真或读取原始数据
fig = plot_paper_sensitivity(Visible="on");

% 联合验证追加到更大样本量，输出到单独的 joint_N3000 目录
% run_paper_sensitivity(JointN=3000);

% 无并行工具箱
% run_paper_sensitivity(UseParallel=false);

% 单独验证统计函数
runtests('D:/lzx/MatlabCode/codexagent_MC/paper_sensitivity/tests');
```

若改变 JointScenarios、Seed、RelativeHalfRange 或源代码，应指定新的联合试验输出目录。缓存严格检查配置及核心源文件，拒绝混用不同实验。默认生成器仅在所有目标完成后提交一个配置的结果。

## 文件

- `run_paper_sensitivity.m`：完整流程入口。
- `run_joint_sensitivity.m`：拉丁超立方设计与可恢复的联合仿真。
- `extend_joint_static.m`：保持全部动态结果及原静态样本，统一追加静态样本。
- `analyze_paper_sensitivity.m`：逐目标读入、配对区间、联合同时区间及物理方向核对。
- `sensitivity_min_width.m`：不施加单调包络的首次相邻交点插值。
- `plot_paper_sensitivity.m`：六子图导出 PNG、矢量 PDF、MATLAB FIG。
- `results/oat_minimum_width.csv`：最小角宽、实际网格达标点、区间及可达比例。
- `results/secondary_requirements.csv`：替代概率要求，包含所有未达标结果。
- `results/joint_contrasts.csv`：全部联合参数、配对效应及点态/同时区间。
- `results/physical_checks.json`：未经单调修正的物理方向检查。
- `results/paper_sensitivity.mat`：自包含绘图数据。
- `figures_final/paper_parameter_sensitivity.png/.pdf/.fig`：最终主图。
- `manuscript_final.md`：完成运行后整理的定量结果、图注及中英文论文段落。
- `verify_paper_sensitivity.m` 与 `results/verification.json`：原始结果重建、拉丁超立方分层、样本保留、完整性和产物检查。

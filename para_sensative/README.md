# 参数敏感度模块

已有原始数据时，在工程目录的 MATLAB 命令行依次运行：

```matlab
para_predraw
para_draw
```

三个入口分别负责：`para_data` 生成原始数据、`para_predraw` 预处理、
`para_draw` 直接读取预处理 MAT 绘图。后续只调整图形样式时运行 `para_draw` 即可。

当前目录已有数据的每个光束、参数值、发散角、动态/静态模式均使用 N=200；发散角为
`5:5:400` mrad，小宽度固定为 1 mrad。动态和静态使用各自的一组公共随机样本，
所有参数值、光束和发散角之间配对复用，以减少 Monte Carlo 抖动。

用户已将 `para_data.m` 的默认 N 改为 1000，这不影响处理或绘制现有 N=200 数据。
重新生成 N=1000 时请使用新输出目录，不能直接混入现有 N=200 目录。

| 参数 | 文章基准 | 默认五个取值 |
|---|---:|---|
| 单脉冲能量 | 150 μJ | 75、112.5、150、187.5、225 μJ |
| 反射率 | 0.5 | 0.25、0.375、0.5、0.625、0.75 |
| 目标速度 | 30 m/s | 15、22.5、30、37.5、45 m/s |
| α | 1.5×10⁻⁵ m⁻¹ | 0.75、1.125、1.5、1.875、2.25（×10⁻⁵ m⁻¹） |

一次仅改变一个参数，始终 β=α/50。每列横轴以文章基准为中心，使用实际物理值。

## 三行的统计定义

1. **所需发散角**：默认三个共同要求为 `(15 s, 70%)`、`(30 s, 80%)`、
   `(60 s, 90%)`。概率分母是全部 N 次试验，时间为首次成功窗口末脉冲的告警时间。
   完全沿用 `equal_wD.m`/`minimum_thetaD.m` 的 `cummax` 单调包络和首次线性插值，
   κ=环光最小发散角/线光最小发散角。两束中任一不能达到要求时为 NaN，不外推。
   这三组要求在工程已有 `equal_wD.mat` 的 EA 基准数据中均可达到；新 N=200 数据
   的可达性仍以实际结果为准。
2. **盲区**：沿用 `blind_snr_EA.m` 的上半球均匀体积静止目标、一个完整扫描周期、
   SNR 判断与 `variable` 窗口。B=1−成功数/N，不能用运动目标漏检率代替。
   选择环光盲区严格较小之后的第一个非零盲区交点（线性插值）。两根柱分别对应
   线光、环光的这个**共同发散角，因此等高**。两条右轴曲线为交点盲区率，以及
   交点以下采样发散角上 `max((B_L−B_R)/B_L)`。分母为零不参与；无交点为 NaN。
   全部候选交点、最大降幅位置均保存。因为采用静止目标定义，速度列的盲区结果相同。
3. **脉冲利用**：柱图默认在共同的 **100 mrad** 下统计成功样本平均首次探测距离。
   距离按原 `advantage_data.m` 的 0.5 m 距离分箱口径记录，由目标位置和告警时间
   直接计算，避免变能量/变大气时误用原固定参数反散射反演表。
   正且有限的信号光子数计为有效脉冲；连续有效脉冲组成一簇，算法与
   `advantage_data.m` 相同。两条右轴曲线分别是各发散角下条件均值之比的最大值，
   即 `max(meanPulses_R/meanPulses_L)` 和 `max(meanClusters_R/meanClusters_L)`。
   两个最大值的位置可不同，都保存对应发散角。默认不套用 `advantage_draw.m`
   中人为的线光脉冲数 0.6 显示修正；如需复现该显示口径，可显式设置为 0.6。

动态仿真直接调用现有 `MC_line_snr.m`、`MC_ring_snr.m`，使用 `MC_EA.m` 的
`fixed` 窗口、`MC_func.m` 的轨迹加载/边界圆/初始相位范围，以及 `init_UAV.m`
的目标采样规则。给两个内核增加了可选第 15 参数 `[E,rho,alpha,beta]`，其默认值
与此前完全一致，缓存键包含全部物理参数。没有替换高斯照明或 SNR 算法。

## 可调整的调用

```matlab
cfg = para_data(DryRun=true);       % 只检查配置，不生成数据、不启动仿真
para_data(N=1000,OutputDirectory="para_sensative_N1000");
para_predraw(DataDirectory="para_sensative_N1000");
para_draw(DataDirectory="para_sensative_N1000");
para_predraw(Requirements=[15 .7;30 .8;60 .9],RangeWidthMrad=150);
para_draw                         % 使用上一步保存的指标和标签
para_predraw(RangeAtCrossing=true); % 距离柱改在各参数的盲区交点处插值
para_predraw(LinePulseDisplayFactor=0.6);
para_draw(Export=false);           % 只显示图，不导出 PNG/FIG
para_predraw(RebuildCache=true);   % 明确要求重建原始统计缓存
```

`Factors` 控制相对于基准的取值，默认 `[0.5 0.75 1 1.25 1.5]`。
`BeamPathDirectory` 默认 `G:/BeamVEC_NEW`，需具备全部对应的 5 kHz fast `.npy`
轨迹文件及 `readNPY`。生成前会检查全部轨迹。支持串行与线程/进程并行；
`UseParallel=false` 可禁用并行。

## 文件与断点续跑

- `para_manifest.mat`：参数、文件索引、随机种子、运行配置及内核源代码记录。
- `raw/energy_s01/PARA_MC_1par_EA_LINE_D100d1mrad.mat`：动态原始数据示例。
- `raw/energy_s01/PARA_STATIC_MC_1par_EA_LINE_D100d1mrad.mat`：静态原始数据示例。
- `raw/baseline/`：四列重复的基准值、速度列的静态结果共享，避免重复仿真。
- `para_summary.mat/.csv`：各发散角曲线、三类指标、选择的交点与最大值位置。
- `para_plot_data.mat`：预处理生成的绘图输入；含全部绘图数值、参数标签和分析设置。
- `para_statistics_cache.mat`：每个原始文件的紧凑统计缓存，包括成功告警时间和脉冲数等均值。
- `para_sensitivity.png/.fig`：完整 3×4 图，后两行各有一个左轴和两个独立右轴。

原始 MAT 中保留 `MC_EA.m` 的 N×1 `detect_L/R`（logical）、`first_time_L/R`、
`first_encounter_time_L/R`、`effective_pulses_L/R`，增加 `meta` 和
`first_detection_range_L/R`。只研究线/环两束，不额外计算点光。

重新运行相同配置会跳过完整文件；不兼容配置或内核修改会拒绝混用，应换输出目录。
`MaxStepIndex` 仅用于开发时截断测试，默认 Inf；绘图会拒绝截断数据。

## 绘图速度与进度

`para_predraw` 负责原始数据/缓存读取、三类指标计算和处理结果保存，不打开图片。
`para_draw` 只读取一个 `para_plot_data.mat`，随后绘制 3×4 图并导出 PNG/FIG。
它不扫描原始目录，不检查原始文件时间戳，不读取清单或统计缓存，也不会自动
调用预处理。各阶段会在命令行输出提示。

第一次运行 `para_predraw` 需要处理原始数据：每个文件只读取一次；静态盲区文件只读取探测标志和
元数据，动态文件统计首次成功窗口。读取时约每 2 秒报告进度，每 100 个文件
保存一次统计缓存。因此首次统计被中断后，再运行也能复用已处理的部分。

之后再次运行 `para_predraw` 会检查生成配置、原始文件大小与修改时间，复用紧凑缓存。
修改性能要求或距离取值不需要重新读取全部光子窗口；某个原始文件
发生变化时，只重新读取该文件。若外部工具修改内容但人为保留大小和时间戳，
请使用 `para_predraw(RebuildCache=true)` 强制重建。

生成了新数据、改变性能要求或统计口径后，需要主动重跑 `para_predraw` 更新
绘图 MAT，再运行 `para_draw`。纯绘图始终使用上次保存的数值和对应标签，不会
自动判断原始数据是否更新。将 `para_plot_data.mat` 单独复制到其他目录，也能用
`para_draw(DataDirectory="该目录")` 绘图，不需要附带原始数据。

`Export=false` 跳过 PNG/FIG 导出，适合快速查看和调整图形。如果界面已出图但
程序仍运行，可根据 `[3/3]` 提示判断是否处于图片导出阶段。

2026-09-23 拆分前，对当前 N=200、4800 个原始 MAT（约 3.15 GiB）的实际测量：
首次统计与绘图 156.1 秒，缓存重绘 4.3 秒，缓存重绘并导出 PNG/FIG 13.4 秒。
缓存文件约 3.4 MB。修改前后在构造数据上逐项比对结果一致，实际数据首次统计
与缓存统计也完全一致；另已验证仅修改一个原始文件时只重读该文件。

2026-09-24 拆分后，复用现有缓存完成预处理约 1.3 秒，生成的绘图 MAT 为
5158 字节。已验证三类指标及各发散角曲线与拆分前完全一致，并在只包含绘图 MAT、
没有原始文件/清单/缓存的目录完成独立绘图和 PNG/FIG 导出验证。实际绘图耗时
还包含 MATLAB 图形初始化与文件导出开销。

## 验证记录

已通过 33 项现有回归测试，并检查默认/显式物理参数一致性、能量和反射率缩放、
缓存切换、串行数据生成、断点续跑和 N=2 的完整动态轨迹/静态扫描周期及线程池。
无探测、零盲区、无交点保持 NaN 的处理及最终 `print` 图片导出也已在独立
MATLAB 进程中验证；布局预览见 `verification_synthetic/final_layout.png`。
`verification_smoke/`、`verification_full_cycle/` 为小规模开发验证数据。
`verification_synthetic/` 是检验 κ、交点、比值与图形布局的**人工构造数据**，
不是仿真结果，也不能用于论文。初次交付时未运行默认 N=200 的全部参数仿真；
后续已使用用户生成的完整 N=200 数据完成统计缓存和 PNG/FIG 导出验证。

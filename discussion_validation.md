# Discussion 仿真验证记录

日期：2026-09-08。MATLAB：R2025a。
共享代码 SHA-256：b2fe836b3ea9ff6f74ff925929fc74ca6f4a76a788acdbedb7e47c3f20d52f97。

## 实施范围

新增 5 个场景入口 tmp_1.m 至 tmp_5.m、22 个 +discussion 公共模块和一个测试类；未修改原有仿真或论文原稿。当前实现以 ROE 为完整系统验证场景，方向方案还包含局部相交模型。

## 执行记录

| 场景 | 完整系统样本数 | 其它样本数 | 输出图数 |
|---|---:|---|---:|
| tmp_1 | 128 / 工况 | 10 因素各两端点、9 个交互点、3 条双光束宽度响应 | 1 |
| tmp_2 | 128 / 工况 | 局部相交 1024；航向 8 点×2 速度；线光场旋转 6 点 | 2 |
| tmp_3 | 1024 / 工况 | 4 个发散角，3 种决策规则，3 个闭合性工况 | 2 |
| tmp_4 | 128 / 工况 | 接收器、背景、抖动、缺口、相位、积分窗对照 | 1 |
| tmp_5 | 训练 512；验证 1024 | 10 个宽度，6 个角域预算，7×10 个时间/性能要求 | 1 |

这些是实际运行的预览数据。尚未执行全部 paper 模式默认 N=3000、5 mrad 网格的正式扫描。关键输出的配置和完整样本数保存在 MAT 文件中。

## 一致性与单元测试

直接对照当前 MC_line_snr / MC_ring_snr 的 64 次配对试验，探测判定及首次探测时间全部相同，最大时间误差 0 s。

tests/discussionSimulationTest.m 共 8 项，最终运行全部通过：

1. originalLineParity
2. originalRingParity
3. annularRotationInvariance
4. accumulationDominatesSingle
5. nonmonotonicFrontier
6. unreachableRemainsMissing
7. zeroDiscordanceStillHasUncertainty
8. photonSolidAngleNormalization

测试中的原内核一致性依赖当前 G:/BeamVEC_NEW 扫描库。公共包使用 PathFixture 加入仓库路径，不污染永久 MATLAB 搜索路径。

## 静态检查与依赖

28 个新增 MATLAB 文件由 Code Analyzer 检查，无报告问题。requiredFilesAndProducts 识别 MATLAB 和 Parallel Computing Toolbox；后者仅在 cfg.useParallel=true 时使用。原扫描库读取需要用户环境已有的 readNPY。CI、随机数、积分和统计表不需要 Statistics and Machine Learning Toolbox。

## 图件质量

已生成并检查 7 套图，每套包含矢量 PDF、600 dpi PNG 和 FIG。图件采用统一字号、颜色、线型、坐标标注与白色背景；修正了半环消融横坐标标签的拥挤以及柱状图图例遮挡。

复现文件：discussion_results/preview/tmp_1/ 至 tmp_5/。共享源码和配置被保存用于缓存核对；外部扫描库本身未被打包，若替换扫描库需使用新的缓存目录。

## 主要推断边界

- 数据是随机轨迹下期望 SNR 的可行比例，不是已校准虚警率的实际检测器 ROC。
- 配对区间为点态或指定单一有限网格比较的同时区间，不能对所有扫描因素自动宣称全局显著。
- 有限网格最小可行宽度不等于连续最优值；未达标与左边界截断均单独标记。
- 保留了强背景、低回波和过大发散角下的负结果。


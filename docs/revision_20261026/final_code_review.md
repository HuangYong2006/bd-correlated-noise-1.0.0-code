# 最终只读代码审查（2026-10-26）

审查结论：未发现阻碍当前科学修订交付的 Critical 或 Important 问题。该结论针对本工作区当时的未提交改动和新生成数据，不表示批准合并、发布或投稿。

审查依据是实际读取的代码、`git diff`、105 个新原始数据单元、派生表与压缩包；未沿用先前审稿结论。未修改源码、数据或 Git，只新增本报告。

| 等级 | 结果 |
|---|---|
| Critical | 无可核验发现。 |
| Important | 无可核验发现。 |
| Minor | 未列出需要本轮修复的交付问题。 |

实际追踪的数据流与检查如下：

- `code/bd_core.py:460`：每个被滤波的行映射到一个 run/site 的完整时间记录，`eta[rows // L, :, rows % L]` 保持 `(run,time,site)`。`code/bd_core.py:354`、`:370`、`:518` 保存、读取并标记 `spectral-run-time-site-v2`，缺失标记的旧存档读为 `unversioned`。
- `code/diagnostics.py:163`、`:196`：Hellinger 的外层重采样单位是完整 run，重复计算包含信号与参考 split-null 的完整 excess；条件分割散布另存。`code/diagnostics.py:318` 的训练、阈值校准和 null run 互不重叠；`:386` 将同一 run 的晚期特征先平均，并明确森林误差为参考集合重分割散布。`code/analyze.py:169` 与 `code/feature_ablation.py:304` 没有将多个快照展平为独立样本。
- `code/revision_analysis.py:28`、`:40`、`:63`、`:99`：先对每个 run 的晚期 W² 平均，再对 ensemble 取 RMS；尺寸斜率独立重采样各尺寸的 run，相邻尺寸斜率沿用同一批抽样；谱数据必须通过 corrected-layout guard。新分析将零参数 excess 明确定义为零，未从自比较误差构造物理结论。
- `code/export_revision_tables.py:10`、`:53` 从当前 summary/noise 数据生成正文数字；`code/build_overleaf.py:43` 按明确文件清单打包。当前压缩包没有依赖旧谱分析或旧图。

在 `D:/TOOLS/anaconda3/python.exe -B` 中执行了不落盘的独立核验，全部断言通过：

1. 逐一加载 105 个新 raw，核对 SHA-256、参数、64-run 数组形状、最终 profile 时间、宽度时间单调性、有限非负 W²、BD 更新规则和谱类型。每个 campaign 正好是声明的 35 个单元。
2. 重新计算全部单元的 W、p_flat、晚期窗口比，以及它们的 500 次 run-bootstrap SD；与 summary 数值完全一致。
3. 重新计算全部 21 组 500 次 bootstrap 尺寸拟合和相邻尺寸斜率；与 summary 完全一致。
4. 重新计算全部单元的 Hellinger 中心值、split-null 均值与散布；重新计算每个 generator 在 θ=0.25、L=32/512 的共 6 个外层 bootstrap 误差；与存档值完全一致。未独立重算其余 84 个非零参数单元的完整 200-draw Hellinger 误差。
5. 两个 CSV 与 summary 一致；在禁止导出写入的情况下执行数字生成部分，所得 `revision_results.tex` 与现有文件逐字一致。
6. zip 恰含 13 个声明文件；每个成员与 `manuscript/overleaf/` 完全一致，其正文、参考文献、数字宏、中文摘要、PDF/SVG 图也与当前源文件完全一致。figure_sources 的 summary/noise 哈希、noise_validation 的三个源文件哈希均匹配；八个生产噪声验证单元通过，故障对照拒绝结果也匹配。
7. 另直接生成 θ=0.25、L=17、T=16384、3 runs 的两种谱场，改变 row/run chunk 后样本逐元素一致。测得 power-law 的时间/空间 lag-1 相关为 `0.240145/-0.000694`，fGn 为 `0.268510/-0.001485`。
8. Hellinger 的独立二项式例子得到外层 SD `0.08201861`，与条件 signal split SD `0.10867490` 分离；森林重复快照不改变 run 数或结果，三组参考 run 的索引交集为空。

本次未重新运行会写入临时文件的完整 pytest 套件，也未重新编译或渲染 LaTeX；48 项测试通过属于主执行流程的验证记录，不冒充本次独立审查所得。此次对包的验证是源文件及字节一致性验证。

超出本轮软件审查范围的事项：进一步放大 L/T、精细 θ 网格、估计 θ_c、证明渐近标度、重跑历史 campaign、核验历史作者实现、公开数据存储与作者投稿信息。缺少这些研究或发布环节不在本报告中列为代码缺陷。当前正文已将 late-time/finite-size 结果与渐近物理断言区分。

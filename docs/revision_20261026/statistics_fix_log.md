# Independent-run statistics correction

2026-10-26。Task 2 完成。修改范围为 `code/diagnostics.py`、`code/analyze.py`、`code/feature_ablation.py`，新增 `tests/test_diagnostics_runs.py`。未改主文、生产噪声代码、绘图代码、冻结评审或原始 campaign，未提交 Git。

## 实现

- 森林接收 `(runs, features)` 或 `(runs, snapshots, features)`。后者先按 run 平均末期 snapshots，一个独立 run 只有一个 vector。两个调用入口先按 runs 匹配数量，再交给同一实现，消除 flatten snapshots 的伪重复。
- 参考 runs 分为 60% training、20% threshold-calibration、其余 independent-null。特征选择和标准化只看 training；null FPR 直接在未参与训练及阈值选择的 null runs 上测量。
- `dPhi` 是有符号的目标与 held-out null 比例差，`dPhi_raw` 保留相同值。`z_anom` 是标准化 score effect，不是 significance z。森林 `*_err` 是 reference-run resplit dispersion，不是完整 run-bootstrap sampling standard error；未观察到散布时返回 NaN，避免虚假零误差。
- 返回实际独立样本数、train/calibration/null 数、calibration resolution、实际 null FPR 和条件二项区间。区间是各次固定森林下 95% Clopper-Pearson 区间的包络，不是平均森林估计量的无条件置信区间。64 runs 时为 38/12/14，校准分辨率 1/13，不能宣称精准 2% FPR。0/14 null events 的条件区间上界约 0.2316。
- null 或 score 方差为零时，标准化量为 NaN，不用 epsilon 产生巨大 z。样本不足以使三组各至少有两个 runs 时，返回完整 NaN 结果及明确 status。θ=0 自参考时仅对独立 null 组打分，不混入训练 runs。
- 两类 Hellinger 函数的 `dH_err` 改为外层按完整独立 runs 有放回重采样、重算完整 `mean(Hsignal)-mean(Hnull)` 的 SD。同一外层样本内 signal 和 null 共享同一参考 resample。内部使用共同随机 partitions 降低额外 Monte Carlo 噪声。旧 conditional signal split scatter 另存。
- onset 仅保留描述用途。文档和 summary 明示逐点 Gaussian 近似忽略共享参考的跨 θ 协方差；未知逐点误差时 onset uncertainty 为 NaN，移除固定 `zH_err=1`。没有添加新物理 FSS 框架。旧分析 alpha error 也明确标为 log-fit residual error，不能冒充 run bootstrap；根的新主分析独立处理该量。

## 公共接口及元信息

```python
discrete_hellinger_excess(ref_pmf, th_pmf,
                          n_splits=120, seed=0, n_boot=200)
```

既有 `dH`、`dH_err`、`Hsig`、`Hnull_mean`、`Hnull_std`、`z`、`H_full` 保留。新增 `Hsig_split_std`、`dH_split_std`、`n_ref_runs`、`n_theta_runs`、`n_half`、`n_bootstrap`、`error_method`。`n_boot=0` 仅算中心量与 null floor，误差为 NaN。`z` 是相对于条件 null scatter 的描述量，零 null 方差时为 NaN。

自参考 baseline 应传同一个数组对象，例如 `discrete_hellinger_excess(ref, ref, ...)`，使外层 bootstrap 保持同一 ensemble 身份。原分析新增 `--hell-boot`，默认 200。主 pmf 和 terrace pmf 使用外层 bootstrap；已知不适合物理解释的 histogram controls 只保留中心量和条件 null，误差明确未估计。summary 保存语义，CSV 保存实际独立样本数。

## RED / GREEN

解释器为 `D:/TOOLS/anaconda3/python.exe`。

1. 首先新增六个回归测试，原实现六个均 RED。一个可独立枚举的八 run pmf 例子，其外层 bootstrap 理论 SD 为 `0.07676388065406728`，旧 conditional split SD 为 `0.12247043492572626`。其他测试验证零 null 方差的 infinity、没有独立 null、常数 score 仍返回有定义 z、没有独立 run 元信息，以及小样本空 calibration 错误。
2. 第一轮修复后五项 GREEN。oracle 测试的 200 次 bootstrap Monte Carlo 波动略超预设容差，因此该核验显式使用 1000 次，得到 SD 约 `0.08202`，落入预设容差。没有放宽容差，生产默认仍为 200。
3. 再增加自参考森林 baseline 与未知 onset error 两项，修复前均 RED，分别得到 `dPhi=-0.01111` 和未知误差被报告成 `0.0`；修复后八项 GREEN。
4. 增加真实 analyze 与 feature_ablation CLI 集成检查，使用临时合成的 12 runs × 4 snapshots，两个入口均记录 12 independent runs，而非 48 snapshots，并记录 3 independent-null runs。统计新测试共九项。
5. 本任务最后一次完整项目测试为 **48 passed**，无 warning。命令为 `python -m pytest -q --tb=short --basetemp <workspace>/docs/revision_20261026/statistics_pytest_20261026_3`。此前一次遇到系统 pytest 临时目录不可写，改用工作区；另一次遇到根正在实现的三个 revision-analysis RED 测试，根完成后完整套件通过。对本任务文件的 `git diff --check` 通过。

64 runs、49 bins、120 splits、200 outer bootstraps 的单次 pmf 检查约 0.55 秒，只是函数级速度检查，没有启动 campaign。

## 统计边界

森林是补充敏感性检查，跨 reference split 散布没有包含目标 ensemble 的全部采样误差，不能解释为完整置信区间。小 calibration/null set 无法高精度标定 2% 事件率。旧 onset 没有跨 θ 协方差模型，其物理解释已移除。Hellinger SD 是有限 outer bootstrap 与有限 inner partitions 的数值估计，应报告实际设置。

上一轮工具在持久化本日志及清理临时测试目录前，因自动审批服务额度限制未执行；本次继续后优先写入日志。没有通过其他途径绕过拒绝。

继续后的复核再次得到九个统计测试全部通过，两个 CLI 的 `--help` 正常。还移除了旧 onset caller 用中位数填补未知误差的步骤，并把 `--if-boot` 帮助明确成 reference-run resplits。

对根代理的新 `revision_analysis.py` 仅做只读检查，未发现本次范围内的具体实现 bug。额外以各 size 中心值恰好落在 α=0.6 直线、但 run 内有散布的数据检验，返回 `alpha=0.6000000000000008`、`alpha_err=0.08097335406193323`，确认误差来自 runs 而不是拟合残差。`require_corrected_layout` 接受 v2 谱布局并拒绝 unversioned powerlaw 数据。本检查不等于重跑其完整 campaign 分析。

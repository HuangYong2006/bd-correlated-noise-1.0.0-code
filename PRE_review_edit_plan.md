# PRE 投稿前审阅与修改任务单

## 0. 评估范围

本文基于以下完整论文 PDF 进行预投稿审阅：

- `C:/Users/MSI/Downloads/overleaf_bd_paper__4_.pdf`
- 论文题目：`Generator bias in long-range-correlated growth simulations and the determination of the critical temporal-correlation exponent of ballistic deposition`
- 主题：长时间相关噪声、FFGN、fractional Gaussian noise、ballistic deposition、KPZ 普适性、有限尺寸标度和 onset detection。

本任务单面向后续论文修改 agent。修改时只处理本文相关内容，不进行无关重构。

## 1. 投稿判断

### 结论

**可以考虑投稿 Physical Review E，但不建议当前版本直接投稿。**

本文具有 PRE 潜力，核心价值是指出：标准 FFGN 实现可能在名义参数 `theta -> 0` 时仍保留显著时间相关性，从而把噪声发生器偏差和有限尺寸检测偏差误判为 KPZ 普适类变化。

当前主要风险不是题目不适合 PRE，而是：

1. 对既有文献的概括有些过强；
2. 对 `theta_eff`、FFGN 与 exact fGn 的关系需要更严谨限定；
3. Isolation Forest 与 Hellinger diagnostic 的“独立性”表述过强；
4. 当前有限尺寸范围不足以测定真正的 `theta_c`；
5. 统计校准和已有文献实现的核查还应更完整；
6. 论文仍有占位符、排版空格和参考文献格式问题。

## 2. 修改目标

修改后应满足以下成功标准：

- 明确本文测得的是“发生器偏差和 onset-detection 方法学问题”，而不是已经确定了 KPZ 的真实临界点。
- 将结论限定在本文直接研究的 `1+1` 维 parallel-update ballistic-deposition 模型、所用发生器和可访问尺寸范围内。
- 把 FFGN 的物理参数、等效参数和 universality-class 参数严格区分。
- 让所有主要结论都能由图表、误差条和有限尺寸分析直接支持。
- 删除所有占位符和投稿前排版问题。

## 3. 必须修改项

### P0-1. 降低过强的文献概括

全文搜索并审查以下类型的表述：

- `used by every one of these studies`
- `every published simulation`
- `the reported thresholds are properties of the FFGN noise`
- `what is no longer open`
- `every morphological and dynamical signature`

如果没有逐篇核对代码、参数、归一化和二值化实现，不要使用 `every`、`all`、`no longer open` 等绝对表述。

推荐改写为：

> For the `1+1`-dimensional parallel-update ballistic-deposition model studied here, the FFGN-based trends reported in the implementations that we could verify are quantitatively reproduced by the generator bias. They are not reproduced by the corrected spectral generator over the accessible system sizes.

摘要、Discussion 和 Conclusion 中都要保持这一限定。

### P0-2. 严格限定 `theta_eff`

在正文第一次定义 `theta_eff` 的位置明确写出：

- `theta_eff` 是对二值化噪声相关函数进行 shape matching 得到的诊断参数；
- 它不是 FFGN 的严格重参数化；
- FFGN 与 exact fGn 在中间 lag 和低频权重上仍可能不同；
- 因此不能直接把 `theta_eff` 等同于该模型的 RG 耦合参数。

推荐加入以下意思的句子：

> We use `theta_eff` only as a shape-matched diagnostic for the binarized correlation function. It is not assumed to be an exact reparameterization of the dynamical renormalization-group coupling or of the universality class.

同时检查全文是否把以下三者混为一谈：

1. nominal `theta`；
2. shape-matched `theta_eff`；
3. continuum/RG 理论中的 temporal-correlation exponent。

### P0-3. 明确本文没有测出 `theta_c`

摘要、Results、Discussion、Conclusion 必须统一使用以下逻辑：

- 本文没有测得真正的 `theta_c`；
- 在当前 `L <= 512` 的范围内，`theta <= 0.25` 的 crossover 没有被分辨；
- 结果说明 published onset values 不能在未校验噪声发生器和 detection floor 的情况下解释；
- 本文结果与 `theta_th = 1/6` 一致，但不能确认该理论预测。

避免把以下表述写成确定性结论：

- `the crossover lies beyond L ≈ 512`，除非给出明确拟合或理论依据；
- `should become resolvable at L ≳ 4096`，除非标明这是外推或研究计划；
- `the RG prediction survives`，建议改为 `is not contradicted within the accessible window`。

推荐统一表述：

> No crossover is resolved within the accessible window. The present data therefore do not determine `theta_c`.

### P0-4. 核查 FFGN 实现和已有文献

增加一个表格，逐篇列出 Ref. [9]-[14] 中与噪声实现有关的信息：

| Reference | Generator | Parameters | Binarization | Normalization | Verification basis |
|---|---|---|---|---|---|

至少核查：

- Eq. (10) 是否与 Mandelbrot FFGN 实现一致；
- `zeta_t` 的分布和中心化方式；
- `W_n` 的归一化；
- `a=6, B=2` 是否与被比较工作的实现一致；
- AR component 截断策略；
- 是否在零点二值化；
- 是否进行了额外的标准化或方差重标度。

无法确认的文献必须写成 `not verified from the published description`，不能归入“所有研究”。

### P0-5. 修改 diagnostic 的独立性表述

全文不要把 Hellinger distance 和 Isolation Forest 称为严格统计独立的两个诊断，除非增加正式证明。

建议将：

- `independent diagnostics`
- `use disjoint information`

改为：

- `complementary diagnostics`
- `diagnostics based on different feature channels`

补充以下分析之一，最好全部完成：

1. 两种 diagnostic 输出之间的 Pearson 或 Spearman correlation；
2. 删除 increment-related features 后重新运行 Isolation Forest；
3. 删除 terrace-related features 后重新运行 Isolation Forest；
4. 报告 feature ablation 对主要结论的影响。

### P0-6. 补充统计方法

Methods 或 Appendix 中明确写出：

- bootstrap 的重采样单位是 independent realization；
- 是否对 snapshot 或时间序列进行任何重采样；
- 每个 `(theta, L)` 单元的 realization 数量；
- Hellinger split-half null 的 split 数量和有效样本量；
- `chi^2` 拟合的自由参数和自由度；
- onset threshold 的具体定义；
- slope significance 的置信区间或检验方法；
- 多个 `theta` 同时检测时是否考虑 multiple testing；
- `sub-milli error bars` 与 bootstrap errors 的区别。

所有 `2.7 sigma`、`significant`、`well outside the errors` 等说法应尽量配套给出数值误差或置信区间。

### P0-7. 加强有限尺寸证据，或降低结论强度

优先增加以下模拟：

- `theta = 0.20, 0.25, 0.30`；
- `L > 512`，至少增加一个或两个更大尺寸；
- 使用 exact fGn filter；
- 保持 realization 数量和统计流程一致。

如果不能增加更大尺寸，则必须：

- 把结论限定为 `within L = 32...512`；
- 删除或弱化对 `L ≈ 4096` 的确定性外推；
- 在摘要中明确写 `the crossover threshold remains unresolved`。

## 4. 摘要修改建议

当前摘要信息密度过高、结论过强、数字过多。建议重写为 180-220 词左右，按以下顺序：

1. 说明 temporally correlated KPZ/BD simulations 的 threshold disagreement；
2. 说明 standard FFGN 在 nominal `theta=0` 不产生 white noise；
3. 说明该偏差如何影响 morphology 和 roughness；
4. 说明 corrected spectral generator 的结果；
5. 明确 fixed-`theta` finite-size test 没有在当前窗口内测出 `theta_c`；
6. 最后给出方法学结论：噪声发生器和 detection floor 必须先校验。

摘要中避免：

- `every published simulation`；
- `what is no longer open`；
- 把 `theta_eff` 写成真实临界参数；
- 过多列出 `0.0002 ± 0.0007`、`0.575 ± 0.012` 等次级数字。

## 5. 文章结构建议

建议把主线收缩为三条：

### 主线一：FFGN 在 `theta -> 0` 不连续

- correlation slope 可以正确，但 correlation amplitude 不正确；
- nominal `theta=0` 仍有显著 temporal correlation；
- N-convergence 证明不是简单截断误差。

### 主线二：发生器偏差改变 BD 的观测结果

- `p_flat` 在 FFGN 与 corrected generator 中趋势相反；
- roughness exponent 和 saturated width 的变化显著不同；
- 不能把 FFGN track 的趋势直接解释为 universality-class change。

### 主线三：onset detection 受到 size-dependent detection floor 影响

- fixed-significance onset 会随 `L` 漂移；
- fixed-`theta` scale dependence 是更可靠的判断；
- 当前尺寸范围不足以测出真正的 `theta_c`。

第 5.4 节的五条 methodological lessons 建议压缩成一个短 subsection，或者移入 Supplemental Material。

## 6. 语言和版式修复

### 必须修复

- `Acknowledgments` 中的占位符：
  - `The authors thank ... for discussions.`
  - `Computations were performed on ....`
- 正文中大量缺少空格的组合，例如：
  - `correlationamplitude`
  - `long-rangephase`
  - `atθ=0`
  - `everypublished`
- 参考文献中的作者名、期刊名、页码、连字符和空格格式。
- 统一 `theta`, `theta_eff`, `theta_c`, `theta_th` 的排版和文字说明。

### 建议修复

- 标题可缩短为：

  > Noise-generator bias in temporally correlated ballistic-deposition simulations

- Figure 2 和 Figure 3 的图例、误差条和符号略密，应提高字体并减少图内文字。
- Figure 4 的 profile offset、纵轴含义和 bar height 在图注中说明得更明确。
- 投稿前使用 REVTeX 4.2 的 `reprint` 或 `twocolumn` 格式重新检查版式。
- 完整算法、伪代码、参数和额外验证图建议放入 Supplemental Material。

## 7. 推荐的关键句改写

### 关于文献覆盖范围

不建议：

> The noise generator used by every one of these studies ...

建议：

> The implementations that we could verify use a Mandelbrot-type FFGN construction. We find that this construction does not approach white noise continuously at nominal `theta=0`.

### 关于 `theta_eff`

建议：

> The shape-matched `theta_eff` is used as a diagnostic for the binarized correlation function, not as an exact reparameterization of the dynamical RG coupling.

### 关于临界点

建议：

> We do not determine a crossover exponent or a critical value `theta_c` in the present size range. Instead, we show that the commonly used onset construction is confounded by generator mismatch and a system-size-dependent detection floor.

### 关于 RG 预测

建议：

> The result is not inconsistent with the RG prediction `theta_th=1/6`, but the present simulations do not provide an independent confirmation of that prediction.

## 8. 修改后的验证清单

### 科学结论

- [ ] 所有主要结论都有对应图表或表格。
- [ ] 所有超出 `L=32...512` 的表述都被标记为外推、预期或未来工作。
- [ ] `theta_eff` 没有被当作严格 RG 参数使用。
- [ ] 没有无证据地声称所有已有研究使用完全相同的 FFGN 实现。
- [ ] 没有把“没有观察到 crossover”写成“证明没有 crossover”。

### 统计与数值

- [ ] bootstrap 单位和样本数写清楚。
- [ ] null calibration 写清楚。
- [ ] onset threshold 的定义写清楚。
- [ ] slope、误差条和显著性判断一致。
- [ ] exact fGn filter 的验证结果可复现。
- [ ] feature ablation 或 diagnostic correlation 已报告，或删除“independent”措辞。

### 投稿格式

- [ ] 删除 acknowledgments 占位符。
- [ ] 修复缺少空格和不可读的 LaTeX 输出。
- [ ] 参考文献格式统一。
- [ ] 使用 PRE/REVTeX 模板重新生成最终 PDF。
- [ ] 所有图在双栏版式下仍可读。
- [ ] 代码、原始数据和分析脚本的存档链接真实有效。

## 9. 推荐投稿定位

投稿时建议把文章定位为：

> A numerical and methodological study showing that noise-generator validation and fixed-parameter finite-size tests are necessary before interpreting apparent KPZ crossover thresholds in temporally correlated ballistic deposition.

不建议把文章定位为：

> A definitive determination of the critical temporal-correlation exponent of KPZ growth.

## 10. 最终建议

完成 P0 项后再投稿 PRE。若无法增加更大尺寸，仍然可以投稿，但必须把文章明确定位为“数值方法学诊断和有限尺寸检测问题”，并在标题、摘要、Discussion、Conclusion 中统一承认 `theta_c` 尚未被确定。


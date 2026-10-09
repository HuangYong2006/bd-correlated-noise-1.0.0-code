# PRE 投稿前最新评估

## 评估对象

- 文件：`C:/Users/MSI/Downloads/overleaf_bd_paper__6_.pdf`
- 页数：17 页
- 当前稿件标题：`Noise-generator bias in temporally correlated ballistic-deposition simulations`
- 目标期刊：Physical Review E

## 总体判断

**当前版本已经基本达到 PRE 送审水平，但建议完成最后一轮小修后再提交。**

上一轮提出的主要科学风险已经得到实质处理。当前剩余问题主要是措辞强度、少数内部一致性和最终投稿信息，而不是研究主线或方法框架需要重写。

## 已经改好的内容

本版本已经完成以下重要改进：

- 增加了 Table I，对相关文献的噪声发生器、参数、二值化和归一化进行逐项核查。
- 明确说明不能声称所有文献使用完全相同的 FFGN 实现。
- 明确区分 nominal `theta`、shape-matched `theta_eff` 和 RG 理论中的 temporal-correlation exponent。
- 将 Hellinger distance 和 Isolation Forest 称为 `complementary diagnostics`，不再声称它们统计独立。
- 增加 Pearson/Spearman 相关性分析和 feature-ablation 分析。
- 增加 bootstrap、null calibration、自由度、onset 定义和显著性判断的附录说明。
- 摘要已经明确说明当前研究没有测得真正的 `theta_c`。
- 对 `L >= 4096` 的表述已经标注为外推，而非已经完成的模拟结果。
- 结论中已经加入 `within the accessible window`、`not confirmed` 等必要限定。
- 版式已经改为更接近 PRE 的双栏格式。

## 投稿前必须完成的修改

### P0-1. 修改摘要中的绝对表述

当前摘要包含：

> A theta-continuous spectral generator carrying the exact fractional-Gaussian-noise spectrum reproduces none of these features.

这句话过强，因为正文仍显示 corrected generator 在较大 `theta` 时存在较弱的 roughness 变化和有限尺寸 crossover 迹象。

建议改为：

> A theta-continuous spectral generator carrying the exact fractional-Gaussian-noise spectrum reproduces these trends only much more weakly and does not reproduce the published onset behavior within the accessible size range.

或：

> The corrected generator strongly suppresses these trends and does not reproduce the reported onset behavior for the accessible system sizes.

### P0-2. 限定 FFGN 结论的适用范围

全文检查以下表述：

- `no member of the family reaches the short-range KPZ regime at all`
- `The noise generator used throughout this literature`
- 任何使用 `all`、`every`、`throughout`、`none` 的句子。

建议将：

> no member of the family reaches the short-range KPZ regime at all

改为：

> For the FFGN construction and parameter choices implemented here, no nominal `theta >= 0` reaches the short-range KPZ regime.

Figure 1 图注中的：

> The noise generator used throughout this literature does not connect continuously to uncorrelated noise.

建议改为：

> The FFGN construction examined here does not connect continuously to uncorrelated noise.

### P0-3. 修改 Conclusion 中的 `alone`

当前结论包含：

> the reported thresholds of 0.20 and 0.25 are quantitatively reproduced by the generator bias alone

这句话可能被理解为已经排除了真实物理 crossover。当前数据更直接支持的是：在可访问尺寸范围内，发生器偏差足以产生这些 threshold-like trends。

建议改为：

> the reported thresholds of 0.20 and 0.25 can be quantitatively accounted for by the generator bias within the accessible size range.

或：

> the generator bias is sufficient to reproduce the reported threshold-like behavior within the accessible size range.

### P0-4. 解决 `96 realizations` 与 Table III 的一致性

正文中有：

> with run counts matched at 96 across every L and theta

但 Appendix Table III 显示：

- FFGN：`256` at `L=64...512`
- power-law spectral：`192` at `L=32...256`，`96` at `L=512`
- exact-fGn spectral：`192` at `L=32...256`，`96` at `L=512`

必须明确正文中的 `96` 是：

1. 所有 corrected-generator 数据的实际 realization 数；还是
2. onset drift 分析中从每个 cell 抽取的共同子样本数。

如果是第二种情况，建议正文改为：

> For the fixed-significance onset comparison, we subsampled to 96 realizations per `(theta,L)` cell across all sizes.

Appendix Table III 也应加入相同说明，避免审稿人认为正文和附录矛盾。

### P0-5. 收紧 Figure 3 的显著性表述

正文和图注中类似：

> the deviation becomes significantly positive only at `theta >= 0.30`

Appendix C 已经说明没有进行 multiple-comparison correction，因此建议改为：

> Positive slopes become systematic only for `theta >= 0.30`; individual two-sigma flags are not interpreted independently because no multiple-comparison correction is applied.

这样可以与 Appendix C 完全一致。

### P0-6. 替换最终投稿信息

第一页仍然存在占位符：

- `First A. Author`
- `Second B. Author`
- `Department, Institution, City, Country`

投稿前必须替换为真实信息，并检查：

- 通讯作者和邮箱；
- 作者单位；
- 日期是否应删除；
- 致谢内容；
- 代码和数据存档链接。

## 建议但不一定阻塞投稿的修改

### 1. 统一措辞

建议全文统一使用以下逻辑：

- `nominal theta`：发生器输入参数；
- `theta_eff`：二值化相关函数的 shape-matched diagnostic；
- `theta_c`：模型 crossover 的待测参数；
- `theta_th`：理论 RG threshold。

不要在同一段落中交替使用 `theta` 而不说明具体含义。

### 2. 继续修复 LaTeX 粘连

双栏版式已经明显改善，但仍建议全文搜索并修复类似：

- `atθ`
- `withL`
- `theθ`
- `notwhite`
- `correlationamplitude`

这些问题在 PDF 中会降低正式投稿观感。

### 3. Figure 1 图注

把：

> The noise generator used throughout this literature ...

改成仅针对本文实现的表述，避免图注的范围比正文更宽。

### 4. Appendix C 的多重比较说明

当前说明已经比较诚实。建议在 Figure 3 图注中同步提醒读者，避免单独的 `2 sigma` 圆圈被理解成经过多重比较校正的显著结果。

### 5. 第 16 页版式

Appendix C 结束后第 16 页留白较多。不是科学问题，但可以通过调整分页或让参考文献提前开始来改善最终版式。

## 不建议再做的改动

以下内容当前不建议继续扩张：

- 不要在没有新大尺寸模拟的情况下宣称已经确定 `theta_c`。
- 不要把 `theta_eff` 重新解释成严格 RG 参数。
- 不要恢复 `independent diagnostics` 的说法。
- 不要继续增加方法学列表，除非移入 Supplemental Material。
- 不要把本文扩展成对所有 KPZ、所有 FFGN 实现的普遍定理。

## 最终检查清单

### 科学结论

- [ ] 摘要中的 `reproduces none of these features` 已改弱。
- [ ] 所有 FFGN 结论都限定到本文实际验证的实现和参数。
- [ ] Conclusion 中的 `alone` 已改为 `can be accounted for` 或等价表达。
- [ ] `theta_eff` 没有被当作严格 RG 参数。
- [ ] 明确说明当前研究没有测得真正的 `theta_c`。
- [ ] 所有超出 `L=32...512` 的结论都标注为外推或未来工作。

### 统计和数值

- [ ] `96 realizations` 与 Table III 的关系已经解释。
- [ ] bootstrap 的重采样单位和数量保持一致。
- [ ] null calibration 和 onset definition 保持一致。
- [ ] Figure 3 的显著性措辞与 multiple-comparison 说明一致。
- [ ] exact-fGn filter 的结果和正文数值相互一致。

### 投稿格式

- [ ] 替换真实作者、单位、通讯作者和邮箱。
- [ ] 删除或确认日期格式。
- [ ] 删除所有占位符。
- [ ] 修复 LaTeX 文字粘连。
- [ ] 检查图注中的范围表述。
- [ ] 确认代码、数据和分析脚本的存档链接。
- [ ] 最终用 PRE/REVTeX 模板重新生成 PDF。

## 最终建议

完成以上 P0 项后，可以投稿 PRE。

当前最合适的投稿定位是：

> A numerical and methodological study showing that noise-generator validation and fixed-parameter finite-size tests are necessary before interpreting apparent KPZ crossover thresholds in temporally correlated ballistic deposition.

不建议把文章定位为：

> A definitive determination of the critical temporal-correlation exponent of KPZ growth.

最终判断：

**当前稿件已经从“需要大幅修改后再考虑投稿”进入“完成最后一轮措辞、一致性和投稿信息检查后即可投稿 PRE”的阶段。**


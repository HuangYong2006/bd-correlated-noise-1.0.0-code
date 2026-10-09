# PRE 投稿就绪度与最终排版检查（overleaf_bd_paper__8_.pdf）

## 结论

**这篇稿件的科学主线、篇幅和基本 APS/REVTeX 版式已经达到可以向 Physical Review E 投稿的水平。**

但当前 PDF 还不建议原样上传。需要先完成下面的 3 项提交前修正：

1. 删除或改正首页的空日期标记 `(Dated:)`；
2. 修正 PDF 元数据中的作者占位符 `Author list to be filled in`；
3. 补充通讯作者邮箱，并在投稿系统中核对作者与单位信息。

这些是投稿材料和排版完整性问题，不代表需要重做论文科学内容或主要模拟。

---

## 一、对当前 PDF 的核对结果

| 项目 | 当前状态 | 判断 |
|---|---|---|
| 页面数量 | 17 页 | PRE 没有对常规研究论文设置固定页数上限；篇幅偏长但有充分方法、附录和统计细节，原则上可接受 |
| 页面尺寸 | Letter，612 × 792 pt | 符合 APS/REVTeX 常见投稿尺寸 |
| 版式 | REVTeX/APS 双栏 | 符合 PRE 初投稿文本的常见格式 |
| 图表 | 图、表嵌入正文，编号和图注完整 | 未见明显越界、遮挡或不可读问题 |
| 参考文献 | [1]–[24] 均出现 | 内容完整；[24] 单独落在第 17 页，建议优化分页 |
| 首页作者 | 页面显示 `Yong Huang` | 页面显示正常，但 PDF 元数据仍是占位符 |
| 首页单位 | `College of Physics, Chengdu University of Technology, Chengdu, China` | 基本完整；投稿前确认正式英文校名 |
| 通讯作者信息 | 当前首页未见邮箱 | 投稿前应补充 `\email{...}` 或按投稿系统要求填写 |
| 日期 | 首页显示 `(Dated:)` | 不应保留空日期标记 |
| PDF 元数据 | `Author: Author list to be filled in` | 必须改为真实作者信息 |

---

## 二、必须在上传前完成（P0）

### 词数要求

- PRE 常规研究论文没有固定的全文英文词数上限；17 页（含附录和参考文献）本身不会因为词数而不合格。
- APS 对摘要通常要求保持在 250 个英文词以内。当前摘要约 **187 词**，处于安全范围内。
- 标题、正文、图注、表注和参考文献目前没有需要额外填写的总词数限制；仍应保持内容紧凑，避免重复论证。
- 投稿系统若单独限制摘要字符数，应以系统实时提示为准；若系统要求比 250 词更短，优先压缩摘要而不是删减方法和统计验证。

### P0-1. 修正作者和 PDF 元数据

在 LaTeX 源文件中检查以下位置：

```latex
\\author{Yong Huang}
\\affiliation{College of Physics, Chengdu University of Technology, Chengdu, China}
\\email{真实通讯作者邮箱}
```

同时检查 `hyperref` 或 PDF 信息设置中是否仍有类似内容：

```latex
pdftitle = {...}
pdfauthor = {Author list to be filled in}
```

应改为真实作者，例如：

```latex
\\hypersetup{
  pdftitle={Noise-generator bias in temporally correlated ballistic-deposition simulations},
  pdfauthor={Yong Huang},
  hidelinks
}
```

如果作者不止一人，`pdfauthor` 应列出全部作者；不要保留任何 `to be filled in`、`Author list` 或模板占位符。

### P0-2. 修正空日期

当前首页的 `(Dated:)` 是明显的模板残留。二选一：

- 需要显示日期：使用 `\\date{\\today}` 或明确的投稿日期；
- 不需要显示日期：按所用 REVTeX 模板的方式关闭日期，并重新编译，确认首页不再出现空的 `(Dated:)`。

不要直接用 PDF 编辑器覆盖文字；应修改 `.tex` 后重新编译。

### P0-3. 补充通讯作者邮箱

PRE 投稿系统通常要求通讯作者联系方式。建议在作者信息中加入：

```latex
\\email{your-real-email@example.com}
```

若作者和通讯作者不是同一人，按 APS 模板标明对应作者。邮箱不要使用示例地址。

---

## 三、建议在最终编译前完成（P1）

### P1-1. 合并最后一页参考文献

第 17 页目前只有参考文献 [24]，页面大部分为空白。这不会导致 PRE 退稿，但显得稿件未做最终分页整理。优先尝试：

1. 检查参考文献条目中是否有不必要的换行、过长 DOI 或 URL；
2. 使用 APS 推荐的 DOI/URL 排版方式；
3. 适当减少参考文献段前后间距或修正 bibliography 设置；
4. 重新编译，目标是让 [24] 与第 16 页合并，或至少避免单独占用一页。

不要为了强行压缩而明显缩小正文字号、图注或参考文献字号。

### P1-2. 隐藏超链接彩色边框

如果最终 PDF 仍显示 DOI、arXiv 或期刊链接的彩色边框，建议使用：

```latex
\\usepackage[hidelinks]{hyperref}
```

或在已有 `hyperref` 设置后使用：

```latex
\\hypersetup{hidelinks}
```

这是视觉清理项，不改变链接可用性。

### P1-3. 做一次源文件级文本检查

全文搜索可能的 LaTeX 粘连或缺空格，例如：

```text
atθ
withL
theθ
notwhite
correlationamplitude
exactpmf
```

PDF 的双栏抽取文本会自然地产生部分断词，因此应以源 `.tex` 和页面视觉结果为准；不要仅根据 `pdftotext` 的断词结果批量修改。

### P1-4. 检查参考文献元数据

逐条确认：作者顺序、重音符号、期刊名、卷号、文章号/页码、年份、DOI 和 arXiv 编号。特别检查较新的 2024–2026 年文献是否已经正式发表，以及预印本信息是否与期刊信息重复。

---

## 四、篇幅和内容定位建议

### 当前篇幅是否过长？

17 页（含附录和参考文献）对 PRE 常规数值/方法论文并不自动构成问题。当前篇幅主要来自：

- 发生器推导和校准；
- 两种诊断方法及其 null calibration；
- 固定参数有限尺寸检验；
- 统计协议和附录证明。

这些内容与论文的可复现性和核心方法论结论直接相关，因此不建议为追求页数而大幅删减。

如希望让主文更紧凑，可将部分实现细节和完整统计协议放入 Supplemental Material；这属于投稿策略优化，不是当前投稿的必要条件。

### 推荐投稿定位

建议把稿件定位为：

> A numerical and methodological study showing that noise-generator validation and fixed-parameter finite-size tests are necessary before interpreting apparent KPZ crossover thresholds in temporally correlated ballistic deposition.

不要把标题、摘要或 cover letter 写成“已经精确测定了临界相关指数 `theta_c`”。当前稿件已经正确限定：在可访问的 `L = 32...512` 范围内没有测得真正的 `theta_c`。

---

## 五、最终编译和上传前验收

后续 agent 按以下顺序执行：

1. 修改作者元数据、邮箱和日期；
2. 重新编译至少两次，使交叉引用、图表编号和参考文献稳定；
3. 检查首页不再出现 `(Dated:)` 或作者占位符；
4. 检查 PDF 信息中的 `Author` 为真实作者；
5. 检查第 17 页是否仍只剩 [24]；若仍存在，确认这是可接受的分页而不是引用格式异常；
6. 逐页快速查看：图注、表注、公式、附录标题、页码、双栏底部是否有截断；
7. 确认代码、数据、分析脚本的公开存档链接真实可访问；
8. 在 APS 投稿系统中核对：作者顺序、通讯作者、单位、标题、摘要、关键词和 cover letter；
9. 上传最终 PDF 前，删除所有模板占位符和内部注释。

## 最终判定

完成 P0-1 至 P0-3 后，**这篇稿件可以作为 PRE 的投稿文本上传**。P1 项目是最终排版和专业度提升，除非编译后出现新的截断、缺字或引用错误，否则不需要因为 17 页篇幅本身继续推迟投稿。

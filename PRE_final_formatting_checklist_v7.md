# PRE 第 7 版最终排版与投稿检查清单

## 适用文件

- PDF：`C:/Users/MSI/Downloads/overleaf_bd_paper__7_.pdf`
- 目标期刊：Physical Review E
- 当前判断：科学内容已达到可以投稿的程度，本清单只处理最终排版和投稿材料。

## 结论

第 7 版可以投稿 PRE。提交前建议完成以下排版修正，但这些问题不再构成科学内容上的投稿阻塞项。

## P0：提交前必须检查

### 1. 参考文献分页

当前第 17 页只有参考文献 `[24]`，页面大部分为空白。建议通过调整参考文献排版，使 `[24]` 与第 16 页合并。

可采用以下任一方式：

- 适当减小参考文献条目之间的垂直间距；
- 减小参考文献字体一个等级；
- 微调 `\bibliographystyle` 或参考文献环境的 `\itemsep`；
- 在不影响正文和图表可读性的前提下调整页面分页。

不要为了填满页面而显著缩小正文、图注或公式字体。

### 2. 隐藏参考文献链接边框

当前 PDF 的 DOI、arXiv 和期刊链接显示为青色边框。正式投稿版建议在导言区使用：

```latex
\usepackage[hidelinks]{hyperref}
```

如果模板已经加载 `hyperref`，不要重复加载，只修改选项或使用：

```latex
\hypersetup{hidelinks}
```

目标是保留链接可点击性，同时去除 PDF 中的彩色矩形边框。

### 3. 最终作者和单位信息

确认首页已经使用真实信息：

- 作者姓名；
- 作者单位；
- 通讯作者；
- 通讯邮箱；
- 基金信息，如适用；
- 致谢内容；
- 代码和数据存档链接。

不要保留以下占位符：

- `First A. Author`；
- `Second B. Author`；
- `Department, Institution, City, Country`；
- `...`；
- `....`；
- 未完成的 repository URL。

### 4. 日期设置

确认是否需要保留：

```latex
\date{September 25, 2026}
```

如果这是投稿前内部版本日期，可以删除或改为模板默认设置。不要让内部编译日期误显示为正式投稿日期，除非作者明确需要保留。

## P1：推荐修正

### 5. 检查正文中的文字粘连

从 PDF 文本抽取和视觉检查看，公式附近仍可能出现 LaTeX 字符粘连。全文搜索并检查以下类型：

- `atθ`
- `withL`
- `theθ`
- `notwhite`
- `correlationamplitude`
- `exactpmf`
- `long-rangephase`
- `fixed-θscale`
- `crossoverlength`

这些问题通常来自数学模式和普通文本之间缺少空格。修改时保持数学含义不变，例如：

```latex
at $\theta$
with $L$
the $\theta$ dependence
not white noise
```

### 6. 统一数学符号排版

全文统一检查：

- `\theta_c`：模型 crossover 参数；
- `\theta_{\rm th}`：理论阈值；
- `\theta_{\rm eff}`：shape-matched diagnostic；
- `\theta_{\rm nom}`：发生器 nominal parameter；
- `L=32\ldots512` 或 `L=32\text{--}512` 的统一写法；
- `\alpha`, `\beta`, `z`, `p_{\rm flat}` 的字体和下标。

不要在普通文本中混用 Unicode `θ` 和 LaTeX `\theta`，除非模板或作者风格明确要求。

### 7. 检查图注和表注

重点检查以下内容：

- Figure 1 图注仅针对 `FFGN construction examined here`，不要写成所有文献或所有发生器；
- Figure 3 图注保留“未进行 multiple-comparison correction”的说明；
- Table I 的 `not verified` 说明保持谨慎，不得被理解为“作者没有做该步骤”；
- Table III 明确说明 96 个 realization 是 onset comparison 和 feature-ablation 使用的共同子样本；
- 所有图中坐标轴、单位、误差条和图例在双栏宽度下可读。

### 8. 参考文献格式

检查：

- 期刊名、卷号、页码和年份格式统一；
- DOI 和 arXiv 信息是否符合 APS 格式；
- 参考文献 `[16]` 的 DOI 是否完整并可访问；
- 书籍引用 `[20]` 是否包含出版社和出版地；
- 作者姓名中的重音符号、连字符和空格是否正确；
- 参考文献是否全部在正文中被引用，正文引用是否都能在参考文献表找到。

## P2：提交前 PDF 视觉检查

重新编译后逐页检查：

- [ ] 首页标题、作者、单位和摘要没有溢出；
- [ ] 双栏之间没有文字穿栏；
- [ ] 公式编号连续且没有重叠；
- [ ] 图表没有被裁剪；
- [ ] 图中文字和坐标轴在 100% 显示比例下可读；
- [ ] 表格没有超出栏宽；
- [ ] 页码连续；
- [ ] 参考文献没有孤立的异常空白页；
- [ ] 没有彩色链接边框；
- [ ] 没有黑色方框、缺字、乱码或字体替换；
- [ ] PDF 可以复制文字和搜索关键词；
- [ ] PDF 元数据中的标题和作者正确。

## 推荐 REVTeX 设置

如果使用 REVTeX 4.2，投稿前建议确认文档类和链接设置类似：

```latex
\documentclass[
  aps,
  pre,
  reprint,
  superscriptaddress,
  nofootinbib
]{revtex4-2}

\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{amsmath,amssymb}
\usepackage[hidelinks]{hyperref}
```

具体选项应服从当前稿件的作者单位、脚注和参考文献需求。不要为了形式统一而改变已经正确的物理排版。

## 最终提交材料

除主 PDF 外，确认以下材料已经准备好：

- [ ] LaTeX 源文件；
- [ ] `.bib` 文件或完整参考文献源；
- [ ] 所有图像源文件；
- [ ] Supplemental Material，如有；
- [ ] 代码仓库或数据存档链接；
- [ ] 简短 cover letter；
- [ ] 作者贡献和利益冲突声明，如系统要求；
- [ ] 数据可用性声明，如系统要求。

## Cover letter 定位建议

Cover letter 中建议强调：

> This manuscript identifies a reproducible numerical bias in a commonly used fast fractional-Gaussian-noise construction and shows how generator validation and fixed-parameter finite-size tests affect the interpretation of apparent KPZ crossover thresholds in ballistic deposition.

同时明确：

- 本文没有声称已经精确测定 `theta_c`；
- 本文的重点是噪声发生器验证和 finite-size onset detection；
- 结论限定在所研究的 `1+1` 维 parallel-update ballistic-deposition model 和可访问尺寸范围内。

## 最终命令级检查

在最终源文件目录执行：

```text
1. 清理辅助文件并重新编译两次，确保交叉引用稳定。
2. 检查 PDF 总页数和最后一页是否只剩孤立参考文献。
3. 搜索占位符和明显粘连文本。
4. 检查 PDF 链接是否可点击但无彩色边框。
5. 在 100% 和 200% 缩放下检查所有图、表和公式。
6. 将最终 PDF 与源文件、图像、补充材料和数据链接一起归档。
```

## 最终状态

完成 P0 项后，第 7 版即可作为 PRE 投稿版本提交。剩余工作属于排版和投稿材料整理，不需要再改变论文的科学主结论。


# 时间关联弹道沉积：数据与代码

本仓库共享数值数据和复现软件。本次数据发布不含文稿、论文 PDF、TeX/BibTeX 源码或 Overleaf 包。

## 数据下载

固定版本：[`data-20261119`](https://github.com/HuangYong2006/bd-correlated-noise-1.0.0-code/releases/tag/data-20261119)。

- [当前数值数据包](https://github.com/HuangYong2006/bd-correlated-noise-1.0.0-code/releases/download/data-20261119/bd_strict_binary_data_20261117.zip)：奇偶约束、协方差扰动、联合时长/周期/模态控制、运行种子、冻结代码和分析输入。
- [历史校准数据包](https://github.com/HuangYong2006/bd-correlated-noise-1.0.0-code/releases/download/data-20261119/innovation_controls_data_20261031.zip)：保留初始失败验证、重新校准及最终验证；也嵌套于当前数据包内。
- [SHA-256 校验清单](https://github.com/HuangYong2006/bd-correlated-noise-1.0.0-code/releases/download/data-20261119/SHA256SUMS.txt)。

数据采用 NPZ、JSON 和 CSV 格式，并提供代码、参数、种子、数据来源和复现说明。包内没有论文 PDF、TeX 或 BibTeX 文件。包中包含新实验及当前比较所需的历史输入，不宣称重发了全部历史原始记录。

## 复现与引用

将数据包解压到空目录，从[数值复现说明](docs/strict_binary_20261117/REPRODUCE.md)开始。复现完整历史核验时，仅补充解压嵌套包中缺失的文件，不覆盖新代码和清单。

[严格匹配与敏感性分析汇总](data/analysis/strict_binary_20261117/summary.json)、[时长及谱控制汇总](data/analysis/controlled_20261117/summary.json)随仓库提供。逐图 CSV 位于 `figs/strict_binary_20261117/` 和 `figs/refocus_20261117/`；模拟、分析代码及测试位于 `code/`、`tests/`。参数、种子和重采样单位见复现说明。

现有 85 项测试通过；数据包内部清单及下载文件均提供 SHA-256 校验。

数据引用：Yong Huang, *Numerical data and code for temporally correlated ballistic deposition*, version data-20261119, GitHub (2026), https://github.com/HuangYong2006/bd-correlated-noise-1.0.0-code/releases/tag/data-20261119。

机器可读信息见 [CITATION.cff](CITATION.cff)，现有许可见 [LICENSE](LICENSE)。当前提供固定版本链接，尚未获得 DOI。

## 结论边界

严格二元两时统计匹配下的结果是有限尺度响应。长程比较仍有协方差残差及收敛限制，不据此宣称普适类改变。

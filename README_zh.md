> 2026-10-31 更新：论文、PDF、七张图及完整创新变量对照数据已同步。
> 新增三个数据目录各15份NPZ，共45份；保留原失败验证、重新校准和追加独立验证记录。
> 最终15/15组满足原定0.005精度。最新复现入口见
> [本轮复现说明](docs/innovation_20261031/REPRODUCE.md)。

# 时间关联噪声下的弹道沉积

当前英文稿题目为 **Noise construction and finite-size inference in temporally correlated ballistic deposition**，采用 PRE 格式。中文阅读摘要见 manuscript/abstract_zh.md。

## 最新材料与独立核验

- [最新论文 PDF](manuscript/main.pdf)与[Overleaf 源码包](manuscript/overleaf_bd_paper.zip)。
- [数值补充 ZIP](manuscript/innovation_controls_data_20261031.zip)：包含45份NPZ、
  逐运行种子、重新校准记录、原始失败、分析与图7代码及独立核验脚本。
- [投稿材料总包](manuscript/PRE_submission_materials_20261031.zip)：PDF、源码、数据ZIP与README.TXT。
- [实际ZIP核验报告](docs/innovation_20261031/Supplemental_Material_核验报告.md)。

将补充ZIP解压到独立空目录后，运行
`python code/independent_supplement_audit.py` 可在不导入项目生成器或分析函数的情况下，
复算原失败和最终验证、有效斜率及bootstrap区间，并按显式种子再生校准记录。
本轮还通过了70项测试。原105组生产数据与新增15组参数单元是不同的实验网格。
结论仍限于已验证的协方差精度、有限尺寸与观测窗口。

## 2026-10-26 科学修订

旧生产谱滤波路径存在时间轴与空间轴排列错误，现已修复并增加端到端时空相关测试。原 spec、fgnspec 结果仅在 Git 历史中保留追溯，但不得作为“正确发生器”结果使用。之前的投稿准备状态和发生器偏差结论由本次修订取代。

三条新数据轨道 rev20261026_shared、rev20261026_powerlaw、rev20261026_fgn 均已重算。每条包括 7 个 θ、5 个尺寸（32–512）和每单元 64 次独立运行，共 105 个单元、6720 次运行。新谱数据记录 noise_layout=spectral-run-time-site-v2，新主分析会拒绝旧布局。

新结果显示，三种输入的有限尺寸晚期宽度斜率均随关联增强而上升，平坦对比例均下降。θ=0.25 时，共享 AR 为 0.717±0.008，fGn 谱滤波为 0.678±0.007。误差为独立 run bootstrap 标准差，不包含有限时间和噪声截断系统误差。本研究没有测出 θ_c，也没有证明其他作者的实现错误。旧稿关于恢复平坦 KPZ 分支、形貌变化方向相反、定量解释文献阈值的论断已撤下。

## 入口

- docs/REPRODUCE.md：当前复现命令。
- docs/DATA.md：数据字段及旧数据边界。
- code/bd_core.py：模拟器及已修正的时空场装配。
- code/validate_production_noise.py：最终生产接口的相关验证。
- code/revision_analysis.py：按独立 run 分析宽度、有效斜率及局域 pmf。
- code/export_revision_tables.py：从源表导出论文数值。
- code/revision_figures.py：四张生产图及源文件记录。
- code/submission_controls.py：配对分析、时长对照及第五张图。
- manuscript/main.tex、revision_results.tex、submission_results.tex：当前论文源文件。
- manuscript/overleaf_bd_paper.zip：重建后的 Overleaf 上传包。
- docs/revision_20261026/：修复、测试与视觉检查记录。

旧图表、旧分析结果、临时脚本及重复展开目录已从仓库当前版本移除，清理清单见 docs/CLEANUP.md。历史版本可通过 Git 历史找回。正式投稿仍需作者核对新结论和通信邮箱。

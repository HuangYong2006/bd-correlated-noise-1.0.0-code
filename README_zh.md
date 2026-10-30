> 2026-10-30 更新：当前论文为 `manuscript/main.tex`，PDF 为 `manuscript/main.pdf`。
> 本轮增加配对 bootstrap、四组固定输入记录时长对照以及第五张图。
> 三组修正版生产数据共105个单元已纳入仓库，其余历史原始数据仍排除。
> 最新复现步骤见 `docs/submission_20261030/REPRODUCE.md`。

# 时间关联噪声下的弹道沉积

当前英文稿题目为 **Noise construction and finite-size inference in temporally correlated ballistic deposition**，采用 PRE 格式。中文阅读摘要见 manuscript/abstract_zh.md。

## 2026-10-26 科学修订

旧生产谱滤波路径存在时间轴与空间轴排列错误，现已修复并增加端到端时空相关测试。原 spec、fgnspec 数据保留追溯，但不得作为“正确发生器”结果使用。之前的投稿准备状态和发生器偏差结论由本次修订取代。

三条新数据轨道 rev20261026_shared、rev20261026_powerlaw、rev20261026_fgn 均已重算。每条包括 7 个 θ、5 个尺寸（32–512）和每单元 64 次独立运行，共 105 个单元、6720 次运行。新谱数据记录 noise_layout=spectral-run-time-site-v2，新主分析会拒绝旧布局。

新结果显示，三种输入的有限尺寸晚期宽度斜率均随关联增强而上升，平坦对比例均下降。θ=0.25 时，共享 AR 为 0.717±0.008，fGn 谱滤波为 0.678±0.007。误差为独立 run bootstrap 标准差，不包含有限时间和噪声截断系统误差。本研究没有测出 θ_c，也没有证明其他作者的实现错误。旧稿关于恢复平坦 KPZ 分支、形貌变化方向相反、定量解释文献阈值的论断已撤下。

## 入口

- docs/REPRODUCE.md：当前复现命令。
- docs/DATA.md：数据字段及旧数据边界。
- code/bd_core.py：模拟器及已修正的时空场装配。
- code/validate_production_noise.py：最终生产接口的相关验证。
- code/revision_analysis.py：按独立 run 分析宽度、有效斜率及局域 pmf。
- code/export_revision_tables.py：从源表导出论文数值。
- code/revision_figures.py：四张修订图及源文件记录。
- manuscript/main.tex、revision_results.tex：当前论文源文件。
- manuscript/overleaf_bd_paper.zip：重建后的 Overleaf 上传包。
- docs/revision_20261026/：修复、测试与视觉检查记录。

旧上传包保存在 manuscript/archive/。旁边的 overleaf_bd_paper__1_.pdf 是旧稿，应重新编译当前源文件。正式投稿仍需作者核对新结论和通信邮箱。

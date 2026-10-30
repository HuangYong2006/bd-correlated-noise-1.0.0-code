# 修订稿最终证据一致性检查

检查日期 2026-10-26。此次为只读终检，不是三份互盲全审稿，也不继承先前模拟审稿结论。

## 结论

**Critical 0，Important 0。** 在本轮规定的证据范围内，未发现交付前必须修复的公式错误、数字矛盾、旧谱数据混入、图注与估计量错配，或将有限窗口结果升级为渐近物理结论的残留跳跃。两项 Minor 均为非阻断的表述澄清。

这一结论仅表示修订后的论证与指定证据相符，不表示已经测得渐近指数、相图或物理临界点，也不表示独立复现过原始长时模拟。

## 检查范围和方法

- 权威正文为 `manuscript/main.tex` 与 `manuscript/revision_results.tex`。四图为 `figs/revision_20261026/fig1_generator` 至 `fig4_morphology`，逐图查看 PNG 并核对绘图代码及源清单。
- 源表为 `data/analysis/revision_20261026/summary.json`、`noise_validation.json`。推导依据为 `docs/revision_20261026/theory_notes.md`。
- 只读核对了源表声明的全部 105 个 raw NPZ。逐一核验 SHA-256，读取宽度序列与最终剖面，重算各 cell 的宽度、平坦对比例、相邻时间窗比、Hellinger excess、条件 null 散布与 outer-run bootstrap 标准差。
- 依照已声明随机种子重算了 21 组五尺寸斜率及 84 组相邻尺寸斜率，包括各自的 500 次 bootstrap 误差。重算使用现有分析函数和原始数据，是数据到报告的一致性复算，不是另一套独立算法实现。
- 核对了 26 个标量 LaTeX 宏和表中 7 行数字，并从 noise validation 的逐 run 记录重算均值、SE 和容差判定。
- 文献核对只使用现有本地 `.lit` 中的 Squizzato–Canet、Alés–López、Rodríguez-Fernández 等、Jiang–Xia 原文。没有新检索，未将未读的 Wang–Xia 全文实现当作已经验证。
- 没有修改正文、图、代码或数据，没有阅读旧模拟审稿报告，没有新建子代理，没有重跑生产模拟。本次唯一写入为本报告。

## 已通过的关键检查

| 检查对象 | 核对结论 | 证据位置 |
|---|---|---|
| 网格与样本量 | 3 个输入 × 7 个 θ × 5 个 L = 105 cell，每 cell 64 runs，共 6720 runs。105 个 key 唯一，全部 raw 哈希匹配；各 raw 均含 32 个剖面。 | `main.tex` L304–319；`summary.json` L2–32、L792–2788、L2789 起的 sources |
| 时长 | AR 两端时长 2715、173779；谱输入 2880、174960，均与 raw metadata 一致。正文明确仅名义时长因子匹配，没有声称每个输入的精确 T 完全相同。 | `main.tex` L310–317；`summary.json` 中各 cell 的 T；`code/bd_core.py` L405–419 |
| 版本与图源 | 所有谱 raw 均为 `spectral-run-time-site-v2`。四图清单中的 summary/noise 哈希分别与当前源表完全匹配，未发现引用旧谱归档。 | `main.tex` L256–260；`code/revision_analysis.py` L35–37；`figure_sources.json` L2–5；全部谱 raw metadata |
| 宽度估计量 | Eq. (width) 与实际“run 内最后 20% 已存宽度平方取均值，再跨 run 平均后开方”一致。相邻时间窗比及 bootstrap 保留整条 run 的方法也一致。105 cell 的 W、W_err、p_flat、p_flat_err、late_ratio、late_ratio_err 全部复算匹配。 | `main.tex` L319–340；`code/revision_analysis.py` L28–32、L105–121 |
| 五尺寸与相邻尺寸斜率 | 21 组拟合和 84 组 pair 估计及误差全部复算匹配。θ=0.25 的 shared AR 为 0.7168428 ± 0.0082416，fGn-spectrum 为 0.6778455 ± 0.0066551，宏与表正确舍入为 0.717 ± 0.008、0.678 ± 0.007。 | `summary.json` L179–212、L683–716；`revision_results.tex` L2–4、L28–35；`main.tex` L399–443 |
| 两谱输入差异 | 全网格最大五尺寸斜率差为 0.0046787318，发生在 θ=0.30，宏正确舍入为 0.0047。正文明确两谱输入共享 Gaussian base streams，没有将相近结果当作独立重复。 | `revision_results.tex` L11；`main.tex` L315–317、L432–438；`summary.json` 的 21 组 exponents |
| 时间窗稳定性 | 最终/前一时间窗 RMS 宽度比最小 0.935019995、最大 1.053007437，正文和宏正确。正文同时保留“不是饱和证明”和初始化、时长、cutoff 的系统误差边界。 | `summary.json` L1679、L1964；`revision_results.tex` L22–23；`main.tex` L334–340、L444–447、L526–532 |
| 局部形貌 | 三个输入在每一个给定 L 上，观测到的 p_flat 均随所列 θ 严格下降。L=256 的四个引用值与宏一致；强相关时 shared/fGn 的 late-time W 比为约 1.30。没有旧的反向形貌论断残留。 | `main.tex` L449–459；`revision_results.tex` L15–21；`summary.json` L2180–2197、L2751–2767 及 shared 对应 cell |
| Hellinger 定义与不确定度 | 单步约束加周期性确实给出三点 pmf 和 Eq. (localH)。全部 ΔH 与误差复算匹配。200 次 outer independent-run bootstrap 重复 signal-minus-null 估计，保留同一次 resample 中共享的 reference；120 次条件分割的 null SD 是另存的量。θ=0 被明确设为参考定义下的零。 | `main.tex` L187–199、L342–378；`code/diagnostics.py` L163–221；`code/revision_analysis.py` L101–119；`summary.json` L27–30 |
| Null floor | fGn-spectrum 的 L=32、512 参考 null 均值为 0.0154622709、0.0035852800，比例 4.3127094；宏 0.0155、0.0036、4.31 正确。正文没有把首个显著性 crossing 认作物理 onset。 | `summary.json` L2133、L2209；`revision_results.tex` L12–14；`main.tex` L485–498 |
| 谱场验证 | 8 个验证 cell 的逐 run 均值、SE 和“5 SE + 0.003”容差均复算通过。L=32、θ=0.25 的 fGn 测量/target 为 0.2704129/0.2701902；同样样本的轴置换控制给出时间 0.0399911、空间 0.2632012，与宏一致。源表所载三个代码哈希仍匹配当前文件。 | `main.tex` L264–301、L382–390；`noise_validation.json` L772–958、L1525–1701、L1766–1777；`revision_results.tex` L24–27 |
| 四图图注 | Fig. 1 的八 run SE、finite-period target 与 raw AR covariance 区分正确。Fig. 2 包括全部四个尺寸 pair，横轴为几何中点。Fig. 3 使用 p_flat bootstrap 与 ΔH outer bootstrap。Fig. 4 的单条示例按最终宽度接近 run 中位数选取，谱为全部 64 条最终剖面均值 ± 1 SE，未做幅度归一、平滑、伪 floor 或谱指数拟合。 | `main.tex` L293–301、L424–429、L463–470、L475–482；`code/revision_figures.py` L75–117、L120–163、L182–218；`figure_sources.json` L38–58 |
| AR 公式与极限 | stationary covariance 保留 Q_nm 交叉项；shared Q_nm=1/12 和 independent diagonal 情形正确区分。cold-start covariance 和 shared transfer-function spectrum 一致。θ=0 的 log(τ)/τ 机制仅在 continuum-grid intermediate-range 近似下讨论，未误当有限 N 纯幂尾。Gaussian arcsine identity 未套用 uniform-driven AR。 | `main.tex` L201–233、L265–274、L392–397、L574–602；`theory_notes.md` L13–62、L64–126 |
| 局部量与大尺度结论 | Appendix 的 alternating/block 反例在相同局部 pmf 下确有不同宽度尺度。正文据此明确局部 ΔH 既不能独立排除、也不能建立 universality crossover。没有从 single-step 约束跳到“局部 pmf 足够识别粗糙度”。 | `main.tex` L198–199、L355–378、L518–534、L604–614；`theory_notes.md` L139 |

## 概念与文献阈值

1. `main.tex` L145–154 将 θ=1/6 的 continuum SR fixed-point stability、约 0.25 的 global/spectral 指数分裂与 faceting、height-distribution marker、emergence observable 分开。现有本地原文支持这种区分。
   - `.lit/squizzato_canet.txt` L1636–1643、L1677–1705 支持 SR/LR fixed-point 稳定性语境。
   - `.lit/ales_lopez2019.txt` L387–407 明示 L=2048、4096、8192，以及 θ_th=0.25±0.03 处 global/spectral roughness 分裂。
   - `.lit/ae8f20.txt` L399–425、L659–685 将约 0.2 的增长区高度分布 marker 与其他宏观/饱和区 marker 区分。主文没有把这些数字当作同一物理定义的矛盾估计。
2. `main.tex` L163–171、L217–221、L493–498、L536–541 撤回了对未验证历史代码和历史阈值来源的归因。对 Wang–Xia 的现有断言只限于 Hellinger/novelty 诊断已有先例，没有声称已重建其 detector、样本或判定规则。其全文实现本轮未核验。
3. `main.tex` L510–516、L543–554 已明确撤回“谱生成器恢复 θ≤0.25 的平坦 KPZ 分支”以及“上升主要由 shared AR 造成”的旧主张。当前数据确实支持所列窗口中三输入斜率同向上升，而非 asymptotic exponent 或新 θ_c。
4. `.lit/fragility.txt` L192–200、L504–522 已讨论 model/discretisation sensitivity。`main.tex` L153–154、L536–541 承认这一先例，未把所有微观敏感性、谱滤波或 Hellinger 方法据为首次发现。

## Critical / Important

本轮没有可核实的 Critical 或 Important 条目，也没有基于缺少新相图、更多尺寸或未取得全部历史代码而新增阻断要求。

## Minor

### FMR-m1　图内 fGn 简写可与全文命名统一

- 严重度 Minor。非阻断。
- 证据为 `main.tex` L248–254 明确使用 “fGn-spectrum filter”，但 Fig. 1(a) 标题、Figs. 2–4 图例仍仅写 “fGn”，对应 `code/revision_figures.py` L16–18、L89–90、L161。
- 图注和正文已经正确披露有限周期、DC removal 和 alias truncation，因此没有实质误导。可在下次导出时将简写统一为 “fGn spectrum” 或在图例中使用已定义缩写，避免图单独传播时被理解为精确无限 fGn 样本。

### FMR-m2　推导注记的离散 fGn 协方差可写明 lag 域

- 严重度 Minor。非阻断，涉及附属注记，未进入当前主文公式。
- 证据为 `theory_notes.md` L128–132。给出的差分式含 `(τ−1)^(2H)`，同时另列 C_fGn(0)=1，但没有在差分式旁明确 τ≥1。
- 可补上整数 τ≥1，或将三个底数全部取绝对值，写成 `0.5*(|τ+1|^(2H)−2|τ|^(2H)+|τ−1|^(2H))`。这样可避免在 τ=0 或负 lag 直接代入非整数次幂时的形式歧义。现有主文和数值目标没有据此出现错误。

## 范围边界

- 未重新估计 continuum θ、动态指数、渐近 roughness/spectral exponents、faceting threshold 或 universality class。当前手稿自己也没有声称完成这些任务。
- 未检验所有历史作者私有代码或完整 detector，实现归因保持在本仓库范围内。
- 未开展新文献搜索、完整参考文献逐字段审计或 Wang–Xia 全文复现。上述未覆盖事项不得被转述为“已验证”。
- 核验了已有容差的数值使用与报告一致性；没有独立审计 “predeclared” 一词所隐含的声明时间顺序。
- 检查了四图内容及其数据链，没有在本轮重新编译最终 PDF，也未审计中文摘要、投稿信、Overleaf 包或其他旁支文件。最终版式和交付包同步由集成流程负责。
- 对图源清单的追溯检查确认当前源表哈希与已列 raw 哈希一致；本轮没有逐字节独立重生成 PDF/SVG。因此“图一致”包括视觉检查和现有生成链核对，不声称进行了新的可重复渲染试验。

# BD 模型：长程时间关联噪声下的弹道沉积

> 本文件是 [`README.md`](README.md) 的中文对照说明，面向中文读者概述仓库内容、结论与使用方式。
> 完整的技术细节、字段定义与复现步骤以英文文档为准：
> [`docs/REPRODUCE.md`](docs/REPRODUCE.md)、[`docs/DATA.md`](docs/DATA.md)、
> [`docs/CODE_MAP.md`](docs/CODE_MAP.md)。

---

## 这个仓库是什么

(1+1) 维**弹道沉积（ballistic deposition, BD）**模型在**时间上长程关联的噪声**驱动下的
纯数值研究代码与数据。核心问题：文献中报告的临界时间关联指数 `theta_c`（KPZ 标度让位于
长程标度的地方）是否真的可以在现有系统尺寸下测到。

**结论是测不到，而原因是噪声发生器本身有问题。**

> **结论边界（请务必注意）**：本工作**没有给出 `theta_c` 的测量值**。
> 在 `L = 32 … 512` 范围内，使用关联结构正确的发生器时，`theta <= 0.25` 没有可分辨的
> crossover；重整化群预言 `theta_th = 1/6` 既未被证实也未被否定。
> 已确立的结论是：发生器偏差足以解释文献报告的阈值；并且标准的"固定显著性判据找
> onset"方法，在探测下限依赖系统尺寸时无法定位 crossover。

## 模型

关联函数（Squizzato & Canet 2019, Eq. 4）与离散模型（Alés & López 2019; Wang & Xia 2025;
Jiang & Xia 2026），并行更新、平坦初始条件、周期边界：

```
D(omega, k) = D0 + D_theta * omega^(-2 theta)   <=>   <eta(x,t) eta(x',t')> = 2 D delta_{x,x'} |t - t'|^(2 theta - 1)

h(x, t+1) = max[ h(x,t) + eta'_theta(x,t),  h(x-1,t),  h(x+1,t) ]
eta'_theta(x,t) = 1 if eta_theta(x,t) > 0 else 0
```

## 主要结论

1. **FFGN 发生器在 `theta = 0` 处不连续（核心发现）。**
   文献普遍使用 Mandelbrot 的快速分数高斯噪声（FFGN）。其归一化关联幅度基本与 `theta`
   无关（`A_FFGN ≈ 1.27`），而精确分数高斯噪声的幅度 `2 theta (theta + 1/2)` 在
   `theta = 0` 处消失。因此 `FFGN(theta -> 0+)` 并不是白噪声：实测
   `C(1)/C(0) = 0.195`，且 `C(tau) ∝ tau^-1` 尾巴。该效应**内禀**，不是截断假象
   （`N = 23 / 32 / 48` 个 AR 分量下拟合幅度为 `1.3558 / 1.3559 / 1.3559`）。
   由于模型要对噪声做二值化，标定必须在阈值化之后进行
   （`C_bin = (2/pi) arcsin C`）：标定后**整条名义 `theta` 轴都落在 RG 阈值 `1/6` 之上**，
   从未采样到短程 KPZ 区。

2. **模型是单步的（结构性结论）。** 0/1 噪声下从平坦初始条件出发 `|Delta h| <= 1` 恒成立，
   于是增量分布支撑集为 `{-1, 0, +1}`；而周期格子上 `sum_x Delta h = 0` 精确成立，故
   每个构型都有 `p_up = p_down`，`skew(Delta h) ≡ 0`（机器精度）。分布只有一个自由参数
   `p_flat`。因此本模型中的"faceting"不可能来自 overhang。

3. **标准 `theta_c` 定位方法是混淆的（自我更正）。** 在名义 FFGN 坐标下两个诊断给出
   `theta_c ≈ 0.0002`，但那正是发生器在 `theta = 0` 处的间断，不是物理测量。在修正坐标下，
   固定显著性 onset `theta*(L)` 随 `L` 单调下移，且强制过零的有限尺寸标度拟合形式上极好
   ——**但这是假象**：split-half 零假设下限随每轮列数（即 `L`）下降，所以判据本身就会漂移。

4. **决定性检验是固定 `theta`。** 真正的 crossover 要求偏差随 `L` **增大**；实测在
   `theta <= 0.25` 时 `DeltaH(theta, L)` 与 `L` 无关，而在 `theta >= 0.30` 时显著为正。
   即 crossover 是真实的，但在较小的 `theta` 处位于 `L ≈ 512` 之外。

5. **指数与诚实的误差棒。** 残差法误差把不确定度低估约 8 倍；改用按轮 bootstrap 后，
   修正发生器给出的 `alpha` 平坦（`theta = 0…0.25` 内约 0.49–0.54），而 FFGN 给出
   0.63–0.70。用**精确 fGn 谱密度**做的独立验证闭合了唯一可能的系统性来源，两个修正
   滤波器处处一致到 **0.008** 以内。负 `theta` 对照（`-0.2, -0.1, -0.05`）回到 `alpha = 1/2`。

6. **诊断方法：四个候选观测量中三个被数据否决。** 采用**精确离散 pmf**（不做直方图分箱）
   的 Hellinger 距离配合 split-half 零校正作为主诊断；直方图化的离散分布、以及任何保留
   涨落幅度的观测量都被证明会产生虚假的 `theta = 0+` 转变。第二个独立诊断
   （Isolation Forest，基于无量纲标度/形貌特征）在修正坐标下正确报告"无重构"。

## 仓库结构

```
code/          模型、诊断、分析、作图（扁平模块目录，import 依赖 code/ 在 sys.path 上）
data/raw/      每个 (theta, L) 一个 .npz；data/analysis/ 为派生分析表
figs/ctp/      论文四张图（pdf/svg/tiff + 对齐与碰撞审计 json）
manuscript/    论文源文件（REVTeX 4.2，目标期刊 Physical Review E）与 Overleaf 打包
docs/          复现步骤、数据结构、代码清单、历史遗留说明
tests/         核心不变量的单元测试
pack_repo.ps1  生成上传用压缩包
```

详细文件清单见 [`docs/CODE_MAP.md`](docs/CODE_MAP.md)。

## 快速开始

```powershell
pip install -r requirements.txt

# 先做验证（几分钟）
python code/validate_model.py

# 小型演练：1 个尺寸 x 3 个 theta x 64 轮
python code/run_campaign.py --tag demo --rule bd --L 64 --thetas 0 0.1667 0.40 --runs 64 --t-mult 15
python code/analyze.py --tag demo --rule bd --checkpoints 4

# 单元测试
python -m pytest tests -q
```

完整复现（含三个正式 campaign 的精确参数、运行时间与预期输出）见
[`docs/REPRODUCE.md`](docs/REPRODUCE.md)。

## 数据

`.npz` 文件按 `bd_p<theta>_L<L>.npz` 命名（`+`→`p`、`-`→`m`、`.`→`d`）。每个文件包含
`profiles`（各检查点的界面高度）、`w_series`（宽度时间序列，注意存的是 `W^2`）、
`stats`、`check_times`、`w_times`，以及**完整的 `SimSpec`**——后者是每个 cell 参数的权威记录。
字段逐项说明、campaign 清单、以及若干已知坑（含 `manifest` 只记录最后一次运行这一条）
见 [`docs/DATA.md`](docs/DATA.md)。

## 许可与引用

MIT 许可，见 [`LICENSE`](LICENSE)。引用信息见 [`CITATION.cff`](CITATION.cff)（其中
ORCID、仓库链接、发布日期等待填写项已用 TODO 标出）。

## 主要参考文献

Medina, Hwa, Kardar & Zhang, PRA **39**, 3053 (1989) · Lam, Sander & Wolf, PRA **46**, R6128
(1992) · Ma & Ma, PRE **47**, 3738 (1993) · Katzav & Schwartz, PRE **70**, 011601 (2004) ·
Fedorenko, PRB **77**, 094203 (2008) · Alés & López, PRE **99**, 062139 (2019) · Squizzato &
Canet, PRE **100**, 062143 (2019) · Song & Xia, PRE **103**, 012121 (2021) · Song & Xia,
J. Stat. Mech. (2021) 073203 · Hu, Hao & Xia, Physica A **619**, 128744 (2023) ·
Rodríguez-Fernández et al., PRE **110**, 024104 (2024) · Wang & Xia, PRE **111**, 024124
(2025) · Jiang & Xia, Commun. Theor. Phys. **78**, 105602 (2026)

# Task 1 核心谱噪声修复记录

日期为2026年10月26日。按用户已授权的 implementation_plan.md Task 1 在当前项目修改。没有创建worktree、提交git、重跑campaign或修改冻结评审及旧数据。

## 修改范围与决策

- 修改 code/bd_core.py、新增 tests/test_spectral_noise.py，并仅修改 tests/test_bd_core.py 中两处错误物理断言的docstring。
- 原生产代码把长度T的滤波行写入 `(run,T,L)` 的连续reshape，导致相关样本沿空间相邻。现在每行直接写入 `eta[run,:,site]`，仍返回C连续的 `(run,T,L)`、uint8数组。不引入另一个完整噪声场副本。
- 保留RNG初始化和样本生成顺序，测试要求更改行分块和run分块不改变同一种子的样本。旧版错误布局的空间场不会与修复后相同，这是本次修复的预期变化。
- simulate_spectral 返回 `noise_layout="spectral-run-time-site-v2"`。save_result保存此字段，load_result恢复为Python字符串。无此字段的旧归档明确读为 `"unversioned"`，不会自动升级为校正版。没有改变FFGN权重、共享创新、theta=0分支或更新器的数值路径。
- 清理文献实现完全相同、形状参数等于RG指数、有限DC移除记录严格为白噪声、未经限定的精度百分比等过强说明。谱函数数值公式未变。
- Ruling。依用户明确指令在当前目录工作，未执行技能中另建worktree、自动commit等通用建议。测试临时目录显式放在本任务目录内，因为系统默认pytest临时路径无法创建。

## RED

先写真实生产API测试，没有用一维校准器替代 generate_spectral_eta。

首次命令

```powershell
& 'D:/TOOLS/anaconda3/python.exe' -m pytest tests/test_spectral_noise.py -q
```

观察到4项时空相关性测试失败，6项分块和种子测试通过。另3项归档测试因系统临时目录创建失败而报setup错误，不能据此判断元数据功能。

修正测试运行环境后，仍未改生产代码，再运行

```powershell
& 'D:/TOOLS/anaconda3/python.exe' -m pytest tests/test_spectral_noise.py -q --tb=short --basetemp=docs/revision_20261026/.pytest_core_red
```

有效RED为7失败、6通过。4项真实时空相关性失败的时间一步相关如下，预期均大于0.15。

| spectrum | L | 旧生产输出时间一步相关 |
|---|---|---|
| powerlaw | 3 | 0.123563 |
| fgn | 3 | 0.130304 |
| powerlaw | 17 | 0.048242 |
| fgn | 17 | 0.050968 |

两项simulate_spectral归档测试及一项旧归档识别测试均因缺少noise_layout字段失败。所有失败来自被测功能，而非导入或测试语法问题。

## GREEN

修复布局、增加元数据后运行

```powershell
& 'D:/TOOLS/anaconda3/python.exe' -m pytest tests/test_spectral_noise.py -q --tb=short --basetemp=docs/revision_20261026/.pytest_core_green
```

13项全部通过。覆盖两种谱、L=1/3/7/17、单次和多次run、不整除的行/run分块、种子重复与种子变化，以及新归档往返和旧归档识别。时空统计测试使用theta=0.25、T=16384、3次实现、rows_per_chunk=5、每个run批次最多2次实现，含末尾不完整批次。

按完全相同统计参数直接读取修复后真实API，得到

| spectrum | L | 时间一步相关 | 空间最近邻相关 |
|---|---|---|---|
| powerlaw | 3 | 0.236634 | 0.000285 |
| powerlaw | 17 | 0.240144 | -0.000696 |
| fgn | 3 | 0.265553 | -0.001018 |
| fgn | 17 | 0.268509 | -0.001487 |

这些是固定种子的回归检查，不是生产活动的普适类或渐近指数验证。

## 基线与全套核查

在清理docstring后运行本任务加原基线

```powershell
& 'D:/TOOLS/anaconda3/python.exe' -m pytest tests/test_bd_core.py tests/test_spectral_noise.py -o addopts='' -q --tb=short --basetemp=docs/revision_20261026/.pytest_core_verified
```

结果为 `27 passed in 3.87s`。`git diff --check -- code/bd_core.py tests/test_bd_core.py` 通过。

并行开发期间已运行一次全套

```powershell
& 'D:/TOOLS/anaconda3/python.exe' -m pytest -q --tb=short --basetemp=docs/revision_20261026/.pytest_core_suite
```

当时其他代理刚写入理论任务的RED测试，得到8项 test_ffgn_covariance.py 失败，分别为单AR模式的2种shared参数、双AR模式的2种shared参数、脉冲响应有限和、归一化缩放不变性的2种shared参数，以及signed-noise-floor测试。其中前7项因待实现的ar_sum_covariance不存在失败，第8项是旧振幅拟合对负相关取绝对值的已知RED。其余测试通过。这不是全套最终GREEN，也没有越权修改理论任务文件；已将情况通知根代理，需所有并行任务完成后再次运行全套。

## 交接

已通知根代理核心接口可用于新tag重算。新旧数据应按noise_layout及新tag区分，不能把旧谱结果当作修复后结果。当前修复不声称解决有限记录、谱截断、交叉尺度或既有作者的实际算法问题。

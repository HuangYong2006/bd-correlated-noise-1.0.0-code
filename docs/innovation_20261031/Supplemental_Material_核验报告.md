# Supplemental Material 实际文件核验

核验对象：`innovation_controls_data_20261031.zip`。本报告针对实际ZIP解压内容，不以论文声明代替文件核验。

ZIP SHA256：`25f690993fd5914f7643b511ff3573f148d92da2259cd904c221fbd1333ca2df`。

## 包内实际内容

- 45份NPZ：15组六输入源数据、15组原始替代控制、15组最终控制。原失败文件未被删除。
- `DATA_INDEX.json`：列出初始失败、重新校准、追加验证的具体文件及数组字段。
- `SEEDS.json`：列出源过程、基线替代过程、重新校准、追加验证的逐运行种子或SeedSequence参数，以及bootstrap种子。
- `refinement_report.json`和`additional_validation_report.json`：两次方案的采样参数与前后上界。
- 数值分析、原始生成、重新校准、追加验证、图7重绘及两个核验脚本；依赖说明、许可、精确生产代码快照。
- 73个负载文件均有SHA256；另有MANIFEST.json，共74个ZIP成员。

## 直接复算的初始失败

| θ | L | 从原始运行记录复算的上界 | 原始归档 |
|---|---|---|---|
| 0.4 | 32 | 0.010739237 | `data/raw/innovation_binary_20261031_v3/theta0.40_L32.npz` |
| 0.4 | 64 | 0.005850207 | `data/raw/innovation_binary_20261031_v3/theta0.40_L64.npz` |
| 0.4 | 512 | 0.005309561 | `data/raw/innovation_binary_20261031_v3/theta0.40_L512.npz` |

## 重算范围与结果

- 独立脚本只使用NumPy与标准库，不导入项目生成器或分析函数。
- 重算全部15组最终协方差上界，均低于保留的0.005标准。
- 用独立最小二乘实现重算18组源输入斜率与bootstrap区间、3组替代对照差异及区间；重算全部四点矩差。
- 从拟合系数独立重建两组重新校准谱；每组按SEEDS.json再生第1条64位置校准记录，与归档协方差的最大差为0。这里的种子抽查各为1条，不声称本轮重新生成了全部Monte Carlo样本。
- 追加L512验证之外，原拟合谱、界面及原验证数组与保留的原文件逐元素一致。
- 使用随包项目分析代码重新运行全部分析后，summary的全部字段与归档一致。

## 你可以自行复核

在一个空目录解压数据ZIP，按README安装依赖后运行：

```text
python code/independent_supplement_audit.py
python code/verify_submission_data.py
```

第一条生成`INDEPENDENT_AUDIT.json`，显示初始失败、最终上界、斜率/区间与种子再生核验结果。完整时序可由逐运行种子及生成代码重建。完整重分析、重新生成校准和图7的命令见`README_REPRODUCE.txt`。

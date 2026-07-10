# dCS 引力波约束计算代码库

[English](README.md)

本仓库用于预测引力波观测对 dynamical Chern–Simons（dCS）引力理论的约束能力，包含探测器噪声模型、高阶模频域波形、dCS 振幅与相位修正、信噪比（SNR）、Fisher 矩阵、事件筛选、参数扫描和论文绘图流程。

当前分析主要关注的物理量为

```text
Delta sqrt(alpha_bar_dCS)
```

单位为 km。代码支持单探测器和多频段探测器组合，并包含观测时间 `T_tc` 与质量比 `q` 的专项敏感性分析。

## 主要功能

- 计算 TianQin、LISA、Taiji、DECIGO、BBO、ET、CE 和 aLIGO 的噪声功率谱密度。
- 通过 LALSuite 生成 `IMRPhenomXHM` 高阶模频域波形。
- 为 `(2,2)`、`(2,1)`、`(3,3)`、`(3,2)` 和 `(4,4)` 模加入 dCS 振幅与相位修正。
- 计算单探测器和双探测器的 SNR、Fisher 矩阵与协方差矩阵。
- 批量扫描质量、质量比、距离和事件参数。
- 根据可探测性与 dCS 有效性条件筛选事件。
- 分别研究 `T_tc` 和 `q` 对 dCS 约束的影响。
- 输出适用于论文正文或补充材料的 PDF 和 PNG 图像。

## 目录结构

```text
dCS/
├── README.md
├── README.zh-CN.md
├── Sn/
│   ├── Sn.py                  # 全天平均探测器噪声曲线
│   ├── Sn_NSA.py              # 非全天平均探测器噪声曲线
│   └── data/                  # ET、CE 和 aLIGO 灵敏度数据
├── math/
│   ├── fim.py                 # 波形、dCS 修正、SNR 和 Fisher 核心
│   ├── swsh.py                # 自旋权球谐函数
│   └── z_to_D.py              # 红移到距离的换算
├── getresults/
│   ├── onedetector.py         # 单探测器质量网格计算
│   ├── twodetectors.py        # 双探测器多频段计算
│   ├── filter.py              # 事件样本 dCS 约束
│   ├── snr.py                 # 事件样本 SNR
│   ├── readresult.py          # 结果读取和数组处理工具
│   ├── plot.py                # 质量—距离等高线绘图
│   ├── tc_sensitivity.py      # 观测时间敏感性计算
│   ├── q_sensitivity.py       # 质量比敏感性计算
│   └── plot_sensitivity.py    # 敏感性结果绘图
├── script/
│   ├── contour.py             # 质量网格扫描入口
│   ├── events.py              # 事件约束批处理入口
│   ├── events_snr.py          # 事件 SNR 批处理入口
│   ├── tc_sensitivity.py      # 观测时间扫描入口
│   ├── q_sensitivity.py       # 质量比扫描入口
│   └── plot_sensitivity.py    # 敏感性绘图入口
└── data/
    └── <YYYYMMDD>/
        ├── results/           # CSV 和文本结果
        └── fig/               # PDF 和 PNG 图像
```

## 环境依赖

- Python 3
- NumPy
- SciPy
- Matplotlib
- tqdm
- psutil
- LALSuite（`lal` 和 `lalsimulation`）

波形计算依赖 LALSuite，因此推荐使用 conda 环境：

```bash
conda create -n dcs python=3.11
conda activate dcs
conda install -c conda-forge lalsuite numpy scipy matplotlib tqdm psutil
```

如果服务器上的 Matplotlib 缓存目录不可写，可以设置临时缓存目录：

```bash
export MPLCONFIGDIR=/tmp/matplotlib
```

## 快速开始

以下命令均在仓库根目录运行。

### 绘制探测器灵敏度曲线

```bash
python Sn/Sn.py
```

### 单探测器与多频段质量扫描

```bash
python script/contour.py
```

运行前请检查脚本底部的日期标签、进程数、探测器列表和质量比。结果首先写入 `data/<date>/results/`，随后在 `data/<date>/fig/` 中生成图像。

### 计算事件样本的 dCS 约束

```bash
python script/events.py
```

该流程分批读取源参数，计算 dCS 约束，并输出满足有效性条件和 SNR 条件的事件。

### 只进行事件 SNR 筛选

```bash
python script/events_snr.py
```

该流程不计算 Fisher 矩阵，只计算 SNR 并筛选 `SNR >= 8` 的事件。

## 敏感性分析

仓库包含两个单变量扫描，用于检验 dCS 约束对固定分析参数选择的依赖。

### 观测时间

```bash
python script/tc_sensitivity.py
```

默认源参数：

```text
M_total = 100 Msun
D_l     = 100 Mpc
q       = 1.2
x1      = 0.40
x2      = 0.20
iota    = 0.9 rad
```

扫描的观测时间：

```text
1 天、1 周、1 个月、3 个月、0.5 年、
1 年、2 年、3 年、4 年、5 年
```

输出：

```text
data/<YYYYMMDD>/results/tc_sensitivity.csv
```

### 质量比

```bash
python script/q_sensitivity.py
```

默认源参数：

```text
M_total = 100 Msun
D_l     = 100 Mpc
T_tc    = 0.25 yr
x1      = 0.40
x2      = 0.20
iota    = 0.9 rad
```

扫描范围为 `q = 1.2, 1.3, ..., 8.9`，结果写入：

```text
data/<YYYYMMDD>/results/q_sensitivity.csv
```

扫描从 `q = 1.2` 开始，是因为 `q = 1` 对应物理边界 `eta = 0.25`，此时当前针对 `eta` 的前向有限差分会越过允许范围。

### 绘制敏感性分析结果

```bash
python script/plot_sensitivity.py
```

也可以指定已有结果的日期目录：

```bash
python script/plot_sensitivity.py --date 20260506
```

主要输出为：

```text
data/<YYYYMMDD>/fig/sensitivity_combined.pdf
data/<YYYYMMDD>/fig/sensitivity_combined.png
```

左图展示 dCS 约束和低频截止随观测时间的变化，右图展示 dCS 约束和 SNR 随质量比的变化。

## 输入与输出约定

事件输入文件每行包含 5 个由空格分隔的参数：

```text
Dl  m1  m2  Tc  theta_s
```

| 参数 | 含义 | 单位 |
|---|---|---|
| `Dl` | 光度距离 | Mpc |
| `m1`, `m2` | 两个分量的质量 | 太阳质量 |
| `Tc` | 距离并合的剩余观测时间 | s |
| `theta_s` | 当前传入 `iota` 的角参数 | rad |

敏感性 CSV 包含源参数、有效频率上下限、频率间隔、dCS 参数误差、有效性尺度和 SNR。论文图中使用的约束量为：

```python
d_sqrt_alpha_dCS_km = sigma_alpha_dCS_bar_square**0.25 * 1e-3
```

## Fisher 参数顺序

`math/fim.py::get_h_par` 使用以下参数顺序：

```text
0  M
1  eta
2  D_l
3  t_c
4  phi_c
5  x1
6  x2
7  alpha_dCS_bar_square
```

因此，`array_delta[7]` 对应 `alpha_dCS_bar_square` 的 1-sigma 误差。

## 探测器与多频段约定

- 探测器频段、扫描范围、噪声函数和调节后的有限差分步长统一注册在 `math/fim.py::get_detecror_parameters` 中。
- 地面探测器灵敏度通过 `Sn/data/` 中的数据文件插值得到。
- 在灵敏频段之外，探测器噪声曲线会被有意设置为非常大的功率谱密度。因此，在当前双探测器多频段权重中，其逆噪声贡献可以忽略。
- 当前的有限差分步长和频率间隔均为针对本任务设置的分析参数。

## 特定机器上的配置

部分脚本保留了生产计算所用机器的路径和运行参数：

- `Sn/Sn.py` 默认使用 `~/code/dCS/Sn/data`。
- `Sn/Sn_NSA.py` 中保留了 `/home/chexinyi/code/dCS/...` 路径。
- `script/events.py` 和 `script/events_snr.py` 使用 `/home/...` 与 `/data/...` 作为批量输入输出路径。
- 日期标签和多进程数量在各脚本入口中配置。

如果需要迁移到其他机器，请先检查这些设置。

## 主要调用关系

```text
script/contour.py
    -> getresults/onedetector.py
    -> getresults/twodetectors.py
    -> math/fim.py
    -> Sn/Sn.py
    -> getresults/plot.py

script/events.py
    -> getresults/filter.py
    -> math/fim.py
    -> getresults/readresult.py

script/events_snr.py
    -> getresults/snr.py
    -> math/fim.py
    -> getresults/readresult.py

script/tc_sensitivity.py
    -> getresults/tc_sensitivity.py
    -> math/fim.py

script/q_sensitivity.py
    -> getresults/q_sensitivity.py
    -> math/fim.py

script/plot_sensitivity.py
    -> getresults/plot_sensitivity.py
    -> data/<YYYYMMDD>/results/*.csv
```

## 扩展分析

- 新增探测器时，先在 `Sn/Sn.py` 中实现噪声曲线，再在 `math/fim.py::get_detecror_parameters` 中注册频段和数值参数。
- 修改 Fisher 参数顺序时，需要同步更新所有使用 `array_delta` 的模块，尤其是对索引 `7` 的引用。
- 修改输出格式时，需要同步更新 `getresults/readresult.py` 中对应的读取函数和敏感性绘图代码。

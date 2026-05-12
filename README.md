# dCS 引力波约束计算代码库

本仓库用于计算不同引力波探测器下 dynamical Chern-Simons gravity（dCS）修正参数的约束能力。代码包含探测器噪声谱、频域波形、dCS 修正、SNR、Fisher 矩阵、协方差矩阵、批量扫描、事件筛选和论文图绘制流程。

当前主要研究对象是不同探测器或探测器组合对

```text
Delta sqrt(alpha_bar_dCS)
```

的约束精度，以及这些约束对观测时间 `T_tc`、质量比 `q` 等固定参数选择的敏感性。

## 功能概览

- 计算 TianQin、LISA、Taiji、DECIGO、BBO、ET、CE、aLIGO 等探测器的噪声功率谱密度。
- 调用 `lalsimulation` 的 `IMRPhenomXHM` 高阶模频域波形。
- 在波形中加入 dCS 振幅修正和相位修正。
- 计算单探测器和双探测器联合观测下的 SNR、Fisher 矩阵和协方差矩阵。
- 批量扫描总质量、质量比、距离和事件参数。
- 针对审稿意见，计算固定源参数下 `T_tc` 和 `q` 的敏感性。
- 生成可用于论文正文或补充材料的 PDF/PNG 图像。

## 目录结构

```text
dCS/
├── README.md
├── Sn/
│   ├── Sn.py                  # 全天平均噪声谱；可绘制灵敏度曲线
│   ├── Sn_NSA.py              # 非全天平均噪声谱
│   └── data/                  # ET、CE、aLIGO 等插值用灵敏度数据
├── math/
│   ├── fim.py                 # 波形、dCS 修正、SNR、Fisher、协方差核心
│   ├── swsh.py                # 自旋权球谐函数
│   └── z_to_D.py              # 红移到宇宙学距离换算
├── getresults/
│   ├── onedetector.py         # 单探测器质量网格计算
│   ├── twodetectors.py        # 双探测器联合网格计算
│   ├── filter.py              # 事件样本 dCS 约束计算
│   ├── snr.py                 # 事件样本 SNR 计算
│   ├── readresult.py          # 结果文件读取、排序、重排工具
│   ├── plot.py                # 原有等高线绘图工具
│   ├── tc_sensitivity.py      # 观测时间 T_tc 敏感性计算
│   ├── q_sensitivity.py       # 质量比 q 敏感性计算
│   └── plot_sensitivity.py    # T_tc/q 敏感性结果绘图
├── script/
│   ├── contour.py             # 质量网格扫描入口
│   ├── events.py              # 事件样本 dCS 约束批处理入口
│   ├── events_snr.py          # 事件样本 SNR 筛选入口
│   ├── tc_sensitivity.py      # T_tc 敏感性扫描入口
│   ├── q_sensitivity.py       # q 敏感性扫描入口
│   └── plot_sensitivity.py    # 敏感性结果绘图入口
└── data/
    └── <date>/
        ├── results/           # CSV 或 txt 计算结果
        └── fig/               # PDF/PNG 图像输出
```

## 环境依赖

核心依赖：

- Python 3
- `numpy`
- `scipy`
- `matplotlib`
- `tqdm`
- `psutil`
- `lalsuite` / `lal` / `lalsimulation`

普通 Python 包可用 `pip` 安装：

```bash
python -m pip install numpy scipy matplotlib tqdm psutil
```

`lalsuite` 建议通过 conda 安装：

```bash
conda install -c conda-forge lalsuite
```

说明：

- 没有 `lalsuite` 时，`math/fim.py` 中的波形计算无法运行。
- 绘图脚本使用 `matplotlib` 的 `Agg` 后端，适合无图形界面的服务器。
- 如果服务器上 Matplotlib 缓存目录不可写，可设置：

```bash
export MPLCONFIGDIR=/tmp/matplotlib
```

## 快速运行

### 1. 绘制探测器灵敏度曲线

```bash
python Sn/Sn.py
```

输出灵敏度曲线文件到 `Sn/` 目录。

### 2. 扫描单探测器和双探测器质量网格

```bash
python script/contour.py
```

该入口会调用：

```text
getresults/onedetector.py
getresults/twodetectors.py
getresults/plot.py
```

结果写入：

```text
data/<date>/results/
data/<date>/fig/
```

运行前请检查脚本底部的 `date`、`ncpu`、探测器列表和质量比列表。

### 3. 对事件样本计算 dCS 约束

```bash
python script/events.py
```

该脚本读取服务器路径下的事件参数文件，计算事件样本的 dCS 约束，并把满足有效性条件的事件写入 `data/<date>/validity*/`。

### 4. 对事件样本只做 SNR 筛选

```bash
python script/events_snr.py
```

该脚本只计算事件 SNR，并筛选 `SNR >= 8` 的事件。

## 敏感性分析

为回应审稿人关于固定参数选择的意见，仓库新增了两个单变量敏感性计算。

### 1. 观测时间 T_tc 敏感性

```bash
python script/tc_sensitivity.py
```

默认固定参数：

```text
M_total = 100 Msun
D_l     = 100 Mpc
q       = 1.2
x1      = 0.40
x2      = 0.20
iota    = 0.9 rad
```

扫描：

```text
T_tc = 1 day, 1 week, 1 month, 3 months, 0.5 yr, 1 yr, 2 yr, 3 yr, 4 yr, 5 yr
```

输出：

```text
data/20260506/results/tc_sensitivity.csv
```

### 2. 质量比 q 敏感性

```bash
python script/q_sensitivity.py
```

默认固定参数：

```text
M_total = 100 Msun
D_l     = 100 Mpc
T_tc    = 0.25 yr
x1      = 0.40
x2      = 0.20
iota    = 0.9 rad
```

扫描：

```text
q = 1.2, 1.3, ..., 8.9
```

输出：

```text
data/20260506/results/q_sensitivity.csv
```

两个 CSV 都包含：

```text
detector
M_total_Msun
D_l_Mpc
q 或 T_tc_year
eta
x1, x2, iota
f_low_Hz
f_high_Hz
f_pn_low_Hz
f_isco_Hz
deltaF_Hz
sigma_alpha_dCS_bar_square
d_sqrt_alpha_dCS_km
alpha_dCS_bar_square_max_km
SNR
```

其中 `d_sqrt_alpha_dCS_km` 是论文图中直接展示的 dCS 约束精度，单位 km。

## 绘制论文图

敏感性分析结果跑完后，使用：

```bash
python script/plot_sensitivity.py
```

该脚本读取：

```text
data/20260506/results/tc_sensitivity.csv
data/20260506/results/q_sensitivity.csv
```

并输出：

```text
data/20260506/fig/sensitivity_constraints.pdf
data/20260506/fig/sensitivity_constraints.png
data/20260506/fig/sensitivity_diagnostics.pdf
data/20260506/fig/sensitivity_diagnostics.png
```

推荐论文正文使用：

```text
data/20260506/fig/sensitivity_constraints.pdf
```

该图包含两栏：

- `Delta sqrt(alpha_bar_dCS)` 随 `T_obs` 的变化；
- `Delta sqrt(alpha_bar_dCS)` 随 `q` 的变化。

诊断图：

```text
data/20260506/fig/sensitivity_diagnostics.pdf
```

用于展示：

- `f_low` 随 `T_obs` 的变化；
- SNR 随 `q` 的变化。

这张图适合放在补充材料或回复审稿人的说明中。

## 输入与输出格式

事件输入文件通常每行包含 5 列：

```text
Dl  m1  m2  Tc  theta_s
```

含义：

- `Dl`：光度距离，单位 Mpc。
- `m1, m2`：双星质量，单位太阳质量。
- `Tc`：并合前观测时间，单位秒。
- `theta_s`：当前代码中传入 `iota` 的角参数，单位弧度。

常见输出列：

- 输入事件参数；
- 8 个 Fisher 参数的不确定度；
- 总质量；
- 光度距离；
- dCS 参数有效性上界；
- SNR。

## 核心参数约定

`math/fim.py::get_h_par` 中 Fisher 参数顺序为：

```text
M
eta
D_l
t_c
phi_c
x1
x2
alpha_dCS_bar_square
```

因此结果中常用的：

```python
array_delta[7]
```

对应 `alpha_dCS_bar_square` 的 1-sigma 误差。绘图时通常转换为：

```python
d_sqrt_alpha_dCS_km = array_delta[7] ** 0.25 * 1e-3
```

## 注意事项

- `Sn/Sn.py` 中噪声谱数据目录默认为 `~/code/dCS/Sn/data`。如果仓库路径变化，需要同步修改。
- `Sn/Sn_NSA.py` 中仍保留服务器绝对路径 `/home/chexinyi/code/dCS/Sn/data/...`，本地运行前需要检查。
- `script/events.py` 和 `script/events_snr.py` 中有硬编码服务器路径，例如 `/home/chexy/...`、`/data/chexy/...`，换机器运行前需要修改。
- `script/q_sensitivity.py` 和 `script/tc_sensitivity.py` 当前最多使用 30 个进程。
- `q=1` 是 `eta=0.25` 的物理边界，前向差分会越界；当前 q 敏感性扫描从 `q=1.2` 开始，避免该问题。
- 当前仓库中 `.DS_Store` 是 macOS 自动生成文件，不参与计算。

## 代码调用关系

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
    -> data/20260506/results/*.csv
```

## 开发建议

- 新增探测器时，先在 `Sn/Sn.py` 中增加噪声谱函数，再在 `math/fim.py::get_detecror_parameters` 中注册频段、质量范围、距离范围和差分步长。
- 修改 Fisher 参数顺序时，需要同步检查 `math/fim.py::get_h_par`、`getresults/filter.py`、`getresults/onedetector.py`、`getresults/twodetectors.py`、`getresults/tc_sensitivity.py`、`getresults/q_sensitivity.py` 和绘图脚本中对第 8 个参数的使用。
- 如果结果文件格式变化，需要同步更新 `getresults/readresult.py` 和 `getresults/plot_sensitivity.py`。

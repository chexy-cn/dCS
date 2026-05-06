# dCS 引力波计算代码库

本仓库用于计算不同引力波探测器下的信噪比、Fisher 矩阵、协方差矩阵，以及 dynamical Chern-Simons（dCS）修正参数的约束结果。代码按“噪声谱、波形/Fisher 核心、结果处理、批量脚本”拆分，便于在服务器上批量运行并生成图像。

## 主要功能

- 计算 TianQin、LISA、Taiji、DECIGO、BBO、ET、CE、aLIGO 等探测器的噪声功率谱密度。
- 调用 `lalsimulation` 的 `IMRPhenomXHM` 高阶模频域波形。
- 在波形中加入 dCS 振幅与相位修正。
- 计算单探测器和双探测器联合观测的 SNR、Fisher 矩阵和协方差矩阵。
- 批量扫描质量、质量比、距离和事件参数，并输出文本结果。
- 根据结果绘制约束等高线、SNR 阈值线和有效性区域。

## 目录结构

```text
dCS/
├── README.md                  # 项目说明
├── Sn/
│   ├── Sn.py                  # 全天平均噪声谱；绘制灵敏度曲线
│   ├── Sn_NSA.py              # 非全天平均噪声谱
│   └── data/                  # ET、CE、aLIGO 等插值用灵敏度数据
├── math/
│   ├── fim.py                 # 波形、dCS 修正、SNR、Fisher、协方差核心
│   ├── swsh.py                # 自旋权球谐函数
│   └── z_to_D.py              # 红移到距离的宇宙学换算
├── getresults/
│   ├── onedetector.py         # 单探测器网格计算
│   ├── twodetectors.py        # 双探测器联合网格计算
│   ├── filter.py              # 事件样本的 dCS 约束计算
│   ├── snr.py                 # 事件样本的 SNR 计算
│   ├── readresult.py          # 结果文件读取、排序、重排工具
│   └── plot.py                # 结果图像绘制
└── script/
    ├── contour.py             # 扫描参数并绘制等高线的批量入口
    ├── events.py              # 事件样本 dCS 约束批处理入口
    └── events_snr.py          # 事件样本 SNR 筛选批处理入口
```

## 运行环境

代码依赖 Python 科学计算栈和 LALSuite：

- Python 3
- `numpy`
- `scipy`
- `matplotlib`
- `tqdm`
- `psutil`
- `lalsuite` / `lal` / `lalsimulation`

建议在服务器或本地虚拟环境中安装依赖：

```bash
python -m pip install numpy scipy matplotlib tqdm psutil
```

`lalsuite` 通常需要通过 conda 或系统环境安装，例如：

```bash
conda install -c conda-forge lalsuite
```

## 基本用法

### 1. 绘制探测器灵敏度曲线

```bash
python Sn/Sn.py
```

脚本会使用 `Sn/data/` 中的灵敏度数据，对 ET、CE、aLIGO 等地面探测器做插值，并在 `Sn/` 目录下输出灵敏度曲线图。

### 2. 扫描单探测器和联合探测器结果

```bash
python script/contour.py
```

该脚本会：

- 对多个探测器和质量比进行网格扫描；
- 写入 `data/<date>/results/`；
- 调用 `getresults.plot` 生成合并图像。

运行前请在脚本底部确认：

- `date` 是否是目标输出目录；
- `ncpu` 是否适合当前机器；
- 探测器列表和质量比列表是否符合本次任务。

### 3. 对事件样本计算 dCS 约束

```bash
python script/events.py
```

该脚本读取服务器路径下的事件参数文件，计算单探测器或 TianQin+地面探测器联合观测的约束，并把满足有效性条件的事件写入 `data/<date>/validity*/`。

### 4. 对事件样本只做 SNR 筛选

```bash
python script/events_snr.py
```

该脚本只计算事件 SNR，并筛选 `SNR >= 8` 的事件。

### 5. 固定源参数时的敏感性分析

审稿人若要求说明固定观测时间、质量比等参数对结果的影响，可以运行两个单变量敏感性脚本：

```bash
python script/tc_sensitivity.py
python script/q_sensitivity.py
```

`script/tc_sensitivity.py` 固定总质量 `100 Msun`、光度距离 `100 Mpc`、质量比 `q=1.2`、自旋和倾角，只改变 `T_tc = 0.25, 0.5, 1, 5 yr`，输出 `data/<date>/results/tc_sensitivity.csv`。

`script/q_sensitivity.py` 固定总质量 `100 Msun`、光度距离 `100 Mpc`、观测时间 `T_tc=0.25 yr`、自旋和倾角，只扫描 `q = 1 ... 10`，输出 `data/<date>/results/q_sensitivity.csv`。

两个输出文件都会记录 `f_low`、`f_high`、`deltaF`、SNR、`sigma_alpha_dCS_bar_square` 和 `d_sqrt_alpha_dCS_km`，可直接用于画敏感性曲线或写入回复审稿人的表格。

## 输入与输出约定

事件输入文件每行通常包含 5 列：

```text
Dl  m1  m2  Tc  theta_s
```

含义为：

- `Dl`：光度距离，单位 Mpc。
- `m1, m2`：双星质量，单位太阳质量。
- `Tc`：并合前观测时间，单位秒。
- `theta_s`：倾角或角参数，按当前公式传入 `iota`。

主要输出文件为逗号或制表符分隔文本。常见列包括：

- 输入事件参数；
- 8 个 Fisher 参数的不确定度；
- 总质量；
- 光度距离；
- dCS 参数有效性上界；
- SNR。

## 重要注意事项

- `math/fim.py` 依赖 `lalsimulation`，没有 LALSuite 时只能做语法检查，不能完整运行波形计算。
- 当前脚本中部分服务器路径是硬编码路径，例如 `/home/chexy/...`、`/data/chexy/...`。换机器运行前需要检查并修改这些路径。
- `Sn/Sn.py` 使用当前用户家目录下的 `~/code/dCS/Sn/data` 作为噪声谱数据目录；如果仓库位置变化，需要同步修改路径或改成相对路径。
- `Sn/Sn_NSA.py` 仍保留服务器绝对路径 `/home/chexinyi/code/dCS/Sn/data/...`，适合原服务器环境；本地运行前需要改为实际路径。
- 批量脚本会开启多进程，运行前应根据机器资源调整 `ncpu` 或 `SLURM_CPUS_PER_TASK`。
- 绘图脚本使用 `matplotlib` 的 `Agg` 后端，适合无图形界面的 Linux 服务器。

## 代码入口关系

```text
script/contour.py
    -> getresults/onedetector.py 或 getresults/twodetectors.py
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
```

## 开发建议

- 新增探测器时，优先在 `Sn/Sn.py` 增加噪声谱函数，再在 `math/fim.py::get_detecror_parameters` 中注册频段、质量范围、距离范围和差分步长。
- 修改 Fisher 参数顺序时，需要同步检查 `math/fim.py::get_h_par`、`getresults/filter.py`、`getresults/onedetector.py`、`getresults/twodetectors.py` 和绘图脚本中对 `array_delta[7]` 的使用。
- 若结果文件格式变化，需要同步更新 `getresults/readresult.py` 中对应的读取函数。

# dCS 论文绘图代码

[English](README.md)

本项目公开范围仅为论文 *Constraining the dynamical Chern-Simons gravity with
future gravitational-wave detectors* 中图 1–6 的计算、绘图代码及必要依赖。
不包含旧结果、预先计算的 CSV、图片、调试/验证脚本和族群表格计算流程。
catalog 生成见独立的 [population 项目](https://github.com/chexy-cn/population)。

## 安装

```bash
conda env create -f environment.yml
conda activate dcs
python -m pip install --no-deps --no-build-isolation -e .
```

## 重算论文图

批量计算需在计算节点运行：

```bash
python -m dcs noise                      # 图 1：噪声曲线
python -m dcs scan all --workers 8        # 图 2–6：重算数值
python -m dcs plot                       # 从新计算的 CSV 出图
```

默认输出为上海时区当日的 `results/<YYYYMMDD>/data/` 和 `fig/`，这些运行产物
不上传。也可用 `--date`、`--results-root` 指定目录。
[逐图复现说明](docs/reproducibility.md) 给出每张图的单独命令。

可用 `--quick` 试跑缩小网格，`--points` 调整采样数。已有扫描 CSV 默认禁止覆盖；
同配置可用 `--resume`，质量扫描需另选输出目录。改参数应选择新的日期/目录。

Slurm 入口保留，节点与分区由提交者决定：

```bash
sbatch --partition=YOUR_CPU_PARTITION --cpus-per-task=8 \
  script/run_sensitivity.slurm component
```

环境默认名为 `dcs`，可设置 `DCS_CONDA_ENV`；可用 `DCS_DATE` 指定批处理日期。

`Sn/data/` 中三份文件是噪声模型必需的输入表，并非历史计算结果。
源码许可证与噪声表上游授权仍需确认。
[科学约定说明](docs/scientific_conventions.md) 如实记录现有代码与论文的差异，
不将可运行流程等同于已验证所有论文数值。

## 贡献者

- [chexy-cn](https://github.com/chexy-cn)：项目维护与研究实现。
- [Codex（OpenAI）](https://openai.com/codex/)：辅助整理仓库、实现论文图流程、
  编写文档和执行验证。

Codex 为 AI 开发工具，科学解释与发布决策由人类作者及维护者负责。

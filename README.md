# dCS paper figures

[简体中文](README.zh-CN.md)

Calculation and plotting code for Figures 1–6 of *Constraining the dynamical
Chern-Simons gravity with future gravitational-wave detectors*.

This source release contains only the figure workflows and their dependencies.
It does not include historical results, generated figures, reference CSVs,
population catalogs, debugging/validation scripts or population-table calculations.
The [population project](https://github.com/chexy-cn/population) is separate.

## Install

```bash
conda env create -f environment.yml
conda activate dcs
python -m pip install --no-deps --no-build-isolation -e .
```

## Compute and plot

Run expensive calculations on an allocated compute node, not a management node:

```bash
python -m dcs noise                     # Figure 1
python -m dcs scan all --workers 8       # Figures 2–6
python -m dcs plot
```

Outputs default to `results/<today in Asia/Shanghai>/data/` (new numerical CSVs)
and `fig/` (PDF/PNG). All outputs are excluded from the source release.
Set `--date YYYYMMDD --results-root PATH` to choose a run directory.
The figure-specific commands are in [docs/reproducibility.md](docs/reproducibility.md).

For a small run: `python -m dcs scan q --detectors CE --quick --workers 2`.
Existing scan CSVs are protected; use `--resume` for an identical scan or a new
run directory. Mass scans require a new output CSV. Calculations record parameters,
dependency versions, source hashes and Git state for the newly generated results.

Slurm resources are site-specific and supplied when submitting:

```bash
sbatch --partition=YOUR_CPU_PARTITION --cpus-per-task=8 \
  script/run_sensitivity.slurm component
```

The environment defaults to `dcs`; override with `DCS_CONDA_ENV`. Optional
`DCS_DATE` selects the batch date. These are CPU jobs even on GPU-named nodes.

## Source layout

```text
dcs/           waveform/Fisher core, figure CLI and scans
getresults/    fixed-source calculations and paper plotting helpers
Sn/            detector noise functions and three required ASD input tables
script/        figure entry points and portable Slurm launcher
docs/          reproduction steps and implemented scientific conventions
```

The numerical environment is pinned. See
[scientific conventions](docs/scientific_conventions.md) for known differences
from the manuscript (fixed inclination, amplitude derivative and network
integration); successful execution alone does not establish exact scientific
agreement with every figure in the paper. Source licensing and upstream noise
table permissions still need author confirmation; no license is inferred here.

## Contributors

- [chexy-cn](https://github.com/chexy-cn): project maintainer and research implementation.
- [Codex (OpenAI)](https://openai.com/codex/): AI-assisted repository organization,
  figure-workflow implementation, documentation and execution checks.

Codex is an AI development tool. Scientific interpretation and release decisions
remain the responsibility of the human authors and maintainer.

# dCS Gravitational-Wave Constraints

[简体中文](README.zh-CN.md)

This repository forecasts constraints on dynamical Chern–Simons (dCS) gravity from gravitational-wave observations. It provides detector noise models, higher-mode frequency-domain waveforms, dCS amplitude and phase corrections, signal-to-noise ratio (SNR) calculations, Fisher-matrix forecasts, event filtering, parameter scans, and publication-ready figures.

The main quantity reported by the analysis is

```text
Delta sqrt(alpha_bar_dCS)
```

in kilometres. The current workflows support both individual detectors and multiband detector combinations, with dedicated sensitivity studies for the observation time `T_tc` and mass ratio `q`.

## Features

- Noise power spectral densities for TianQin, LISA, Taiji, DECIGO, BBO, ET, CE, and aLIGO.
- Higher-mode `IMRPhenomXHM` frequency-domain waveforms generated with LALSuite.
- dCS amplitude and phase corrections for the `(2,2)`, `(2,1)`, `(3,3)`, `(3,2)`, and `(4,4)` modes.
- Single-detector and two-detector SNR, Fisher-matrix, and covariance calculations.
- Mass, mass-ratio, distance, and event-parameter batch scans.
- Event selection using detectability and dCS validity criteria.
- One-variable sensitivity studies for `T_tc` and `q`.
- PDF and PNG figures suitable for papers and supplementary material.

## Repository layout

```text
dCS/
├── README.md
├── README.zh-CN.md
├── Sn/
│   ├── Sn.py                  # Sky-averaged detector noise curves
│   ├── Sn_NSA.py              # Non-sky-averaged detector noise curves
│   └── data/                  # ET, CE, and aLIGO sensitivity data
├── math/
│   ├── fim.py                 # Waveforms, dCS corrections, SNR, and Fisher core
│   ├── swsh.py                # Spin-weighted spherical harmonics
│   └── z_to_D.py              # Redshift-to-distance conversion
├── getresults/
│   ├── onedetector.py         # Single-detector mass-grid calculations
│   ├── twodetectors.py        # Two-detector multiband calculations
│   ├── filter.py              # dCS constraints for event samples
│   ├── snr.py                 # SNR calculations for event samples
│   ├── readresult.py          # Result readers and array utilities
│   ├── plot.py                # Mass-distance contour plots
│   ├── tc_sensitivity.py      # Observation-time sensitivity calculations
│   ├── q_sensitivity.py       # Mass-ratio sensitivity calculations
│   └── plot_sensitivity.py    # Sensitivity figures
├── script/
│   ├── contour.py             # Mass-grid workflow
│   ├── events.py              # Event-constraint workflow
│   ├── events_snr.py          # Event-SNR workflow
│   ├── tc_sensitivity.py      # Observation-time scan
│   ├── q_sensitivity.py       # Mass-ratio scan
│   └── plot_sensitivity.py    # Sensitivity plotting entry point
└── data/
    └── <YYYYMMDD>/
        ├── results/           # CSV and text results
        └── fig/               # PDF and PNG figures
```

## Requirements

- Python 3
- NumPy
- SciPy
- Matplotlib
- tqdm
- psutil
- LALSuite (`lal` and `lalsimulation`)

A conda environment is recommended because the waveform calculation requires LALSuite:

```bash
conda create -n dcs python=3.11
conda activate dcs
conda install -c conda-forge lalsuite numpy scipy matplotlib tqdm psutil
```

If the Matplotlib cache directory is not writable on a server, set a temporary cache location:

```bash
export MPLCONFIGDIR=/tmp/matplotlib
```

## Quick start

Run commands from the repository root.

### Detector sensitivity curves

```bash
python Sn/Sn.py
```

### Single-detector and multiband mass scans

```bash
python script/contour.py
```

Before running, review the date tag, process count, detector lists, and mass-ratio values near the bottom of the script. Results are written under `data/<date>/results/`, followed by figures under `data/<date>/fig/`.

### Event constraints

```bash
python script/events.py
```

This workflow reads batches of source parameters, evaluates dCS constraints, and writes events satisfying the configured validity and SNR criteria.

### Event SNR selection

```bash
python script/events_snr.py
```

This lighter workflow calculates SNR only and selects events with `SNR >= 8`.

## Sensitivity studies

The repository includes two one-variable scans designed to test the dependence of the dCS constraint on fixed analysis choices.

### Observation time

```bash
python script/tc_sensitivity.py
```

Default source parameters:

```text
M_total = 100 Msun
D_l     = 100 Mpc
q       = 1.2
x1      = 0.40
x2      = 0.20
iota    = 0.9 rad
```

Scanned observation times:

```text
1 day, 1 week, 1 month, 3 months, 0.5 yr,
1 yr, 2 yr, 3 yr, 4 yr, 5 yr
```

Output:

```text
data/<YYYYMMDD>/results/tc_sensitivity.csv
```

### Mass ratio

```bash
python script/q_sensitivity.py
```

Default source parameters:

```text
M_total = 100 Msun
D_l     = 100 Mpc
T_tc    = 0.25 yr
x1      = 0.40
x2      = 0.20
iota    = 0.9 rad
```

The scan covers `q = 1.2, 1.3, ..., 8.9` and writes:

```text
data/<YYYYMMDD>/results/q_sensitivity.csv
```

The scan begins at `q = 1.2` because `q = 1` corresponds to the physical boundary `eta = 0.25`, where the current forward finite difference in `eta` would cross the allowed range.

### Plotting the sensitivity results

```bash
python script/plot_sensitivity.py
```

To plot an existing result directory:

```bash
python script/plot_sensitivity.py --date 20260506
```

The main output is:

```text
data/<YYYYMMDD>/fig/sensitivity_combined.pdf
data/<YYYYMMDD>/fig/sensitivity_combined.png
```

The left panel shows the dCS constraint and low-frequency cutoff as functions of observation time. The right panel shows the dCS constraint and SNR as functions of mass ratio.

## Input and output conventions

Event input files contain five whitespace-separated columns:

```text
Dl  m1  m2  Tc  theta_s
```

| Column | Meaning | Unit |
|---|---|---|
| `Dl` | Luminosity distance | Mpc |
| `m1`, `m2` | Component masses | solar masses |
| `Tc` | Time remaining before coalescence | s |
| `theta_s` | Angular parameter passed to `iota` | rad |

Sensitivity CSV files include the source configuration, effective frequency bounds, finite frequency spacing, dCS uncertainty, validity scale, and SNR. The quantity plotted in the paper is

```python
d_sqrt_alpha_dCS_km = sigma_alpha_dCS_bar_square**0.25 * 1e-3
```

## Fisher parameters

`math/fim.py::get_h_par` uses the following parameter order:

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

Consequently, `array_delta[7]` is the one-sigma uncertainty of `alpha_dCS_bar_square`.

## Detector and multiband conventions

- Detector frequency limits, scan ranges, noise functions, and tuned finite-difference steps are registered in `math/fim.py::get_detecror_parameters`.
- Ground-based sensitivity data are interpolated from files in `Sn/data/`.
- Outside their sensitive bands, detector noise curves are deliberately assigned a very large power spectral density. Their inverse-noise contribution is therefore negligible in the two-detector weighting used by the current multiband workflow.
- The current finite-difference settings and frequency spacing are task-specific analysis choices.

## Machine-specific configuration

Several scripts retain paths and runtime settings for the machines on which the production calculations are run:

- `Sn/Sn.py` uses `~/code/dCS/Sn/data` by default.
- `Sn/Sn_NSA.py` contains `/home/chexinyi/code/dCS/...` paths.
- `script/events.py` and `script/events_snr.py` contain `/home/...` and `/data/...` paths for batch input and output.
- Multiprocessing limits and date tags are configured in the script entry points.

Review these values when moving the workflow to another machine.

## Main call graph

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

## Extending the analysis

- To add a detector, first implement its noise curve in `Sn/Sn.py`, then register its frequency range and numerical settings in `math/fim.py::get_detecror_parameters`.
- If the Fisher parameter order changes, update every consumer of `array_delta`, especially index `7` in the calculation and plotting modules.
- If an output format changes, update the corresponding readers in `getresults/readresult.py` and the sensitivity plotting code.

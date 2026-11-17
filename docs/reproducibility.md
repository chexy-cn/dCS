# Recompute paper Figures 1–6

No precomputed results are distributed. Generate numerical inputs from the
retained waveform/noise model, then plot those newly generated CSVs.

| Figure | Compute command | Output stem in dated `fig/` |
| --- | --- | --- |
| 1: detector ASD | `python -m dcs noise` | `detector_sensitivity` |
| 2–3: mass horizons, single/network | `python -m dcs scan mass --workers 8` | `all_detectors`, `TianQin_combined` |
| 4: observation time and mass ratio | `python -m dcs scan tc --workers 8` then `python -m dcs scan q --workers 8` | `sensitivity_combined` |
| 5: component spins | `python -m dcs scan component --workers 8` | `spin_component_sensitivity` |
| 6: inclination | `python -m dcs scan inclination --workers 8` | `inclination_sensitivity` |

After scans run `python -m dcs plot`. Output is PDF and PNG. `scan all` runs all
five scan families; `noise` is separate. Use the same `--date` and `--results-root`
for each command when overriding the default current date/directory. For another
input location, `plot --data-dir PATH` reads its CSVs without copying old results
into the repository.

## Default numerical choices

* Fixed source: M=100 Msun, Dl=100 Mpc, q=1.2, observation duration=0.25 yr,
  x1=0.4, x2=0.2, inclination=0.9 rad; each scan varies its named parameter.
* Time scan: 1 day, 1 week, 1 month, 3 months, 0.5, 1, 2, 3, 4 and 5 years.
* Mass ratio: 78 samples from 1.2 to 8.9.
* Component spins: independent 41-point grids on [-0.9,0.9]; the nonspinning
  constraint is NaN, with the coordinate retained.
* Inclination: 41 points on [0,pi/2].
* Mass: 1001 logarithmic samples over each detector's configured mass support,
  evaluated at a 1-Mpc reference distance for q=1.2 and 8.9.

Mass-horizon distances use exact scaling rho proportional to 1/D and the reported
length uncertainty proportional to D**0.25. This replaces the original raster
contour search; no old manual 600/700-Msun display mask is applied. Rendering
therefore need not be pixel-identical to historical figures. Figure 1 shows
sqrt(Sn), the amplitude spectral density.

Spin contours use powers of ten and an additional 50 km for TianQin/LISA/Taiji/
aLIGO. Stars mark (x1,x2)=(0.4,0.2). The inclination figure marks 0.9 rad.

## Environment and execution

`environment.yml` pins LALSuite 7.22, NumPy 1.26.4, SciPy 1.13.1 and
Matplotlib 3.9.2, with Python 3.11 as the public default. The historical tested
interpreter was Python 3.9.21. Exact BLAS/compiler builds are not locked;
finite-difference Fisher inversions can vary across platforms. Each run records
source hashes, library versions, parameters and Git state in a provenance JSON.

Run full calculations on compute nodes. `--quick` reduces scan sizes; `--points`
changes grid point count (except the fixed time list). For example:

```bash
python -m dcs scan q --detectors CE --quick --workers 2
```

Slurm launchers accept user-selected partitions/cores and require the configured
Conda environment. No particular hostname or personal directory is embedded.
Generated CSVs, plots and logs remain local under ignored `results/`.

Check mandatory detector-input integrity from the repository root:

```bash
sha256sum -c Sn/data/SHA256SUMS
```

See [scientific_conventions.md](scientific_conventions.md) before claiming numerical
agreement with the manuscript. Population realization generation and table
statistics are outside this figure-only release.

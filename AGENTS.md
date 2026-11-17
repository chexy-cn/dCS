# Repository Guidelines

## Publication Scope

Publish only code and required inputs for paper Figures 1–6: detector noise,
mass/network horizons, observation-time and mass-ratio dependence, component spins,
and inclination. Do not add historical CSVs, generated figures, logs, reference
results, debugging/validation scripts or population-table workflows to releases.
Population generation belongs in the separate population project.

## Structure

Scientific primitives are in `dcs/fim.py` and `dcs/swsh.py`. Detector noise
functions and the three required ASD inputs live in `Sn/`. Fixed-source helpers
and figure functions live in `getresults/`; `dcs/cli.py` is the public entry point.
Keep `script/` limited to paper-figure adapters and portable Slurm launchers.

New outputs belong in `results/<YYYYMMDD>/data/` and `fig/`, using today's
Asia/Shanghai date by default. These generated files remain ignored. Record
parameters, input/source hashes and dependency versions with each calculation.

## Environment and Running

Create the environment with `conda env create -f environment.yml`, activate `dcs`,
and install using `python -m pip install --no-deps --no-build-isolation -e .`.
Use `python -m dcs noise`, `python -m dcs scan all --workers 8`, then
`python -m dcs plot`. Individual commands are documented in
`docs/reproducibility.md`.

Run expensive calculations on allocated compute nodes, never the management
node. Use `script/run_sensitivity.slurm` with site-specific resources supplied to
sbatch. Do not embed personal paths, hostnames or scheduler partitions in code.

## Scientific and Coding Conventions

Use four-space indentation, snake_case and explicit units. Preserve the eight
Fisher parameter order with coupling at index 7. SNR uses GR through the detector
band; Fisher stops at ISCO. The nonspinning constraint is NaN. Read
`docs/scientific_conventions.md` before changing amplitude derivatives, mass
conventions, network integration or numerical inversion.

Verify changed figure workflows on a reduced grid on a compute node, inspect CSV
coverage and finite values, and replot the affected figure. Keep checks local;
do not ship task-specific test scripts or outputs. Commit messages should state
the figure or scientific/workflow change in a short imperative summary.

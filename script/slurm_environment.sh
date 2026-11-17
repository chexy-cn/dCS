#!/usr/bin/env bash
# Sourced by Slurm launchers; configure the partition/node at submission time.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the repository root}"
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export DCS_NCPU="${SLURM_CPUS_PER_TASK:-1}"
export MPLCONFIGDIR="${TMPDIR:-/tmp}/dcs-matplotlib-${SLURM_JOB_ID:-local}"
mkdir -p "$MPLCONFIGDIR"
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u
conda activate "${DCS_CONDA_ENV:-dcs}"
set -u
run_date=${DCS_DATE:-$(TZ=Asia/Shanghai date +%Y%m%d)}
log_dir="results/${run_date}/data/logs"
mkdir -p "$log_dir"
exec >"$log_dir/${SLURM_JOB_NAME:-dcs}-${SLURM_JOB_ID:-local}.out" 2>"$log_dir/${SLURM_JOB_NAME:-dcs}-${SLURM_JOB_ID:-local}.err"
echo "Job ${SLURM_JOB_ID:-local} node=${SLURMD_NODENAME:-unknown} workers=$DCS_NCPU"

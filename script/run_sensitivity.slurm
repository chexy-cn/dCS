#!/usr/bin/env bash
# sbatch --partition=<site-partition> --cpus-per-task=8 script/run_sensitivity.slurm component
#SBATCH --job-name=dcs-sensitivity
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm-%j.out
source "${SLURM_SUBMIT_DIR}/script/slurm_environment.sh"
scan=${1:-all}
shift || true
if [[ "$scan" == "plot" ]]; then
    srun --cpu-bind=none python -m dcs plot --date "$run_date" "$@"
else
    srun --cpu-bind=none python -m dcs scan "$scan" --date "$run_date" --workers "$DCS_NCPU" "$@"
    srun --cpu-bind=none python -m dcs plot --date "$run_date"
fi

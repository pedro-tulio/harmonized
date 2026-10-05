#!/bin/bash
#SBATCH --job-name=pyadmd-dist-pff
#SBATCH --partition=gpu_debug
#SBATCH --gpus=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:15:00
#SBATCH --output=slurm-pff-%j.out
#SBATCH --error=slurm-pff-%j.err

set -euo pipefail

module load py_env/miniconda3/25.7.0
module load atomic_simu/charmm/c50b2

unset PYTHONHOME
export PYTHONNOUSERSITE=1

DIST="$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)"
export PYTHONPATH="$DIST/src"

export PYADMD="$HOME/.conda/envs/pyadmd-openmm/bin/pyadmd"

echo "===== pyAdMD distribution validation ====="
echo "DATE       : $(date)"
echo "HOST       : $(hostname)"
echo "DIST       : $DIST"
echo "PYTHONPATH : $PYTHONPATH"

"$HOME/.conda/envs/pyadmd-openmm/bin/python" - <<'PY'
import pyadmd
import pyadmd.simulation.runner
import pyadmd.modes.subspace
print("pyadmd   =", pyadmd.__file__)
print("runner   =", pyadmd.simulation.runner.__file__)
print("subspace =", pyadmd.modes.subspace.__file__)
PY

echo
bash "$SLURM_SUBMIT_DIR/run-pff-charmm.sh"

echo
echo "===== TEST COMPLETED ====="
date

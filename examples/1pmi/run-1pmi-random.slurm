#!/bin/bash
#SBATCH --job-name=pyadmd-1pmi-random
#SBATCH --partition=gpu_debug
#SBATCH --gpus=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-1pmi-%j.out
#SBATCH --error=slurm-1pmi-%j.err

set -euo pipefail

# Directory containing this SLURM script
EXAMPLE="$(cd "$(dirname "$0")" && pwd)"
DIST="$(cd "$EXAMPLE/../.." && pwd)"
INPUTS="$EXAMPLE/inputs"

module load py_env/miniconda3/25.7.0
module load atomic_simu/charmm/c50b2

unset PYTHONHOME
export PYTHONNOUSERSITE=1

# Use the source tree from this distribution
export PYTHONPATH="$DIST/src"

PYADMD="$HOME/.conda/envs/pyadmd-openmm/bin/pyadmd"
PYTHON="$HOME/.conda/envs/pyadmd-openmm/bin/python"

echo "===== pyAdMD 1PMI random example ====="
echo "DATE       : $(date)"
echo "HOST       : $(hostname)"
echo "EXAMPLE    : $EXAMPLE"
echo "DIST       : $DIST"
echo "INPUTS     : $INPUTS"
echo "PYTHON     : $PYTHON"
echo "PYADMD     : $PYADMD"
echo "PYTHONPATH : $PYTHONPATH"
echo

echo "===== Verify imported pyAdMD ====="

"$PYTHON" - <<'PY'
import pyadmd
import pyadmd.simulation.runner
import pyadmd.modes.subspace

print("pyadmd :", pyadmd.__file__)
print("runner :", pyadmd.simulation.runner.__file__)
print("subspace:", pyadmd.modes.subspace.__file__)
PY

echo
echo "===== Starting 1PMI random simulation ====="

RUN="$EXAMPLE/run-random"

rm -rf "$RUN"
mkdir -p "$RUN"
cd "$RUN"

"$PYADMD" run \
  -src OPENMM \
  -psf "$INPUTS/step3_pbcsetup.psf" \
  -pdb "$INPUTS/step5_1.pdb" \
  -rst "$INPUTS/equil-final.xml" \
  -str "$INPUTS/step3_input.str" \
  -m CA \
  -sel "protein or name ZN" \
  -ek 1 \
  -t 250 \
  -nm 7,8,9 \
  --recalc-method random \
  -rep 1 \
  -seed 42 \
  --cycle-steps 50 \
  --nh-frequency 1.0 \
  -r

echo
echo "===== 1PMI random simulation completed ====="
date

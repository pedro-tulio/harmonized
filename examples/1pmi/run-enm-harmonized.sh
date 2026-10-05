#!/bin/bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
INPUTS="$HERE/inputs"
RUN="$HERE/run-enm-harmonized"

PYADMD="${PYADMD:-pyadmd}"

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
  -ek 6 \
  -t 5 \
  -nm 7,8,9 \
  -rep 1 \
  -seed 42 \
  --cycle-steps 50 \
  --nh-frequency 1.0 \
  --recalc-method harmonized \
  --subspace-search-modes 7,8,9,10,11,12 \
  -r

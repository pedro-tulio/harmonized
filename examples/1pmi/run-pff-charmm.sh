#!/bin/bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
INPUTS="$HERE/inputs"
RUN="$HERE/run-pff-charmm"

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
  -m CHARMM \
  -mod "$INPUTS/1pmi-3.mod" \
  -sel "protein or name ZN" \
  -ek 2 \
  -t 1 \
  -nm 7,8,9 \
  -rep 1 \
  -seed 42 \
  --cycle-steps 50 \
  --nh-frequency 1.0

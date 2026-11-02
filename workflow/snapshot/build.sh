#!/usr/bin/env bash
set -euo pipefail
cd /path/to/simc
make FFLAGSA='-O -w -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8 -fallow-argument-mismatch'

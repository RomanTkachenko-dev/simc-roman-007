#!/usr/bin/env bash
set -euo pipefail
cd /home/severus/eic/ACPOSTUMA_SIMC_GFORTRAN
make FFLAGSA='-O -w -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8 -fallow-argument-mismatch'

#!/usr/bin/env bash
set -euo pipefail
output=/path/to/simc-workspace/work/omega_after_rho_fix_20261005
input=/path/to/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/runs/omega_Roman_007/omega_Roman_007.inp
name=omega_Roman_007
for variant in before after; do
    if [ "$variant" = before ]; then cd /opt/simc/roman007-b6025c4; else cd /opt/simc/roman007-rho-jacobian-fix-20261005; fi
    mkdir -p "$output/$variant"
    if [ -e "infiles/$name.inp" ]; then cmp "$input" "infiles/$name.inp"; else cp "$input" "infiles/$name.inp"; fi
    test ! -e "worksim/$name.bin"
    printf '%s\n' "$name" | timeout 240 ./simc > "$output/$variant/simc.log" 2>&1
    test -s "worksim/$name.bin"
    cp "worksim/$name.bin" "$output/$variant/"
    cp outfiles/"$name"* "$output/$variant/"
    cp "infiles/$name.inp" "$output/$variant/"
    tail -n 4 "$output/$variant/simc.log"
done

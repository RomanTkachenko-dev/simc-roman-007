#!/usr/bin/env bash
set -euo pipefail
output=/mnt/c/Users/roman/.codex/.chatgpt-projects/g-p-6a6cd4ee3be0819189a11c2100a48a0c/work/rho_jacobian_fix_20261005
source_old=/opt/simc/roman007-b6025c4
source_fixed=/opt/simc/roman007-rho-jacobian-fix-20261005
input=/mnt/d/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/runs/rho_Roman_007/rho_Roman_007.inp
name=rho_Roman_007
mkdir -p "$output/before" "$output/after" "$output/source_changes"
cp "$input" "$output/$name.inp"
cd "$source_fixed"
make -C CTP FC=gfortran 'CFLAGS=-O -std=gnu89 -w -I/usr/include/tirpc -Df2cFortran -fcommon' > "$output/ctp_build.log" 2>&1
make 'FFLAGSA=-O -w -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8 -fallow-argument-mismatch' > "$output/simc_build.log" 2>&1
test -x simc
cp physics_rho_recoil.f NtupleInit.f "$output/source_changes/"
for variant in before after; do
    if [ "$variant" = before ]; then cd "$source_old"; else cd "$source_fixed"; fi
    if [ -e "infiles/$name.inp" ]; then cmp "$input" "infiles/$name.inp"; else cp "$input" "infiles/$name.inp"; fi
    test ! -e "worksim/$name.bin"
    mkdir -p worksim outfiles runout
    printf '%s\n' "$name" | timeout 240 ./simc > "$output/$variant/simc.log" 2>&1
    test -s "worksim/$name.bin"
    cp "worksim/$name.bin" "$output/$variant/"
    cp outfiles/"$name"* "$output/$variant/"
    cp "infiles/$name.inp" "$output/$variant/"
    tail -n 5 "$output/$variant/simc.log"
done

#!/usr/bin/env bash
set -euo pipefail
base=/mnt/c/Users/roman/.codex/.chatgpt-projects/g-p-6a6cd4ee3be0819189a11c2100a48a0c/work
old="$base/sign_check_2021"
output="$base/sign_comparison_2021_alicia"
cd /opt/simc/roman007-b6025c4
for file in test_eep_h_root.inp extra_hydrogen_root.inp test_pion_h_sign.inp extra_nucpi_sign.inp; do
    if [ -e "infiles/$file" ]; then
        cmp "$old/$file" "infiles/$file"
    else
        cp "$old/$file" "infiles/$file"
    fi
    cp "$old/$file" "$output/$file"
done
for name in test_eep_h_root test_pion_h_sign; do
    test ! -e "worksim/$name.bin"
    printf '%s\n' "$name" | timeout 240 ./simc > "$output/$name.log" 2>&1
    test -s "worksim/$name.bin"
    cp "worksim/$name.bin" "$output/"
    cp outfiles/"$name"* "$output/"
    tail -n 8 "$output/$name.log"
done

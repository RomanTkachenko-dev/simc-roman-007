#!/usr/bin/env bash
set -euo pipefail
output=/path/to/simc-workspace/work/omega_W247_Q235
cd /opt/simc/roman007-b6025c4
name=omega_W247_Q235
test ! -e "infiles/$name.inp"
cp "$output/$name.inp" "infiles/$name.inp"
printf '%s\n' "$name" | timeout 240 ./simc > "$output/simc.log" 2>&1
test -s "worksim/$name.bin"
cp "worksim/$name.bin" "$output/"
cp outfiles/"$name"* "$output/"
tail -n 12 "$output/simc.log"

#!/usr/bin/env bash
set -euo pipefail
repo=/home/severus/eic/ACPOSTUMA_SIMC_GFORTRAN
base=/mnt/d/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN
dest="$base/phi007_integration/smoke_$(date +%Y%m%d_%H%M%S)"
name=xphase_threebody_smoke
mkdir -p "$dest"
cd "$repo"
cp "infiles/$name.inp" "$dest/"
git rev-parse HEAD > "$dest/source_commit.txt"
timeout 180 ./simc <<< "$name" > "$dest/simc.log" 2>&1
test -s "worksim/$name.bin"
(cd util/root_tree && ./make_root_tree <<< "$name") > "$dest/conversion.log" 2>&1
cp "worksim/$name.root" "$dest/"
for file in outfiles/"$name"*; do
    [ ! -f "$file" ] || cp "$file" "$dest/"
done
echo "Smoke-test output: $dest"
tail -n 5 "$dest/simc.log"

#!/usr/bin/env bash
set -euo pipefail
repo=/home/severus/eic/ACPOSTUMA_SIMC_GFORTRAN
base=/mnt/d/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN
name=xphsp3_Roman_007
dest="$base/runs/$name"
if [ -e "$dest/$name.root" ]; then
    echo "Existing three-pion result preserved; inspect it before requesting a new run." >&2
    exit 1
fi
mkdir -p "$dest"
cd "$repo"
git diff --quiet HEAD -- .
git rev-parse HEAD > "$dest/source_commit.txt"
make FFLAGSA='-O -w -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8 -fallow-argument-mismatch' > "$dest/build.log" 2>&1
sha256sum simc > "$dest/executable_sha256.txt"
cp "$base/infiles/$name.inp" "infiles/$name.inp"
cp "infiles/$name.inp" "$dest/"
set +e
timeout 180 ./simc <<< "$name" > "$dest/simc.log" 2>&1
status=$?
set -e
echo "$status" > "$dest/exit_code.txt"
test "$status" = 0
grep -Eq 'Last Event.*1000  successes' "$dest/simc.log"
test -s "worksim/$name.bin"
(cd util/root_tree && ./make_root_tree <<< "$name") > "$dest/conversion.log" 2>&1
cp "worksim/$name.root" "$dest/"
for file in outfiles/"$name"*; do
    [ ! -f "$file" ] || cp "$file" "$dest/"
done
echo "Three-pion 10.60 GeV run saved: $dest"
tail -n 5 "$dest/simc.log"

#!/usr/bin/env bash
set -euo pipefail
base=/path/to/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN
repo=/path/to/simc
stamp=$(date +%Y%m%d_%H%M%S)
failed=0
python3 "$base/prepare_inputs.py"
cd "$repo"
mkdir -p worksim runout outfiles "$base/runs"
(
    cd util/root_tree
    gfortran -O3 -fno-automatic -c make_root_tree.f
    g++ -O3 $(root-config --cflags) -DMYROOT -c froot.c -o froot.co
    gfortran -o make_root_tree make_root_tree.o froot.co $(root-config --libs) -lstdc++
) > "$base/converter_build.log" 2>&1
for channel in omega phi rho eta etap xphsp; do
    name="${channel}_Roman_007"
    dest="$base/runs/$name"
    if [ -e "$dest/simc.log" ]; then
        # Exact fixed roots only; keep all previous evidence recoverable.
        archive="$base/archive/$stamp"
        mkdir -p "$archive"
        mv -- "$dest" "$archive/$name"
    fi
    mkdir -p "$dest"
    cp "infiles/$name.inp" "$dest/"
    echo "Starting $name"
    set +e
    timeout 180 ./simc <<< "$name" > "$dest/simc.log" 2>&1
    status=$?
    set -e
    echo "$status" > "$dest/exit_code.txt"
    if [ "$status" = 0 ] && grep -Eq 'Last Event.*1000  successes' "$dest/simc.log" && [ -s "worksim/$name.bin" ]; then
        (cd util/root_tree && ./make_root_tree <<< "$name") > "$dest/conversion.log" 2>&1
        cp "worksim/$name.root" "$dest/"
    else
        failed=1
    fi
    for file in outfiles/"$name"*; do
        [ ! -f "$file" ] || cp "$file" "$dest/"
    done
    echo "Finished $name: exit $status"
done
exit "$failed"

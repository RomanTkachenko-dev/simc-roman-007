#!/usr/bin/env bash
set -euo pipefail
output=/path/to/simc-workspace/work/sign_check_2021
dest=/opt/simc/bill-master-8b43f9d
if [ ! -d "$dest" ]; then
    mkdir "$dest"
    tar -xf "$output/source_8b43f9d.tar" -C "$dest"
fi
cd "$dest"
make -C CTP FC=gfortran 'CFLAGS=-O -std=gnu89 -w -I/usr/include/tirpc -Df2cFortran -fcommon' > "$output/ctp_build.log" 2>&1
python3 "$output/split_legacy.py"
flags='-O -w -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8 -fallow-argument-mismatch'
gfortran $flags -I. -c legacy_build/trg_deriv.f -o legacy_build/trg_deriv.o
gfortran $flags -I. -c legacy_build/trg_track_rest.f -o legacy_build/trg_track_rest.o
ld -r legacy_build/trg_deriv.o legacy_build/trg_track_rest.o -o trg_track.o
make 'FFLAGSA=-O -w -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8 -fallow-argument-mismatch' > "$output/simc_build.log" 2>&1
cp "$output/test_eep_h_root.inp" "$output/extra_hydrogen_root.inp" infiles/
mkdir -p worksim outfiles runout
printf '%s\n' test_eep_h_root | timeout 180 ./simc > "$output/test_eep_h_root.log" 2>&1
test -s worksim/test_eep_h_root.bin
cp worksim/test_eep_h_root.bin "$output/"
cp outfiles/test_eep_h_root* "$output/"
tail -n 8 "$output/test_eep_h_root.log"

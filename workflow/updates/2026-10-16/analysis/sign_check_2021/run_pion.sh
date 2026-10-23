#!/usr/bin/env bash
set -euo pipefail
output=/path/to/simc-workspace/work/sign_check_2021
python3 "$output/prepare_pion.py"
cd /opt/simc/bill-master-8b43f9d
printf '%s\n' test_pion_h_sign | timeout 180 ./simc > "$output/test_pion_h_sign.log" 2>&1
test -s worksim/test_pion_h_sign.bin
cp worksim/test_pion_h_sign.bin "$output/"
cp outfiles/test_pion_h_sign* "$output/"
tail -n 8 "$output/test_pion_h_sign.log"

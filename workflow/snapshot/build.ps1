$ErrorActionPreference = 'Stop'
# RPC headers exist in the WSL host but not in the current EIC container.
wsl.exe -d Ubuntu -- make -C /home/severus/eic/ACPOSTUMA_SIMC_GFORTRAN/CTP
if ($LASTEXITCODE -ne 0) { throw 'CTP build failed' }
wsl.exe -d Ubuntu -- /mnt/d/EIC/eic-shell -- /mnt/d/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/build.sh
if ($LASTEXITCODE -ne 0) { throw 'SIMC build failed' }

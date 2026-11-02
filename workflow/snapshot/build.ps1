$ErrorActionPreference = 'Stop'
# RPC headers exist in the WSL host but not in the current EIC container.
wsl.exe -d Ubuntu -- make -C /path/to/simc/CTP
if ($LASTEXITCODE -ne 0) { throw 'CTP build failed' }
wsl.exe -d Ubuntu -- /path/to/EIC/eic-shell -- /path/to/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/build.sh
if ($LASTEXITCODE -ne 0) { throw 'SIMC build failed' }

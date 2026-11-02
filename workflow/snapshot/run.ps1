$ErrorActionPreference = 'Stop'
wsl.exe -d Ubuntu -- /path/to/EIC/eic-shell -- /path/to/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/run_tests.sh
if ($LASTEXITCODE -ne 0) { throw 'SIMC tests failed: see run logs' }
wsl.exe -d Ubuntu -- /path/to/EIC/eic-shell -- /path/to/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/analyze.sh
if ($LASTEXITCODE -ne 0) { throw 'ROOT analysis failed' }

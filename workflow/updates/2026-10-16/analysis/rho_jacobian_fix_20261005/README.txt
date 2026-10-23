Recoil-rho corrections and paired reruns, 2026-10-05

Base version: Alicia/Acpostuma b6025c462ef33eeb620bdf4b399cd602c92e89bf.
Full patched source and executable:
\\wsl.localhost\Ubuntu\opt\simc\roman007-rho-jacobian-fix-20261005

Exactly two source changes were made:
1. physics_rho_recoil.f: peep_rho now uses Jttheta instead of Jttheta_fx.
2. NtupleInit.f: removed m=m+1 and the factor label from the eepx/Xphasespace
   branch only. The kaon branch factor entry remains intact.
The main%davejac and main%johnjac assignments were not changed.
All 179 Fortran/include files were compared with the base source archive.
Only the two requested files differ, with exactly the requested edits.
source_changes contains the modified files and a unified patch.

Build: gfortran 15, same modern-compiler options as the unmodified build.
Compilation logs: ctp_build.log, simc_build.log.
Input: rho_Roman_007.inp, byte-identical to the original saved input.
Beam 10.6 GeV; electron SHMS 6.7 GeV/c at 13 degrees;
proton HMS 1.1 GeV/c at 32 degrees. These are the original provisional settings.
No extra W=2.47 or Q2=2.35 retuning was made for this rho test.

before/: rerun with the original unmodified executable.
after/: rerun with the corrected executable.
Each contains native binary output, a converted ROOT h10 tree, the input,
run log, SIMC reports and saved initial random state.
Both runs: 138115 generation attempts, 1000 stored events, identical initial
random state and event kinematics. 573 events pass our usual analysis cuts.

Results among all 1000 stored events:
sigcm: 1000 positive before and after; exactly unchanged event by event.
siglab: 993 negative / 7 positive before; 0 negative / 1000 positive after.
Weight: 993 negative / 7 positive before; 0 negative / 1000 positive after.
After analysis cuts, all 573 corrected sigcm, siglab and Weight values are
positive. The same 573 events pass cuts before and after.

The laboratory-cross-section change agrees with the independently calculated
Jttheta/Jttheta_fx ratio, to a maximum absolute difference of 3.11e-14.
Weight changes by the same ratio. No absolute value or sign flip was applied.

The extra factor label previously shifted the eepx tail column names.
The corrected ntuple has 65 columns rather than 66. Every corrected tail
column was checked against its original underlying data slot.
Old ROOT files retain their old schema; the after ROOT uses the corrected one.

rho_before_after_fix.png and .svg show all 1000 stored events, no additional
analysis cuts or normalization. The histograms use signed native values and
one entry per event. verification_summary.json has complete statistics.

Scope: this demonstrates the requested fixes for the tested recoil-rho input.
No omega rerun or modification of the omega cross-section formula was made.

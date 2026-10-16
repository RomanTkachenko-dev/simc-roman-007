# Roman 007 research archive

Published: 2026-10-16. Main corrected-rho studies were performed on 2026-10-05; earlier studies retain their original filenames and records. No new physics simulation was performed for this publication.

## Current conclusions and limitations

1. **Recoil rho:** replacing `Jttheta_fx` with `Jttheta` in `peep_rho`, together with correcting the eepx/Xphasespace ntuple header, was tested in an isolated build. Both runs stored 1,000 events, of which the same 573 pass cuts. The corrected sample has positive `sigcm`, `siglab`, and `Weight`; `sigcm` itself is unchanged event by event. No absolute-value sign repair was applied.
2. **Omega:** the rho-only change does not change `physics_omega.f`. Paired baseline omega runs store 1,000 events, with 595 passing cuts. All 595 accepted `sigcm` values remain negative; positive laboratory weights do not validate this model. There is no omega Jacobian replacement or new omega parameterization in this archive.
3. **Combined signal:** the rho template now uses corrected event weights, and all seven nonnegative amplitudes were refitted with fixed shapes. The fitted omega count-equivalent contribution decreased by about 0.91%. The fit remains provisional, with substantial residuals and no finite-template-statistics uncertainty treatment. Its amplitudes are not absolute cross sections or physical rates.
4. **Historical 102 events/hour:** the normalization arithmetic is preserved, but the unresolved omega model means this is not a validated prediction. Do not use the historical presentation without this qualification.
5. **Henry comparison:** the plotted exponential curves are reference assumptions, not a validated replacement for the SIMC omega model or an extracted absolute cross section.

## Contents

| Directory | Contents |
|---|---|
| [analysis/rho_jacobian_fix_20261005](analysis/rho_jacobian_fix_20261005) | Exact source patch, modified files, inputs, initial random states, build/run logs, before/after native and ROOT output, verification and plots |
| [analysis/omega_after_rho_fix_20261005](analysis/omega_after_rho_fix_20261005) | Omega control runs before/after the rho-only changes; unchanged omega model |
| [analysis/combined_weighted_rho_20261005](analysis/combined_weighted_rho_20261005) | Updated seven-channel fit, histogram ROOT output, diagnostics, source scripts, Henry reference curves |
| [analysis/omega_cross_section_check](analysis/omega_cross_section_check) | Original signed omega cross-section checks and accepted-event information |
| [analysis/omega_W247_Q235](analysis/omega_W247_Q235) | Nominal W=2.47 GeV, Q2=2.35 GeV2 omega input and run outputs; negative-cross-section diagnostic |
| [analysis/sign_check_2021](analysis/sign_check_2021) | Original 2021 elastic/pion sign tests, not an omega validation |
| [analysis/sign_comparison_2021_alicia](analysis/sign_comparison_2021_alicia) | Identical-input 2021/Alicia comparisons |
| [historical_analysis](historical_analysis) | September overlays, weight/scaling studies, original fit, rate audit, three-pion run and original analysis scripts |
| [reports](reports/README.md) | Historical presentations and links to the bilingual report, with current caveats |

The six original channel runs are already retained under [workflow/snapshot/runs](../../snapshot/runs); they are not duplicated here. The additional three-pion run is under `historical_analysis/runs/xphsp3_Roman_007`.

## Exact source changes

Base: [acpostuma/simc_gfortran at b6025c462ef33eeb620bdf4b399cd602c92e89bf](https://github.com/acpostuma/simc_gfortran/tree/b6025c462ef33eeb620bdf4b399cd602c92e89bf).

- `physics_rho_recoil.f`: `peep_rho=Jttheta*gamma_T*Soding*ntup%sigcm` replaces the fixed-mass-Jacobian expression.
- `NtupleInit.f`: removes only the `m=m+1` and `NtupleTag(m)='factor'` pair in the eepx/Xphasespace branch. The kaon branch is unchanged.
- `physics_omega.f`: unchanged.

The full modified files and a unified patch are in [source_changes](analysis/rho_jacobian_fix_20261005/source_changes). Apply the patch to a separate case-sensitive checkout of the base source, not to the frozen September snapshot or unrelated working changes. The original tested build commands and compiler options are retained in [build_run.sh](analysis/rho_jacobian_fix_20261005/build_run.sh); its historical machine paths need adjustment on another computer.

The preserved patch has Windows line endings. From the separate base checkout, use `git apply --check --ignore-whitespace /absolute/path/to/requested_fixes.patch` first, then `git apply --ignore-whitespace /absolute/path/to/requested_fixes.patch`. The non-mutating check was verified against this repository's unchanged upstream source during publication. This does not change the archived patch bytes.

## Reproducibility and integrity

[COPIED_FILES_SHA256.json](COPIED_FILES_SHA256.json) records the size and SHA-256 of each byte-preserved copied research file. Publication documentation is not included in that copied-file manifest.

These are archived research scripts, not a newly validated portable launcher. They retain original Windows/WSL paths. Reading their ROOT/NumPy files and plotting scripts requires packages such as NumPy, uproot and Matplotlib; some older scripts use ROOT/PyROOT or SciPy. Local package caches are deliberately excluded. ROOT outputs are data/histograms; not every output contains a ROOT canvas or minimizer object.

The original raw Pythia/lAger event samples are external dependencies, not bundled here. Their derived comparison histograms and plots are preserved. Henry's externally supplied raw samples and standalone reference files are also not newly uploaded. No claim is made that every historical script can run on a fresh clone without resolving these dependencies and paths.

## Deliberate exclusions

- Credentials, private correspondence, meeting recordings, transcripts and screen captures of conversations.
- Compiler executables, installed Python packages, build caches, Git internals and duplicate upstream source archives.
- The canceled, incomplete 2014 omega build; it is not a completed physics result.
- Unrelated modifications in the local Bill `simc_gfortran` working tree.

All excluded originals remain untouched locally. The September tag and source snapshot are preserved.

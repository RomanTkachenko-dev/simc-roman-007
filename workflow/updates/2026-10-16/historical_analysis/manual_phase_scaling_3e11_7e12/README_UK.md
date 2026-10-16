# Ручне масштабування фазового простору — без фіту

Лінійна шкала. Відносно `meeting_20260916_weight_applied`:

- 2π: додатково ×3×10^11, зі збереженням попереднього ×4.
- 3π: додатково ×7×10^12.
- Пробне ×10^10 НЕ застосовується поверх цих множників.
- ω, φ, ρ, η, η′ та Pythia не змінені.

Weight застосовано до ω, 2π, 3π. Решта чотири SIMC-компоненти — counts,
як у попередньому графіку. Acceptance cuts збережені для всіх каналів.
Червона крива — арифметична змішана діагностична сума, не переріз.
Новий фіт, нормування до піка чи зміна форми не виконувалися.

`manual_phase_scales_linear.png/pdf` — графік.
`manual_phase_scales.root` — компоненти, сума, Pythia, canvas і metadata.
`summary.json` — множники, cuts, обмеження та перевірки.

Код у батьківській папці: `plot_manual_phase_scales.py`.
Запуск з PowerShell (наявні результати захищені від перезапису):

```
wsl.exe -d Ubuntu -- /mnt/d/EIC/eic-shell -- /mnt/d/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/plot_manual_phase_scales.sh
```

# Omega: Q², W та Bjorken x до і після застосування Weight

Використано існуючий `runs/omega_Roman_007/omega_Roman_007.root`, дерево h10, з початкового коду b6025c4. Пучок 10,60 ГеВ; SHMS: 6,70 ГеВ/c, 13°; HMS: 1,10 ГеВ/c, 32°. Це не новий трьохпіонний smoke-тест. Гілка Білла phi-007 не змінює джерело моделі omega.

Обидва варіанти графіків використовують ті самі 595 подій після cuts:

```text
abs(hsdelta)<8.0 && abs(hsxptar)<0.080 && abs(hsyptar)<0.035
&& ssdelta>-10.0 && ssdelta<22.0
&& abs(ssxptar)<0.040 && abs(ssyptar)<0.024
```

`No Weight`: кожна подія додає 1 до біна. `With Weight`: кожна подія додає своє оригінальне значення Weight. Додаткові коефіцієнти каналів і нормування до Pythia не застосовуються. У всіх цих omega-подій Weight додатний; це не вирішує окремі питання фізичної інтерпретації sigcm/ваг.

Q² та W прочитані з однойменних реконструйованих гілок. x_B обчислено для кожної події як `Q2/(W*W-mp*mp+Q2)`, де mp=0,938272 ГеВ для протона в спокої. Це Bjorken x, не координата xtar. Масу взято зі звіту SIMC.

## Результати

- `omega_Q2_W_x_weight_shapes_after_cuts.png/pdf`: три накладання, синій — без Weight, червоний — із Weight. Кожну гістограму окремо нормовано до одиничної площі, лише для порівняння форми.
- `omega_Q2_W_x_absolute_after_cuts.png/pdf`: шість панелей, верхній ряд — фактичні кількості подій, нижній — суми Weight. Числові підписи нижньої осі Y задано в одиницях 10⁻⁶; нормування до площі немає.
- `omega_Q2_weight_comparison`, `omega_W_weight_comparison`, `omega_xB_weight_comparison`: окремі PNG/PDF із накладаннями.
- `omega_kinematics.root`: початкові гістограми raw/weighted до cuts і після cuts, нормовані shape-гістограми після cuts, полотна та метадані.
- `summary.json`: повні визначення, binning, counts, середні, контрольні суми та перевірки.

| Середнє після cuts | Без Weight | З Weight |
|---|---:|---:|
| Q², ГеВ² | 3,751624 | 3,687783 |
| W, ГеВ | 2,054216 | 2,042599 |
| x_B | 0,528966 | 0,528268 |

Усі 595 подій потрапляють у видимі діапазони. Сума ваг дорівнює 5,1281449263×10⁻⁵. Перевірено counts, суму ваг, Sumw2, нульові underflow/overflow та одиничні площі shape-версій. Вхідний ROOT не змінено.

Запуск із PowerShell:

```powershell
wsl.exe -d Ubuntu -- /mnt/d/EIC/eic-shell -- /mnt/d/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/plot_omega_kinematics.sh
```

Код: `D:\EIC\SIMC_RUNS\ACPOSTUMA_SIMC_GFORTRAN\plot_omega_kinematics.py`.

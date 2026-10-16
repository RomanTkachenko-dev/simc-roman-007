# Розподіли самої гілки Weight для шести SIMC-файлів

Показано всі події, що проходять наведені нижче acceptance cuts, з дерев h10 у runs/<channel>_Roman_007/<channel>_Roman_007.root. Нових SIMC-запусків не виконувалося.

```text
abs(hsdelta)<8.0 && abs(hsxptar)<0.080 && abs(hsyptar)<0.035
&& ssdelta>-10.0 && ssdelta<22.0
&& abs(ssxptar)<0.040 && abs(ssyptar)<0.024
```

Вісь X — оригінальний Weight зі знаком; вісь Y — кількість подій у біні. Кожна відібрана подія додає ОДИНИЦЮ до гістограми. Це розподіл значень ваг, а не missing mass, заповнений вагами. Коефіцієнти каналів (/10, /22,5 тощо), abs(Weight) та нормування до Pythia НЕ застосовувалися.

Степінь десяти у заголовку осі X задає масштаб підписаних чисел. Наприклад, позначка −0,5 на осі Weight [×10⁻⁶] означає Weight = −5×10⁻⁷. Діапазони X різні для кожного каналу та охоплюють усі відібрані події й нуль. Пунктирна вертикаль показує Weight=0. У кожній гістограмі 80 бінів; underflow/overflow дорівнюють нулю.

| Канал | Подій після cuts | Weight < 0 |
|---|---:|---:|
| omega | 595 | 0 |
| phi | 549 | 536 |
| rho | 573 | 569 |
| eta | 584 | 584 |
| etap | 591 | 583 |
| xphsp | 523 | 0 |

`all_six_Weight_after_cuts.png/pdf` — аркуш із шістьма панелями для надсилання керівнику. `<channel>_Weight_after_cuts.png/pdf` — окремі графіки. `weight_distributions.root` містить гістограми з оригінальною числовою віссю Weight, полотна з масштабованими підписами осей та метадані. `summary.json` містить cut, min/max/mean, кількості ваг кожного знака та SHA256 джерел.

Перевірено кількість відібраних і від’ємних подій проти попереднього незалежного аудиту, інтеграли гістограм, похибки та незмінність джерел. Причину від’ємних ваг ці графіки самі по собі не встановлюють.

Запуск із PowerShell:

```powershell
wsl.exe -d Ubuntu -- /mnt/d/EIC/eic-shell -- /mnt/d/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/plot_weight_distributions.sh
```

Код: `D:\EIC\SIMC_RUNS\ACPOSTUMA_SIMC_GFORTRAN\plot_weight_distributions.py`.

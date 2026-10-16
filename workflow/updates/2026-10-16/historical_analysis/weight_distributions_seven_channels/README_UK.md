# Оригінальні Weight для семи SIMC-каналів

`all_seven_Weight_before_analysis_cuts.png/pdf` показує сім окремих розподілів: omega, phi, rho, eta, eta-prime, legacy X-phase-space та новий X-phase-space 3pi.

На цьому огляді використано **всі записані події h10 після SIMC, без додаткових analysis cuts**. Для старих шести файлів це по 1000 подій. Для нового 3pi smoke-тесту — 91 записана подія з 10 000 спроб генерації. Не слід плутати ці числа з кількістю подій після подальших acceptance cuts.

Вісь X — оригінальне значення Weight зі знаком. Степінь десяти в підписі осі задає масштаб чисел. Вісь Y — кількість подій у біні: кожна подія додає 1, а не Weight. Ніяких коефіцієнтів каналів, abs(Weight), нормування чи сумування каналів немає. Діапазони окремі для кожного каналу, охоплюють усі події та нуль; underflow/overflow дорівнюють нулю. Пунктир — Weight=0.

| Канал | Усі записані події | Weight < 0 |
|---|---:|---:|
| omega | 1000 | 0 |
| phi | 1000 | 963 |
| rho | 1000 | 993 |
| eta | 1000 | 1000 |
| eta-prime | 1000 | 973 |
| legacy X-phase-space | 1000 | 0 |
| X-phase-space 3pi | 91 | 0 |

## Різні кінематичні налаштування

Попередні шість файлів: Ebeam=10,60 ГеВ, SHMS/HMS, джерело b6025c4. Новий трьохпіонний smoke-тест: Ebeam=4,70897 ГеВ, SOS/HMS, джерело d778c56, radiation/energy-loss/smearing увімкнено. Це окрема мала тестова вибірка; її форма та масштаб ваг не є прямим фізичним порівнянням із шістьма іншими каналами. Відсутність від’ємних ваг у 91 події не встановлює їх відсутності для всіх можливих налаштувань.

Порівняння шести початкових каналів **після узгоджених acceptance cuts** залишилося в `../weight_distributions_after_cuts/all_six_Weight_after_cuts.png/pdf`. Ці SHMS/HMS cuts не застосовувалися автоматично до SOS/HMS smoke-файлу.

## Файли

- `all_seven_Weight_before_analysis_cuts.png/pdf`: загальний огляд.
- `<channel>_Weight_before_analysis_cuts.png/pdf`: окремі графіки.
- `seven_weight_distributions.root`: гістограми, полотна, метадані. Числова вісь збережених вихідних гістограм містить оригінальні Weight, а масштаб у полотнах лише полегшує читання підписів.
- `summary.json`: статистика знаків, діапазони, min/max/mean, SHA256 джерел.

Перевірено counts і знаки проти попередніх незалежних аудитів, кількості подій у гістограмах, Sumw2 та незмінність семи ROOT-джерел. Нових симуляцій не запускали.

Відтворення в PowerShell:

```powershell
wsl.exe -d Ubuntu -- /mnt/d/EIC/eic-shell -- /mnt/d/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/plot_seven_weight_distributions.sh
```

Код: `D:\EIC\SIMC_RUNS\ACPOSTUMA_SIMC_GFORTRAN\plot_seven_weight_distributions.py`.

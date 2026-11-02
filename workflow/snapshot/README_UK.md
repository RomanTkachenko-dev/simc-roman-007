# SIMC: виконання завдань зі зустрічі

Дата перевірки: 7 вересня 2026. **Технічні тести виконані; фізична валідація ще не завершена.**

## Що зроблено

1. Завантажено повний репозиторій https://github.com/acpostuma/simc_gfortran, commit `b6025c462ef33eeb620bdf4b399cd602c92e89bf`.
2. Створено окрему третю копію. Код Білла та фізичний код Генрі не змінювалися.
3. Зібрано SIMC. Допоміжний CTP зібрано в Ubuntu/WSL, де є RPC headers; Fortran SIMC зібрано і запущено в eic-shell. Для сучасного gfortran додано параметр сумісності `-fallow-argument-mismatch`.
4. Створено шість input-файлів. Реакційний блок кожного збережено з upstream-шаблону, за винятком `ngen=1000`. Усі блоки після experiment перенесено з локального Henry input. Виправлено лише описку в закривальному імені `beamandtergetinfo`.
5. Виконано всі шість запусків. Результати перетворено офіційними `make_root_tree.f` і `froot.c` з цього репозиторію в ROOT, дерево `h10`. Фізичний Fortran-код не змінювався.
6. У ROOT побудовано розподіли до/після заданих cuts: Q², W, missing mass, −t, xptar/yptar/ytar та delta обох плечей, Weight. Збережено PNG, PDF і ROOT-гістограми. Для кожного каналу є raw та weighted версії.
7. Окремо проаналізовано два дерева Pythia. Дублювання lAger перевірено заново за всіма 33 гілками і 112441 подіями; дублікат виключено.
8. Наданий native `.simc` перевірено і повторно запущено через версію Генрі. Після SIMC у ROOT реально 97 записів, після cuts — 50. Лічильник завершення в старому Fortran-лозі не використано як кількість ROOT-подій.

## Де лежать код і результати

Основна папка: `/path/to/EIC\SIMC_RUNS\ACPOSTUMA_SIMC_GFORTRAN`.

| Шлях у цій папці | Призначення |
|---|---|
| `prepare_inputs.py` | Відтворює input-файли з реакційних шаблонів і локального Henry baseline |
| `infiles\omega_Roman_007.inp` | Готовий omega input; поруч інші п'ять каналів |
| `input_audit.json` | Commit, джерела, SHA-256 і перевірка однакових блоків установки |
| `build.ps1`, `build.sh` | Збирання |
| `run.ps1`, `run_tests.sh` | Запуск шести тестів і ROOT-аналіз |
| `analyze_root.py`, `analyze.sh` | Код побудови графіків ROOT і застосування cuts |
| `runs\omega_Roman_007\omega_Roman_007.root` | Події omega після SIMC, дерево h10 |
| `runs\omega_Roman_007\plots\missmass_raw.png` | Missing mass omega до/після cuts, без ваг |
| `runs\omega_Roman_007\plots\missmass_weighted.png` | Той самий розподіл з Weight |
| `runs\<канал>_Roman_007\plots\histograms.root` | Гістограми цього каналу, а не нове дерево подій |
| `analysis_summary.json`, `weight_audit.json` | Кількості, ваги, діапазони й діагностика |
| `source_b6025c4.tar.gz` | Повний архів відстежуваного upstream-коду; нові input-файли лежать окремо |

Робочий Fortran-репозиторій зберігається в Linux, щоб не пошкодити файли, назви яких відрізняються регістром:

`/path/to/simc`

Відкрити його у Провіднику Windows можна, вставивши адресу:

`\\wsl.localhost\Ubuntu\home\severus\eic\ACPOSTUMA_SIMC_GFORTRAN`

Вхідні бінарні записи SIMC також залишаються там у `worksim/`. Копії готових ROOT і всі графіки — на диску D.

Окремі результати отриманих від керівника ROOT-файлів:

`/path/to/EIC\SIMC_RUNS\BACKGROUND_ROOT_ANALYSIS`

Окремий новий запуск Генрі:

`/path/to/EIC\SIMC_RUNS\HENRY_SIMC_FILE_INPUT\runs\meeting_check_20260907_104705`

Оригінальні файли в Downloads не змінювалися. Наданий `.simc` і копія в Henry inputs мають однаковий SHA-256; повний запис є в `BACKGROUND_ROOT_ANALYSIS\analysis_summary.json`.

## Фактичні результати

Позитивний `ngen=1000` означає 1000 успішних SIMC-подій, а не 1000 спроб генерації.

| Канал | Спроб у SIMC-лозі | Записів h10 | Після cuts | Від'ємних Weight до cuts |
|---|---:|---:|---:|---:|
| omega | 136898 | 1000 | 595 | 0 |
| phi | 144429 | 1000 | 549 | 963 |
| rho | 138115 | 1000 | 573 | 993 |
| eta | 163403 | 1000 | 584 | 1000 |
| eta-prime (`etap`) | 133143 | 1000 | 591 | 973 |
| X-phase-space | 24596 | 1000 | 523 | 0 |

Це незважені кількості тестових подій, не очікувані експериментальні yields і не абсолютна ефективність детекторів.

| Отриманий sample | До cuts | Після cuts | Як використано |
|---|---:|---:|---|
| Pythia_557ab.root:Background | 10974403 | 8914150 | Фон, raw та unit-area форми |
| Pythia_557ab.root:Signal | 112441 | 92401 | Сигнал, raw та unit-area форми |
| lAger_2zb.root:Signal | 112441 | — | Точний дублікат попереднього, не додано вдруге |
| Native phi `.simc` через Henry | 371782 вхідних подій → 97 записів ROOT | 50 | Окремий file-input маршрут |

У наданих Pythia/lAger деревах немає Weight. Числа в назвах файлів не використано як непідтверджені коефіцієнти нормалізації. Signal/background-форми нормовані кожна окремо до одиниці в показаному діапазоні; це не фізичне співвідношення сигнал/фон.

## Cuts

```cpp
abs(hsdelta)<8.0 && abs(hsxptar)<0.080 && abs(hsyptar)<0.035
&& ssdelta>-10.0 && ssdelta<22.0
&& abs(ssxptar)<0.040 && abs(ssyptar)<0.024
```

У ROOT для нового acpostuma output:

```cpp
auto f = TFile::Open("omega_Roman_007.root");
auto h10 = f->Get<TTree>("h10");
h10->Draw("Q2>>hQ2(100,0,10)",
 "Weight*(abs(hsdelta)<8.0 && abs(hsxptar)<0.080 && abs(hsyptar)<0.035"
 " && ssdelta>-10.0 && ssdelta<22.0 && abs(ssxptar)<0.040 && abs(ssyptar)<0.024)");
```

**Для цього acpostuma output гілка `t` вже містить −t.** Це прямо видно в `event.f`, рядки 1476–1477, і `results_write.f`, рядок 137. Для фізичного графіка −t тут потрібно малювати `t`, а не `-t`. Значення в сторонніх ROOT-файлах залишено під назвою `t branch`, бо їхній writer не надано.

`missmass` у acpostuma — відновлена signed missing mass із recon Em і Pm (`event.f`, рядки 1459–1461). Негативні значення не замінено модулем; позадіапазонні події враховані в ROOT underflow та JSON. Позначення `mx` у сторонніх ROOT-файлах не перейменовано на підтверджену missing mass без визначення автора.

Гілок `hsxtar`/`ssxtar` немає: xtar не підмінювався ytar чи координатою вершини. У поточному Henry h666 також немає `missmass` та `t`; ці два графіки для нього не вигадувалися. Missing mass із наданих 11-колонкових подій у validation report є **до SIMC**, а не після нього.

## Як запустити з Windows

Відкрити PowerShell. Повторити всі шість тестів і аналіз:

```powershell
& '/path/to/EIC\SIMC_RUNS\ACPOSTUMA_SIMC_GFORTRAN\run.ps1'
```

Попередні папки каналів перед повторним запуском переміщуються в `archive\<дата_час>`, не видаляються. Графіки аналізу ROOT перегенеровуються. `run_tests.sh` відхиляє незавершені SIMC-тести, навіть якщо Fortran повернув код 0.

Лише перебудувати графіки, не повторюючи симуляцію:

```powershell
wsl.exe -d Ubuntu -- /path/to/EIC/eic-shell -- /path/to/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/analyze.sh
```

Лише повторити native `.simc` у версії Генрі, у новій timestamp-папці:

```powershell
wsl.exe -d Ubuntu -- /path/to/EIC/eic-shell -- /path/to/EIC/SIMC_RUNS/HENRY_SIMC_FILE_INPUT/run_meeting_check.sh
```

Повторне збирання, якщо змінено Fortran-код:

```powershell
& '/path/to/EIC\SIMC_RUNS\ACPOSTUMA_SIMC_GFORTRAN\build.ps1'
```

`prepare_inputs.py` щоразу відтворює input-файли. Для іншого production input треба змінити джерело `HENRY` у цьому скрипті та знову перевірити `input_audit.json`; просте редагування згенерованого `.inp` буде перезаписано під час наступного повного запуску.

## Що потрібно погодити перед фізичними результатами

1. **Офіційний production Henry `.inp`.** Поточний файл прямо позначений provisional. У ньому Ebeam=10.60 GeV, SHMS electron P=6.70 GeV/c при 13°, HMS proton P=1.10 GeV/c при 32°. Радіація, втрати енергії та smearing вимкнені. Це не підтверджена конфігурація для вимірювання ширини піків. До 10.67 GeV самовільно нічого не переналаштовувалося.
2. **Знаки Jacobian/моделей.** У phi/rho/eta/eta-prime є від'ємні siglab та Weight при додатних sigcm. Функції використовують signed fixed-Mx Jacobian (`physics_phi.f:195,225`, `physics_eta.f:171,196`); rho-log безпосередньо друкує від'ємний Jacobian. Для omega Weight додатні, але sigcm від'ємна в 992/1000 записів. Тому додатна omega-гістограма сама по собі не доводить коректність моделі. Потрібно перевірити прийняту автором конвенцію знаків і застосовність моделей для цих налаштувань. `abs(Weight)` не застосовувався, фізичні формули не змінювалися.
3. **Two-pion phase space.** Upstream X-phase-space запускається, але не є підтвердженим генератором саме двопіонного кінцевого стану: код перевіряє поріг Mpi0 (`event.f:1094`, `physics_Xphasespace.f:66,81`), і в ROOT є missing mass від 0.136897 GeV. Потрібен конкретний двопіонний input/режим або погоджена реалізація. Заміна назви файла не вирішує цього.
4. **Нормалізація і визначення колонок отриманих ROOT.** Потрібні writer/опис `mx`, `t`, одиниць та зовнішні ваги чи scale factors перед порівнянням фізичних yields.

Перший пункт зустрічі — отримати omega-run і показати його — технічно виконано. Для остаточного фізичного завершення потрібні зазначені рішення керівника. Повідомлень від вашого імені не надсилалося.

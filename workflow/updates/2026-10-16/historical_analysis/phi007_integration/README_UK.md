# Інтеграція гілки Білла phi-007, 08.09.2026

Джерело: https://github.com/billlee77/simc_gfortran/tree/phi-007, commit `d778c567aee4c1a43230130fd83f1653fccfeb0f`.

Робочий Fortran checkout `/home/severus/eic/ACPOSTUMA_SIMC_GFORTRAN` оновлено fast-forward із b6025c4 до d778c56. SIMC успішно перекомпільовано та повністю злинковано в eic-shell. Існуючі шість ROOT-вибірок, inputs і графіки не перезапускали та не замінювали: вони залишаються результатами b6025c4.

У локальному checkout користувацького репозиторію `/home/severus/eic/publish/simc-acpostuma-roman-007` створено гілку `integrate-phi-007-20260908` і merge commit `d3df023398129ef6ab7831da36a46c86b91389c9`. Єдиний конфлікт був у README: збережено документацію обох версій та уточнено походження даних. Фізичні зміни й новий input збігаються з гілкою Білла. Заморожений snapshot та старий тег збережено. **Push на GitHub не виконувався.**

## Що додає гілка

`doing_Xphasespace = 3` — відносна вага фазового простору π⁺π⁻π⁰. Поріг у коді — 2 Mπ + Mπ0 = 414,11696 МеВ; нормування внутрішнього фазового простору взято при масі 1 ГеВ. Індивідуальні чотириімпульси піонів не генеруються. Це не абсолютний фізичний переріз. Моделі omega, phi, rho, eta, eta-prime не змінено, тому гілка не є виправленням їхніх попередніх від’ємних ваг.

Режим 1 збережено як legacy continuum. Старий коментар називав його 2pi; водночас код використовував однопіонний поріг, тому цей коментар сам по собі не є перевіркою повного двопіонного генератора.

## Тест нового input

Без змін запущено `infiles/xphase_threebody_smoke.inp`: пучок 4,70897 ГеВ, SOS/HMS, radiation/energy loss/smearing увімкнено. `ngen=-10000` означає 10 000 спроб. SIMC повідомив Ncontribute=91; перетворений ROOT містить 91 подію h10, 0 від’ємних і 0 нечислових Weight. Це smoke-test, не повна фізична валідація.

Усі результати: `smoke_20260908_075908/`, з input, log, .hist, ROOT та `root_validation.json`. Цей sample не змішувався з шістьма попередніми каналами при 10,60 ГеВ.

Upstream README згадує `tests/test_xphase_threebody.py`, але такого файла немає в дереві d778c56. Ці чисельні тести не запускалися.

Для запитаного порівняння omega використано попередній omega_Roman_007 при 10,60 ГеВ: нові Q²/W/x-графіки розташовані окремо в `../omega_kinematics_weight_comparison/`.

Повторити smoke-test із новою окремою папкою результатів:

```powershell
wsl.exe -d Ubuntu -- /mnt/d/EIC/eic-shell -- /mnt/d/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/test_phi007_smoke.sh
```

# ТЗ: Полное удаление асинхронного сэмплера из CoGAPS

**Ветка:** `132-uncertainty-improvements`
**Дата:** 2026-07-04

## 1. Цель и обоснование

Асинхронный (многопоточный) сэмплер нарушает баланс MCMC (detailed balance):
предложения генерируются и обрабатываются параллельно из очереди, что делает
цепь некорректной с точки зрения марковской теории. Принято решение
**полностью удалить** асинхронный путь, оставив только последовательный
`SingleThreadedGibbsSampler`.

Текущее состояние: async уже частично отключён — в `src/GapsRunner.cpp`
`#include` (стр. 8) и вызов диспетчера (стр. 72–73) закомментированы, поэтому
`chooseSampler` всегда падает в `SingleThreadedGibbsSampler`. Однако async-классы
всё ещё компилируются (перечислены в Makevars) и покрыты тестами, поэтому для
чистого удаления требуются правки ниже.

## 2. Область работ

Удаление затрагивает три группы:
- **A.** Файлы, существующие только ради async → удалить.
- **B.** Файлы, ссылающиеся на async → отредактировать.
- **C.** Вспомогательная параллельность (OpenMP), которая имеет смысл только при
  `maxThreads > 1` → снять как мёртвый код.

---

## 3. Файлы под УДАЛЕНИЕ (async-only)

| Файл | Причина |
|------|---------|
| `src/gibbs_sampler/AsynchronousGibbsSampler.h` | сам асинхронный сэмплер (`#pragma omp parallel for`, стр. 105) |
| `src/atomic/ProposalQueue.h` | очередь бесконфликтных `AtomicProposal` для параллельной обработки |
| `src/atomic/ProposalQueue.cpp` | реализация очереди |
| `src/atomic/ConcurrentAtomicDomain.h` | «OpenMP thread-safe» домен атомов |
| `src/atomic/ConcurrentAtomicDomain.cpp` | реализация |
| `src/atomic/ConcurrentAtom.h` | тип атома только для `ConcurrentAtomicDomain` |
| `src/atomic/ConcurrentAtom.cpp` | реализация |
| `src/cpp_tests/testConcurrentAtomicDomain.cpp` | тест исключительно `ConcurrentAtomicDomain` |

**Оставить (это sync-структуры, НЕ async):** `src/atomic/Atom.{h,cpp}`,
`src/atomic/AtomicDomain.{h,cpp}` — используются `SingleThreadedGibbsSampler`.
`DenseNormalModel` / `SparseNormalModel` (базовые DataModel) async-ссылок не
содержат и общие для обоих сэмплеров.

---

## 4. Файлы под ИЗМЕНЕНИЕ

### 4.1. `src/GapsRunner.cpp`
- Удалить закомментированный `#include "gibbs_sampler/AsynchronousGibbsSampler.h"` (стр. 8).
- В `chooseSampler` (стр. 65–78) убрать ветку `if (params.asynchronousUpdates)`
  (стр. 69–74) целиком; оставить прямой вызов `SingleThreadedGibbsSampler`.
- В `updateSampler` (стр. 201–222) убрать передачу `params.maxThreads` в
  `update()` / `sync()` (стр. 207, 210, 216, 219) — вызовы становятся однопоточными.
- Убрать `calculateNumberOfThreads` (стр. 352–364) и её вызов (стр. 441),
  а также `#include <omp.h>` (стр. 21–23), если больше не используется.
- Убрать присваивание `result.averageQueueLengthA/P = ...getAverageQueueLength()`
  (стр. 474–475) — см. §4.5.

### 4.2. `src/GapsParameters.h` — НЕ МЕНЯЕМ структуру

**Решение: структуру `GapsParameters` не трогаем.** Поля `bool asynchronousUpdates`
(стр. 63) и `unsigned maxThreads` (стр. 44) **остаются** — это сохраняет бинарный
формат сериализации/чекпоинтов без изменений. Поля становятся «мёртвыми»
(нигде не влияют на выбор сэмплера, т.к. диспетчер async удалён), но нейтрализованы:
- `maxThreads(1)` — уже дефолт (стр. 89), оставить.
- `asynchronousUpdates` — сменить дефолт `true` (стр. 108) → **`false`**
  (единственная правка в файле), чтобы флаг не вводил в заблуждение.

### 4.3. `src/GapsParameters.cpp`
- Печать `maxThreads` (стр. 17) и `asynchronousUpdates` (стр. 23) — **оставить**
  (поля существуют, вывод безвреден).

### 4.4. `src/Cogaps.cpp`
- Чтение R-параметров `nThreads`→`maxThreads` (стр. 85) и `asynchronousUpdates`
  (стр. 102) — жёстко зафиксировать в `1` / `false`, **не давая R переопределить**:
  ```cpp
  params.maxThreads = 1;                 // async removed — always single-threaded
  params.asynchronousUpdates = false;    // async removed
  ```
  (либо, если R перестанет передавать эти ключи, — просто убрать чтение и
  положиться на дефолты из §4.2). Согласовать с решением по R-API (§4.8).
- Убрать возврат в R `averageQueueLengthA/P` (стр. 177–178).
- `compiledWithOpenMPSupport_cpp()` (стр. 233–240) — удаляется в рамках варианта B
  (см. §5.3–5.4).

### 4.5. `src/GapsResult.h`
- Удалить поля `float averageQueueLengthA` / `averageQueueLengthP` (стр. 34–35) —
  это диагностика длины очереди async.

### 4.6. Сборка
- `src/Makevars` (стр. 4, `OBJECTS`): убрать `atomic/ConcurrentAtom.o`,
  `atomic/ConcurrentAtomicDomain.o`, `atomic/ProposalQueue.o`.
- `src/Makevars.win` (стр. 13, 15, 16): убрать те же три объектника.
- Проверить `src/Makevars.in`, если объектники прописаны и там.

### 4.7. Тесты C++
- `src/cpp_tests/testSamplerHighLevel.cpp`: убрать
  `#include "../gibbs_sampler/AsynchronousGibbsSampler.h"` (стр. 7) и
  `INIT_SAMPLER(..., AsynchronousGibbsSampler, ...)` (стр. 45–46, 61–62).
- `src/cpp_tests/testSerialization.cpp`: убрать `#include "../atomic/ProposalQueue.h"`
  (стр. 7) и тест-кейс «ProposalQueue Serialization» (стр. 298–347).

### 4.8. R-слой — deprecated-заглушки с предупреждением

**Решение:** аргументы `asynchronousUpdates` и `nThreads` **остаются** в сигнатурах
`CoGAPS` / `scCoGAPS` / `GWCoGAPS` (обратная совместимость со старыми скриптами и
Bioconductor), но **функционально игнорируются** — прогон всегда однопоточный
последовательный (C++ форсит `maxThreads=1` / `asynchronousUpdates=false`, §4.4).
`warning` выдаётся **только при попытке нестандартного значения**, чтобы дефолтные
вызовы молчали.

- `R/CoGAPS.R`:
  - Сменить дефолт `asynchronousUpdates=TRUE` → **`FALSE`** (стр. 92), иначе
    обычный `CoGAPS(data)` будет каждый раз ругаться. `nThreads=1` — дефолт оставить.
  - Добавить в начало тела предупреждение о депрекации:
    ```r
    if (!identical(nThreads, 1) || isTRUE(asynchronousUpdates))
        warning("'nThreads' and 'asynchronousUpdates' are deprecated and ignored; ",
                "CoGAPS now always runs single-threaded (async broke MCMC balance)")
    ```
  - Аргументы больше не влияют на результат; в C++ уходят форсированные `1`/`false`.
- `R/HelperFunctions.R`: предупреждение про `nThreads` (стр. 235–236) —
  убрать/подчинить новому deprecation-warning (не дублировать).
- `R/DistributedCogaps.R`: строки 32–33 (`asynchronousUpdates <- FALSE` /
  `nThreads <- 1`) — **оставить** (они выставляют игнорируемые поля в «тихие»
  значения; путь `callInternalCoGAPS` зовёт C++ напрямую, warning не триггерит).
- roxygen: пометить оба параметра как *deprecated* в `@param`; перегенерировать
  `.Rd` (man/) и `NAMESPACE` при необходимости.

> Позже (отдельным релизом) заглушки можно удалить полностью — тогда уйдут и
> аргументы, и форсирование в `Cogaps.cpp`, и правки тестов из §4.7-bis.

### 4.9. Влияние на распределённый режим (GWCoGAPS / scCoGAPS)

**Вывод: удаление async распределённый режим не ломает.** GWCoGAPS
(genome-wide) и scCoGAPS (single-cell) async **не используют** — каждый воркер
уже запускает чисто последовательный `SingleThreadedGibbsSampler`.

Их параллельность находится на **другом уровне** — над MCMC-цепью, а не внутри неё:
- `R/DistributedCogaps.R` бьёт данные на подмножества (`createSets`, стр. 56) —
  по строкам (гены) для genome-wide, по столбцам (клетки) для single-cell.
- `BiocParallel::bplapply(..., BPPARAM=...)` (стр. 60–61, 68–72, 97–101) запускает
  **независимый CoGAPS на каждом подмножестве в отдельном процессе-воркере**
  (по умолчанию `MulticoreParam`).
- Результаты сшиваются: `findConsensusMatrix` → второй проход с фиксированной
  матрицей → `stitchTogether`.

Это параллельность на уровне **процессов R** (независимые корректные
последовательные цепи на непересекающихся данных) — она **не нарушает** detailed
balance и никак не завязана на OpenMP/async в C++.

Сопоставление двух видов параллельности:

| | Async (удаляем) | Distributed (GWCoGAPS/scCoGAPS) |
|---|---|---|
| Уровень | внутри одной MCMC-цепи | над цепями, разные подмножества данных |
| Механизм | OpenMP (`ProposalQueue`, `ConcurrentAtomicDomain`) в C++ | `BiocParallel::bplapply` в R |
| Корректность MCMC | **нарушает** detailed balance | корректен |
| Внутренний сэмплер | `AsynchronousGibbsSampler` | `SingleThreadedGibbsSampler` |

**Единственная правка здесь:** в `callInternalCoGAPS` (`R/DistributedCogaps.R`)
удалить строки 32–33, которые принудительно выключают async:
```r
allParams$asynchronousUpdates <- FALSE
allParams$nThreads <- 1
```
После удаления полей `asynchronousUpdates`/`nThreads` (§4.8) эти строки будут
ссылаться на несуществующие параметры. Синхронизировать с решением из §4.8
(полное удаление vs. deprecated-заглушки). Логику разбиения на подмножества,
`bplapply`/`BPPARAM`, сопоставление и сшивание паттернов **не трогать**.

### 4.10. Автогенерируемое
- `src/RcppExports.cpp` — **не править вручную**, перегенерируется из R-слоя
  (`Rcpp::compileAttributes`).

---

## 5. Полное удаление OpenMP (вариант B)

**Решение: выпилить OpenMP целиком.** После удаления async единственные оставшиеся
прагмы — потокобезопасность без потоков (мёртвый оверхед), а GWCoGAPS/scCoGAPS от
OpenMP не зависят (их параллелизм — `BiocParallel` на уровне процессов R, §4.9).

### 5.1. Прагмы и потоковые вызовы в C++
- `src/gibbs_sampler/DenseNormalModel.cpp:26` — убрать
  `#pragma omp parallel for num_threads(nThreads)` в `sync()`, оставив обычный цикл.
- `src/data_structures/HybridVector.cpp` — убрать `#pragma omp atomic` (стр. 60,
  65, 77, 82) в `add()` / `set()` и комментарии «can be called from multiple
  concurrent OpenMP threads».
- `src/GapsRunner.cpp` — убрать `#include <omp.h>` (стр. 21–23), функцию
  `calculateNumberOfThreads` (`omp_get_max_threads`, стр. 352–364) и её вызов
  (стр. 441) (уже в §4.1).
- Убрать параметр `nThreads` из сигнатур `sync()`:
  `DenseNormalModel.h:61`, `DenseNormalModel.cpp:20`, `SparseNormalModel.h:27`,
  `SparseNormalModel.cpp:27`; и из `SingleThreadedGibbsSampler::update()`
  (`SingleThreadedGibbsSampler.h:45, 118`).

### 5.2. Инфраструктура сборки и макросы
- `src/utils/GlobalConfig.h` — убрать блок `#ifdef _OPENMP / #define __GAPS_OPENMP__`
  (стр. 12–14) и строку статуса «Compiled with OpenMP» в `configReport`
  (стр. 47–51, оставить только SIMD-репорт или заменить на «OpenMP: disabled»).
- `configure.ac` (стр. 56–64) — убрать `AC_ARG_ENABLE(openmp)`, `AX_OPENMP` и
  добавление `OPENMP_CXXFLAGS` в `GAPS_CXX_FLAGS`/`GAPS_LIBS`. **Перегенерировать
  `configure`** (`autoreconf`/`autoconf`) — файл `configure` автогенерируемый,
  вручную не редактировать (в нём соответствующие блоки на стр. ~653, 1288,
  2921–2933 уйдут при регенерации). Опционально удалить `m4/ax_openmp.m4`, если он
  больше нигде не нужен.
- `src/Makevars.win` — флагов OpenMP не содержит (PKG_CXXFLAGS/LIBS пустые), правок
  по флагам не требуется (только список объектников из §4.6).
- `src/Makevars` — генерируется из `Makevars.in` через `@GAPS_CXX_FLAGS@`; после
  правки `configure.ac` `-fopenmp` туда попадать перестанет автоматически.

### 5.3. R-экспорт статуса OpenMP — см. §5.4 (требует решения)
- `src/Cogaps.cpp::compiledWithOpenMPSupport_cpp()` (стр. 233–240) и
  `#ifdef __GAPS_OPENMP__` в нём.
- `src/RcppExports.cpp` — запись `_CoGAPS_compiledWithOpenMPSupport_cpp`
  (автоген, перегенерируется).
- `R/CoGAPS.R:34–37` — экспортируемая обёртка `compiledWithOpenMPSupport()`.
- `NAMESPACE:13` — `export(compiledWithOpenMPSupport)`.
- `R/CoGAPS.R:105–112` — блок `if (!compiledWithOpenMPSupport()) { ... }` —
  **удалить** в любом случае (заменяется deprecation-warning из §4.8).

> Примечание: `std::thread` / `std::async` / `std::mutex` / `std::atomic` в коде
> отсутствуют — вся параллельность была только на OpenMP-прагмах.

### 5.4. Публичная `compiledWithOpenMPSupport()` — РЕШЕНО: B1 (заглушка)
`compiledWithOpenMPSupport()` — **экспортируемая публичная функция** (в `NAMESPACE`,
с примером в доке). **Решение — B1:** оставить R-обёртку, чтобы не ломать публичный
API, но пусть возвращает `FALSE` (OpenMP больше нет — ответ честный):
```r
#' @return FALSE (OpenMP support removed; CoGAPS runs single-threaded)
compiledWithOpenMPSupport <- function() FALSE
```
- Убрать C++-часть: `compiledWithOpenMPSupport_cpp()` в `Cogaps.cpp` и запись
  `_CoGAPS_compiledWithOpenMPSupport_cpp` в `RcppExports.cpp` (перегенерируется) и
  `R/RcppExports.R:20–22`.
- `NAMESPACE:13` `export(compiledWithOpenMPSupport)` — **оставить**.
- Обновить роксиген `@return` (теперь всегда `FALSE`).

---

## 6. Метод «getAverageQueueLength»

`SingleThreadedGibbsSampler::getAverageQueueLength()` (`SingleThreadedGibbsSampler.h:92–96`)
возвращает `0.f` — это заглушка под интерфейс, нужный только для async-диагностики.
После удаления полей `averageQueueLengthA/P` (§4.5) и их вызовов в GapsRunner (§4.1)
метод можно удалить целиком.

---

## 7. Классы, удаляемые полностью (справка)

| Класс | Где определён | Кто использовал (всё уходит) |
|-------|---------------|------------------------------|
| `AsynchronousGibbsSampler<DataModel>` | `AsynchronousGibbsSampler.h` | `GapsRunner.cpp` (закомм.), `testSamplerHighLevel.cpp` |
| `ProposalQueue` + `struct AtomicProposal` | `ProposalQueue.{h,cpp}` | `AsynchronousGibbsSampler.h`, `ConcurrentAtomicDomain.h` (friend), `testSerialization.cpp` |
| `ConcurrentAtomicDomain` | `ConcurrentAtomicDomain.{h,cpp}` | `AsynchronousGibbsSampler.h`, `ProposalQueue.{h,cpp}`, `testConcurrentAtomicDomain.cpp` |
| `ConcurrentAtom` (+ neighborhood) | `ConcurrentAtom.{h,cpp}` | `ConcurrentAtomicDomain.{h,cpp}`, `ProposalQueue.h`, `AsynchronousGibbsSampler.h` (debug) |

---

## 8. Порядок выполнения (рекомендуемый)

1. Удалить файлы из §3.
2. Отредактировать `GapsRunner.cpp`, `GapsParameters.{h,cpp}`, `GapsResult.h`,
   `Cogaps.cpp` (§4.1–4.5).
3. Снять параллельность и параметры `nThreads` (§5, §6).
4. Обновить `Makevars` / `Makevars.win` / `Makevars.in` (§4.6).
5. Обновить C++-тесты (§4.7).
6. Обновить R-слой (включая `DistributedCogaps.R`, §4.8–4.9) и перегенерировать
   `RcppExports.cpp` (§4.10).
7. Пересобрать пакет и прогнать C++-тесты (`cpp_tests`) и R-тесты.

## 9. Критерии приёмки

- [ ] Проект собирается без ошибок и предупреждений о неразрешённых символах
      (`ConcurrentAtom*`, `ProposalQueue`, `AsynchronousGibbsSampler`).
- [ ] `grep -rn "Asynchronous\|ProposalQueue\|Concurrent" src/` не даёт совпадений
      (кроме, возможно, комментариев в истории).
- [ ] Все C++-тесты проходят; удалённые async-тесты отсутствуют в сборке.
- [ ] R-функции `CoGAPS/scCoGAPS/GWCoGAPS` работают; параметры
      `asynchronousUpdates`/`nThreads` либо удалены, либо помечены deprecated
      (по решению из §4.8).
- [ ] **Parity пройден** по протоколу §11 — результаты бит-в-бит совпадают с
      эталоном, снятым до удаления (async уже был отключён, а изменяемые прагмы
      FP-нейтральны, поэтому ожидается точное совпадение, не «в пределах допуска»).

## 10. Открытые вопросы

1. ~~**R-совместимость:**~~ **РЕШЕНО** — deprecated-заглушки с `warning` только при
   нестандартном значении (`nThreads != 1` или `asynchronousUpdates=TRUE`), см. §4.8.
2. ~~**OpenMP:**~~ **РЕШЕНО** — вариант **B** (полное удаление OpenMP, §5) +
   под-вариант **B1** для публичной `compiledWithOpenMPSupport()` (заглушка → `FALSE`,
   §5.4). Инфраструктура OpenMP **не нужна** для GWCoGAPS/scCoGAPS — их параллелизм
   на уровне процессов R (`BiocParallel`), а не C++-потоков; внутри воркеров
   однопоточность даже желательна (иначе N×T = oversubscription).
3. ~~**Формат результата:**~~ **РЕШЕНО** — `averageQueueLengthA/P` удаляются из
   `GapsResult` (async-диагностика, downstream не читает), см. §4.5.

**Все открытые вопросы закрыты.** Дополнительно зафиксировано:
- Ветка: работа продолжается в `132-uncertainty-improvements` (не отдельная).
- `GapsParameters` структуру не меняем: `maxThreads=1`, `asynchronousUpdates=false`.
- Проверка: обязательный parity-протокол, см. §11.

---

## 11. Протокол parity-проверки (обязательно)

Цель — доказать, что уборка не задела последовательный путь (страховка от
взаимодействия ошибок). Async уже отключён, а снимаемые прагмы FP-нейтральны
(`omp parallel for` при 1 потоке = тот же порядок цикла; `omp atomic` не меняет
значения) — поэтому ожидается **точное бит-в-бит совпадение**.

### 11.1. Снять эталон ДО правок
На текущем `HEAD` (до любых изменений) собрать пакет и прогнать матрицу конфигураций
с фиксированным seed, сохранив результаты в RDS:
```r
library(CoGAPS)
data(GIST)  # GIST.matrix
cfg <- function(tag, ...) {
    r <- CoGAPS(GIST.matrix, seed=42, nIterations=1000, messages=FALSE, ...)
    saveRDS(list(fl=r@featureLoadings, sf=r@sampleFactors,
                 chi=r@metadata$meanChiSq), sprintf("parity_before_%s.rds", tag))
}
cfg("dense")                                   # DenseNormalModel
cfg("sparse",  sparseOptimization=TRUE)        # SparseNormalModel
cfg("unc",     uncertainty=GIST.uncertainty)   # путь с матрицей неопределённости
# distributed (затрагиваются DistributedCogaps.R и R-аргументы):
rsc <- scCoGAPS(GIST.matrix, seed=42, nIterations=1000, messages=FALSE,
                nSets=2, BPPARAM=BiocParallel::SerialParam())
saveRDS(list(fl=rsc@featureLoadings, sf=rsc@sampleFactors), "parity_before_scc.rds")
```
> `SerialParam()` для распределённого прогона — чтобы результат был детерминирован
> и не зависел от числа воркеров/планировщика.

### 11.2. Снять то же ПОСЛЕ правок
Пересобрать пакет после всех изменений, прогнать те же вызовы, сохранить в
`parity_after_*.rds`.

### 11.3. Сверить
```r
for (tag in c("dense","sparse","unc","scc")) {
    a <- readRDS(sprintf("parity_before_%s.rds", tag))
    b <- readRDS(sprintf("parity_after_%s.rds",  tag))
    stopifnot(identical(a$fl, b$fl), identical(a$sf, b$sf))
}
```
Критерий — `identical()` (точное совпадение). Любое расхождение = уборка задела
численный путь → разбираться до коммита в PR.

### 11.4. Плюс стандартные проверки
- `cpp_tests` (Catch) — все проходят.
- `R CMD check` / `devtools::test()` — с учётом правок `test_seed_consistency.R`,
  `test_top_level.R`, `inst/scripts/debugRuns.R` (§4.7-related).
- Проверить оба build-пути, где расходилось раньше (см. issues 9–10):
  `devtools::install_local()` (release, `-O2`), не только `load_all()`.

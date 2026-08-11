# Async sampler removal — implementation report

**Branch:** `132-uncertainty-improvements`
**Date:** 2026-07-06
**Spec:** `remove-async-plan-eng.md` (this report is the retrospective record of
executing that spec)

---

## 1. Summary

The asynchronous (OpenMP multi-threaded) Gibbs sampler was **completely removed**
from CoGAPS. It broke MCMC detailed balance: proposals were generated and applied
in parallel out of a `ProposalQueue`, so the Markov chain was not sampling from the
intended posterior. CoGAPS now always runs the sequential
`SingleThreadedGibbsSampler`.

The change is **behaviour-preserving for all users**: the async dispatch had
already been commented out, so every run already fell through to the sequential
sampler. A four-configuration parity check (below) confirms bit-for-bit identical
results before and after removal.

Scope: **8 files deleted, 28 modified** (1589 insertions / 369 deletions; the large
`configure` insertion count is from regenerating it — see §7).

---

## 2. Rationale

- The async sampler's `ProposalQueue` batched conflict-free proposals and applied
  them in parallel (`#pragma omp parallel for`). Overlapping accept/reject decisions
  within a batch violate detailed balance.
- The path was already dead: `chooseSampler` in `GapsRunner.cpp` had its async
  branch commented out, always constructing `SingleThreadedGibbsSampler`.
- The async classes still compiled into the package and were exercised by tests, so
  a clean removal was still required.
- OpenMP existed **only** to support the async sampler. Distributed CoGAPS
  (GWCoGAPS / scCoGAPS) parallelises at the R process level via `BiocParallel`, not
  C++ threads, and is unaffected.

---

## 3. Files deleted (async-only)

| File | Role |
|------|------|
| `src/gibbs_sampler/AsynchronousGibbsSampler.h` | the async sampler (was already non-compiling in isolation) |
| `src/atomic/ProposalQueue.{h,cpp}` | conflict-free proposal batching for parallel apply |
| `src/atomic/ConcurrentAtomicDomain.{h,cpp}` | OpenMP-thread-safe atom domain |
| `src/atomic/ConcurrentAtom.{h,cpp}` | atom type for the concurrent domain |
| `src/cpp_tests/testConcurrentAtomicDomain.cpp` | test for the concurrent domain |

Kept (sequential structures, shared by the surviving sampler): `src/atomic/Atom.*`,
`src/atomic/AtomicDomain.*`, `DenseNormalModel`, `SparseNormalModel`.

---

## 4. C++ changes

### 4.1 Core driver — `src/GapsRunner.cpp`
- Removed the async `#include`, `#include <omp.h>`, and the `if (asynchronousUpdates)`
  branch in `chooseSampler`.
- `updateSampler` no longer passes `maxThreads` into `update()`/`sync()`.
- Deleted `calculateNumberOfThreads()` (used `omp_get_max_threads()`) and its call.
- Dropped `result.averageQueueLengthA/P = ...getAverageQueueLength()`.

### 4.2 `GapsParameters` — fields removed

> **Update (2026-07-07):** the fields `bool asynchronousUpdates` and
> `unsigned maxThreads` were initially **kept** on the belief they were part of the
> checkpoint serialization format. That was wrong: `GapsParameters::operator<</>>`
> serialize only 11 fields (`seed, nGenes, nSamples, nPatterns, nIterations, alphaA,
> alphaP, maxGibbsMassA, maxGibbsMassP, useSparseOptimization, checkpointInterval`)
> and never touched `maxThreads`/`asynchronousUpdates`. So the format was never
> affected, and both fields (plus their `print()` lines and the forced assignments
> in `Cogaps.cpp`) were subsequently **removed** entirely. The R arguments
> `nThreads`/`asynchronousUpdates` remain as deprecated no-ops (§4.8). Verified by
> the new `[serialization][gapsparameters]` round-trip test and unchanged parity.

### 4.3 Samplers / data structures (§ "variant B" OpenMP strip)
- `SingleThreadedGibbsSampler`: removed `getAverageQueueLength()`; `update(nSteps,
  nThreads)` → `update(nSteps)`.
- `DenseNormalModel::sync` / `SparseNormalModel::sync`: dropped the `nThreads`
  parameter; removed the `#pragma omp parallel for` from the dense AP-transpose loop.
- `HybridVector::add`/`set`: removed the `#pragma omp atomic` directives (thread
  safety with no threads = dead overhead).

### 4.4 `GapsResult.h`
- Removed `averageQueueLengthA` / `averageQueueLengthP` (async queue diagnostics,
  not read by any downstream R code).

### 4.5 OpenMP infrastructure (variant B — full removal)
- `src/utils/GlobalConfig.h`: removed the `#ifdef _OPENMP → #define __GAPS_OPENMP__`
  block and the "Compiled with OpenMP" status line (now "OpenMP: disabled").
- `src/Cogaps.cpp`: removed `compiledWithOpenMPSupport_cpp()`.
- Build files: see §7.

---

## 5. R layer

- **Deprecated stubs (variant chosen):** `nThreads` and `asynchronousUpdates` remain
  arguments of `CoGAPS` / `scCoGAPS` / `GWCoGAPS` for backward compatibility but are
  ignored. A `warning` fires **only** on a non-default value
  (`nThreads != 1 || isTRUE(asynchronousUpdates)`), so default calls stay silent.
  `CoGAPS`'s `asynchronousUpdates` default flipped `TRUE → FALSE`.
- **`compiledWithOpenMPSupport()` (variant B1):** the exported R function is kept
  (public API preserved) but now returns `FALSE` directly; the C++ `_cpp` backend
  and its RcppExports entry were removed.
- `R/HelperFunctions.R`: removed the now-dead "can't run multi-threaded and
  distributed" warning.
- `R/DistributedCogaps.R`: `callInternalCoGAPS` still sets `asynchronousUpdates <-
  FALSE` / `nThreads <- 1` on the inner param list — kept (harmless, sets the
  ignored fields to quiet values).
- roxygen `@param` / `@return` updated; `man/*.Rd` regenerated.
- `RcppExports.{cpp,R}` regenerated (`Rcpp::compileAttributes`).

Distributed mode (GWCoGAPS / scCoGAPS) is unaffected — it parallelises over data
subsets with `BiocParallel::bplapply`, each worker running the sequential sampler.

---

## 6. Tests

- `testConcurrentAtomicDomain.cpp` deleted.
- `testSamplerHighLevel.cpp`: async sampler variants removed (only the two
  `SingleThreadedGibbsSampler` variants remain; the "Sampler Update" case is a
  construction smoke test — calling `update()` there would need `sync()` +
  `extraInitialization()` first, or `mOtherMatrix` is NULL → segfault).
- `testSerialization.cpp`: the async-only "ProposalQueue Serialization" case (which
  also used a stale API) removed; the other serialization cases stay.
- `testDenseGibbsSampler.cpp`: `update(100, 1)` callers updated to `update(100)`.
- R tests: `test_seed_consistency.R` and `test_top_level.R` had their `nThreads`
  variants removed (they tested multi-thread determinism that no longer exists);
  seed-consistency coverage for the dense and sparse samplers is retained.
- `inst/scripts/debugRuns.R`: `asynchronousUpdates` / `nThreads` invocations removed.

C++ Catch suite: **46/46 test cases pass** after removal (unchanged count — the
async test was never in the build list).

---

## 7. Build system — `configure` regenerated (not hand-patched)

`configure.ac` edits: removed the `AC_ARG_ENABLE(openmp)` / `AX_OPENMP` block and its
`OPENMP_CXXFLAGS` injection; removed the three async object files
(`atomic/ConcurrentAtom.o`, `atomic/ConcurrentAtomicDomain.o`, `atomic/ProposalQueue.o`).
`src/Makevars.win` had the same three object files removed.

The generated `configure` was **regenerated from `configure.ac`**, not hand-edited.
Because `automake`/`aclocal` is not installed (only `autoconf` + `autoconf-archive`),
the standard `autoreconf` fails on `aclocal`. Worked around it by running `autoconf`
against a temporary `aclocal.m4` that `m4_include`s the three archive macros
(`ax_openmp.m4`, `ax_compiler_vendor.m4`, `ax_compiler_version.m4`); the temp
`aclocal.m4` is **not** left in the repo. The regenerated `configure` is cleaner than
a hand patch: zero `openmp`/`OPENMP` references (including the residual
`enable_openmp` boilerplate), zero async objects, and the `AX_COMPILER_*` macros now
expand properly (they were literal no-ops in the previously committed `configure`).
This accounts for the large `configure` diff. Build validated end-to-end (RC 0).

> For a standard `autoreconf` workflow later: `brew install automake`.

---

## 8. Verification

All run against the actually rebuilt, async-removed package (an early check
accidentally ran against a stale install because a build had silently failed on a
`update(x,y)` caller — caught and fixed).

- **Build:** `R CMD INSTALL --preclean` → RC 0, no unresolved symbols.
- **Symbol audit:** `grep -rn "Asynchronous\|ProposalQueue\|Concurrent" src/` → empty;
  no `__GAPS_OPENMP__` / `#pragma omp` / `omp.h` remain.
- **C++ tests:** 46/46 Catch cases pass (823904 assertions).
- **Parity (§11 of the spec):** `CoGAPS` on `GIST.matrix`, `seed=42`,
  `nIterations=1000`, compared before/after with `identical()`:

  | Config | featureLoadings | sampleFactors |
  |--------|:---:|:---:|
  | dense | identical | identical |
  | sparse (`sparseOptimization=TRUE`) | identical | identical |
  | uncertainty (`GIST.uncertainty`) | identical | identical |
  | genome-wide distributed (`nSets=2`, `SerialParam`) | identical | identical |

  (CoGAPS was first confirmed bit-reproducible for a fixed seed, so `identical()`
  is a valid criterion.)
- **Deprecation behaviour:** `compiledWithOpenMPSupport()` → `FALSE`; a default
  `CoGAPS()` call emits 0 deprecation warnings; `CoGAPS(..., nThreads=4)` emits
  exactly 1 and still returns a `CogapsResult`.
- **Affected R tests:** `test_seed_consistency.R` (4 passed) and `test_top_level.R`
  (36 passed), 0 warnings, after the `nThreads` cleanup.

---

## 9. Related bug found and fixed during this work

Enabling `testSamplerHighLevel` in the build surfaced a **real, latent product bug**
(unrelated to async): `gaps::min/max(const SparseVector&)` dereferenced the first
element before checking `atEnd()`, segfaulting on an empty (all-zero) sparse vector —
reachable on any `sparseOptimization=TRUE` run whose data has an all-zero gene or
sample (`SparseNormalModel` calls `gaps::max(mDMatrix)` at construction). Root-caused
with a standalone AddressSanitizer harness. Fixed (empty → `0.f`) with a regression
test. This is documented separately as **issue 12** in
`132-LLM-assisted-solved-issues.md` and was landed as its own commit
(`dab57f2b`), together with the test-build enablement commit (`989db34f`) that
predates the async removal proper.

---

## 10. Notes / caveats

- `man/*.Rd` were regenerated with the local **roxygen2 8.0.0**, while `DESCRIPTION`
  pins `RoxygenNote: 7.3.2`. The doc *content* is correct; `DESCRIPTION` and
  `NAMESPACE` were reverted to avoid roxygen version churn and to keep
  `#import(fgsea)` commented out. Regenerating with 7.3.2 later is optional (format
  only).
- `GapsParameters` deliberately keeps the two inert fields to preserve checkpoint
  binary compatibility. A future release could drop the deprecated R arguments and
  those fields together (a breaking change).
- `std::thread` / `std::async` / `std::mutex` / `std::atomic` were never used — all
  parallelism was OpenMP pragmas, now gone.

---

## 11. Changeset (uncommitted at time of writing)

Deleted (8): the five async source pairs/headers listed in §3.

Modified (28): `src/GapsRunner.cpp`, `src/GapsParameters.h`, `src/GapsResult.h`,
`src/Cogaps.cpp`, `src/RcppExports.cpp`, `src/Makevars.win`,
`src/gibbs_sampler/{SingleThreadedGibbsSampler.h, DenseNormalModel.{h,cpp},
SparseNormalModel.{h,cpp}}`, `src/data_structures/HybridVector.cpp`,
`src/utils/GlobalConfig.h`,
`src/cpp_tests/{testDenseGibbsSampler.cpp, testSamplerHighLevel.cpp,
testSerialization.cpp}`, `configure`, `configure.ac`,
`R/{CoGAPS.R, HelperFunctions.R, RcppExports.R}`, `man/*.Rd` (4),
`inst/scripts/debugRuns.R`, `tests/testthat/{test_seed_consistency.R,
test_top_level.R}`.

Ready to commit.

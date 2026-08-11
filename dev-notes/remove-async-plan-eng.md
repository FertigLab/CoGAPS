# Spec: Complete removal of the asynchronous sampler from CoGAPS

**Branch:** `132-uncertainty-improvements`
**Date:** 2026-07-04

## 1. Goal and rationale

The asynchronous (multi-threaded) sampler breaks MCMC detailed balance:
proposals are generated and processed in parallel from a queue, which makes the
chain incorrect from the standpoint of Markov-chain theory. The decision is to
**completely remove** the asynchronous path, leaving only the sequential
`SingleThreadedGibbsSampler`.

Current state: async is already partially disabled — in `src/GapsRunner.cpp` the
`#include` (line 8) and the dispatcher call (lines 72–73) are commented out, so
`chooseSampler` always falls through to `SingleThreadedGibbsSampler`. However, the
async classes are still compiled (listed in Makevars) and still covered by tests,
so a clean removal requires the edits below.

## 2. Scope of work

The removal touches three groups:
- **A.** Files that exist only for async → delete.
- **B.** Files that reference async → edit.
- **C.** Auxiliary parallelism (OpenMP) that only makes sense when
  `maxThreads > 1` → strip as dead code.

---

## 3. Files to DELETE (async-only)

| File | Reason |
|------|--------|
| `src/gibbs_sampler/AsynchronousGibbsSampler.h` | the async sampler itself (`#pragma omp parallel for`, line 105) |
| `src/atomic/ProposalQueue.h` | queue of conflict-free `AtomicProposal`s for parallel processing |
| `src/atomic/ProposalQueue.cpp` | queue implementation |
| `src/atomic/ConcurrentAtomicDomain.h` | "OpenMP thread-safe" atom domain |
| `src/atomic/ConcurrentAtomicDomain.cpp` | implementation |
| `src/atomic/ConcurrentAtom.h` | atom type used only by `ConcurrentAtomicDomain` |
| `src/atomic/ConcurrentAtom.cpp` | implementation |
| `src/cpp_tests/testConcurrentAtomicDomain.cpp` | test exclusively for `ConcurrentAtomicDomain` |

**Keep (these are sync structures, NOT async):** `src/atomic/Atom.{h,cpp}`,
`src/atomic/AtomicDomain.{h,cpp}` — used by `SingleThreadedGibbsSampler`.
`DenseNormalModel` / `SparseNormalModel` (the base DataModels) contain no async
references and are shared by both samplers.

---

## 4. Files to EDIT

### 4.1. `src/GapsRunner.cpp`
- Remove the commented `#include "gibbs_sampler/AsynchronousGibbsSampler.h"` (line 8).
- In `chooseSampler` (lines 65–78) remove the `if (params.asynchronousUpdates)`
  branch (lines 69–74) entirely; keep the direct call to `SingleThreadedGibbsSampler`.
- In `updateSampler` (lines 201–222) remove passing `params.maxThreads` into
  `update()` / `sync()` (lines 207, 210, 216, 219) — the calls become single-threaded.
- Remove `calculateNumberOfThreads` (lines 352–364) and its call (line 441),
  as well as `#include <omp.h>` (lines 21–23), if no longer used.
- Remove the assignment `result.averageQueueLengthA/P = ...getAverageQueueLength()`
  (lines 474–475) — see §4.5.

### 4.2. `src/GapsParameters.h` — DO NOT change the struct

**Decision: leave the `GapsParameters` struct untouched.** The fields
`bool asynchronousUpdates` (line 63) and `unsigned maxThreads` (line 44)
**remain** — this preserves the binary serialization/checkpoint format unchanged.
The fields become "dead" (they no longer affect sampler selection, since the async
dispatcher is removed) but neutralized:
- `maxThreads(1)` — already the default (line 89), keep it.
- `asynchronousUpdates` — change the default `true` (line 108) → **`false`**
  (the only edit in this file), so the flag is not misleading.

### 4.3. `src/GapsParameters.cpp`
- The prints of `maxThreads` (line 17) and `asynchronousUpdates` (line 23) —
  **keep** (the fields still exist, the output is harmless).

### 4.4. `src/Cogaps.cpp`
- Reading the R parameters `nThreads`→`maxThreads` (line 85) and
  `asynchronousUpdates` (line 102) — hard-fix to `1` / `false`, **not letting R
  override them**:
  ```cpp
  params.maxThreads = 1;                 // async removed — always single-threaded
  params.asynchronousUpdates = false;    // async removed
  ```
  (or, if R stops passing these keys — just drop the reads and rely on the defaults
  from §4.2). Reconcile with the R-API decision (§4.8).
- Remove the return to R of `averageQueueLengthA/P` (lines 177–178).
- `compiledWithOpenMPSupport_cpp()` (lines 233–240) — removed as part of variant B
  (see §5.3–5.4).

### 4.5. `src/GapsResult.h`
- Remove the fields `float averageQueueLengthA` / `averageQueueLengthP`
  (lines 34–35) — this is async queue-length diagnostics.

### 4.6. Build
- `src/Makevars` (line 4, `OBJECTS`): remove `atomic/ConcurrentAtom.o`,
  `atomic/ConcurrentAtomicDomain.o`, `atomic/ProposalQueue.o`.
- `src/Makevars.win` (lines 13, 15, 16): remove the same three object files.
- Check `src/Makevars.in` in case the object files are listed there too.

### 4.7. C++ tests
- `src/cpp_tests/testSamplerHighLevel.cpp`: remove
  `#include "../gibbs_sampler/AsynchronousGibbsSampler.h"` (line 7) and
  `INIT_SAMPLER(..., AsynchronousGibbsSampler, ...)` (lines 45–46, 61–62).
- `src/cpp_tests/testSerialization.cpp`: remove `#include "../atomic/ProposalQueue.h"`
  (line 7) and the "ProposalQueue Serialization" test case (lines 298–347).

### 4.8. R layer — deprecated stubs with a warning

**Decision:** the arguments `asynchronousUpdates` and `nThreads` **remain** in the
signatures of `CoGAPS` / `scCoGAPS` / `GWCoGAPS` (backward compatibility with old
scripts and Bioconductor), but are **functionally ignored** — a run is always
single-threaded sequential (C++ forces `maxThreads=1` / `asynchronousUpdates=false`,
§4.4). A `warning` is emitted **only on an attempt to use a non-default value**, so
that default calls stay silent.

- `R/CoGAPS.R`:
  - Change the default `asynchronousUpdates=TRUE` → **`FALSE`** (line 92), otherwise
    a plain `CoGAPS(data)` would warn every time. `nThreads=1` — keep the default.
  - Add a deprecation warning at the top of the body:
    ```r
    if (!identical(nThreads, 1) || isTRUE(asynchronousUpdates))
        warning("'nThreads' and 'asynchronousUpdates' are deprecated and ignored; ",
                "CoGAPS now always runs single-threaded (async broke MCMC balance)")
    ```
  - The arguments no longer affect the result; the forced `1`/`false` go to C++.
- `R/HelperFunctions.R`: the warning about `nThreads` (lines 235–236) —
  remove/subsume it into the new deprecation warning (do not duplicate).
- `R/DistributedCogaps.R`: lines 32–33 (`asynchronousUpdates <- FALSE` /
  `nThreads <- 1`) — **keep** (they set the ignored fields to the "quiet" values;
  the `callInternalCoGAPS` path calls C++ directly and does not trigger the warning).
- roxygen: mark both parameters as *deprecated* in `@param`; regenerate the
  `.Rd` (man/) and `NAMESPACE` if needed.

> Later (in a separate release) the stubs can be removed entirely — then the
> arguments, the forcing in `Cogaps.cpp`, and the test edits go away.

### 4.9. Impact on the distributed mode (GWCoGAPS / scCoGAPS)

**Conclusion: removing async does not break the distributed mode.** GWCoGAPS
(genome-wide) and scCoGAPS (single-cell) **do not use** async — each worker already
runs the purely sequential `SingleThreadedGibbsSampler`.

Their parallelism lives at a **different level** — above the MCMC chain, not inside it:
- `R/DistributedCogaps.R` splits the data into subsets (`createSets`, line 56) —
  by rows (genes) for genome-wide, by columns (cells) for single-cell.
- `BiocParallel::bplapply(..., BPPARAM=...)` (lines 60–61, 68–72, 97–101) launches
  an **independent CoGAPS on each subset in a separate worker process**
  (default `MulticoreParam`).
- The results are stitched: `findConsensusMatrix` → a second pass with a fixed
  matrix → `stitchTogether`.

This is parallelism at the level of **R processes** (independent, correct sequential
chains over disjoint data subsets) — it **does not break** detailed balance and is
in no way tied to OpenMP/async in C++.

Comparison of the two kinds of parallelism:

| | Async (removed) | Distributed (GWCoGAPS/scCoGAPS) |
|---|---|---|
| Level | inside a single MCMC chain | across chains, different data subsets |
| Mechanism | OpenMP (`ProposalQueue`, `ConcurrentAtomicDomain`) in C++ | `BiocParallel::bplapply` in R |
| MCMC correctness | **breaks** detailed balance | correct |
| Inner sampler | `AsynchronousGibbsSampler` | `SingleThreadedGibbsSampler` |

**The only edit here:** in `callInternalCoGAPS` (`R/DistributedCogaps.R`) remove
lines 32–33 that forcibly disable async:
```r
allParams$asynchronousUpdates <- FALSE
allParams$nThreads <- 1
```
After removing the `asynchronousUpdates`/`nThreads` fields (§4.8) these lines would
reference non-existent parameters. Reconcile with the decision in §4.8
(full removal vs. deprecated stubs). Do **not** touch the subset-splitting logic,
`bplapply`/`BPPARAM`, or the pattern matching/stitching.

### 4.10. Auto-generated
- `src/RcppExports.cpp` — **do not edit by hand**, it is regenerated from the R
  layer (`Rcpp::compileAttributes`).

---

## 5. Complete OpenMP removal (variant B)

**Decision: strip OpenMP entirely.** After async removal the only remaining pragmas
are thread safety without threads (dead overhead), and GWCoGAPS/scCoGAPS do not
depend on OpenMP (their parallelism is `BiocParallel` at the R process level, §4.9).

### 5.1. Pragmas and threading calls in C++
- `src/gibbs_sampler/DenseNormalModel.cpp:26` — remove
  `#pragma omp parallel for num_threads(nThreads)` in `sync()`, keeping a plain loop.
- `src/data_structures/HybridVector.cpp` — remove `#pragma omp atomic` (lines 60,
  65, 77, 82) in `add()` / `set()` and the "can be called from multiple concurrent
  OpenMP threads" comments.
- `src/GapsRunner.cpp` — remove `#include <omp.h>` (lines 21–23), the function
  `calculateNumberOfThreads` (`omp_get_max_threads`, lines 352–364) and its call
  (line 441) (already in §4.1).
- Remove the `nThreads` parameter from the `sync()` signatures:
  `DenseNormalModel.h:61`, `DenseNormalModel.cpp:20`, `SparseNormalModel.h:27`,
  `SparseNormalModel.cpp:27`; and from `SingleThreadedGibbsSampler::update()`
  (`SingleThreadedGibbsSampler.h:45, 118`).

### 5.2. Build infrastructure and macros
- `src/utils/GlobalConfig.h` — remove the `#ifdef _OPENMP / #define __GAPS_OPENMP__`
  block (lines 12–14) and the "Compiled with OpenMP" status line in `configReport`
  (lines 47–51, keep only the SIMD report or replace with "OpenMP: disabled").
- `configure.ac` (lines 56–64) — remove `AC_ARG_ENABLE(openmp)`, `AX_OPENMP` and the
  addition of `OPENMP_CXXFLAGS` to `GAPS_CXX_FLAGS`/`GAPS_LIBS`. **Regenerate
  `configure`** (`autoreconf`/`autoconf`) — the `configure` file is auto-generated,
  do not edit by hand (the corresponding blocks around lines ~653, 1288, 2921–2933
  disappear on regeneration). Optionally remove `m4/ax_openmp.m4` if it is no longer
  needed anywhere.
- `src/Makevars.win` — contains no OpenMP flags (PKG_CXXFLAGS/LIBS are empty), no
  flag edits needed (only the object list from §4.6).
- `src/Makevars` — generated from `Makevars.in` via `@GAPS_CXX_FLAGS@`; after editing
  `configure.ac`, `-fopenmp` stops flowing in automatically.

### 5.3. R export of OpenMP status — see §5.4 (needs a decision)
- `src/Cogaps.cpp::compiledWithOpenMPSupport_cpp()` (lines 233–240) and the
  `#ifdef __GAPS_OPENMP__` inside it.
- `src/RcppExports.cpp` — the `_CoGAPS_compiledWithOpenMPSupport_cpp` entry
  (auto-generated, regenerated).
- `R/CoGAPS.R:34–37` — the exported wrapper `compiledWithOpenMPSupport()`.
- `NAMESPACE:13` — `export(compiledWithOpenMPSupport)`.
- `R/CoGAPS.R:105–112` — the `if (!compiledWithOpenMPSupport()) { ... }` block —
  **remove** in any case (replaced by the deprecation warning from §4.8).

> Note: `std::thread` / `std::async` / `std::mutex` / `std::atomic` are absent from
> the code — all parallelism was solely on OpenMP pragmas.

### 5.4. Public `compiledWithOpenMPSupport()` — DECIDED: B1 (stub)
`compiledWithOpenMPSupport()` is an **exported public function** (in `NAMESPACE`,
with an example in the docs). **Decision — B1:** keep the R wrapper so as not to
break the public API, but have it return `FALSE` (there is no OpenMP anymore — an
honest answer):
```r
#' @return FALSE (OpenMP support removed; CoGAPS runs single-threaded)
compiledWithOpenMPSupport <- function() FALSE
```
- Remove the C++ side: `compiledWithOpenMPSupport_cpp()` in `Cogaps.cpp` and the
  `_CoGAPS_compiledWithOpenMPSupport_cpp` entry in `RcppExports.cpp` (regenerated) and
  `R/RcppExports.R:20–22`.
- `NAMESPACE:13` `export(compiledWithOpenMPSupport)` — **keep**.
- Update the roxygen `@return` (now always `FALSE`).

---

## 6. The `getAverageQueueLength` method

`SingleThreadedGibbsSampler::getAverageQueueLength()` (`SingleThreadedGibbsSampler.h:92–96`)
returns `0.f` — it is a stub for an interface needed only for async diagnostics.
After removing the `averageQueueLengthA/P` fields (§4.5) and their calls in
GapsRunner (§4.1), the method can be removed entirely.

---

## 7. Classes removed entirely (reference)

| Class | Where defined | Who used it (all goes away) |
|-------|---------------|-----------------------------|
| `AsynchronousGibbsSampler<DataModel>` | `AsynchronousGibbsSampler.h` | `GapsRunner.cpp` (commented), `testSamplerHighLevel.cpp` |
| `ProposalQueue` + `struct AtomicProposal` | `ProposalQueue.{h,cpp}` | `AsynchronousGibbsSampler.h`, `ConcurrentAtomicDomain.h` (friend), `testSerialization.cpp` |
| `ConcurrentAtomicDomain` | `ConcurrentAtomicDomain.{h,cpp}` | `AsynchronousGibbsSampler.h`, `ProposalQueue.{h,cpp}`, `testConcurrentAtomicDomain.cpp` |
| `ConcurrentAtom` (+ neighborhood) | `ConcurrentAtom.{h,cpp}` | `ConcurrentAtomicDomain.{h,cpp}`, `ProposalQueue.h`, `AsynchronousGibbsSampler.h` (debug) |

---

## 8. Recommended order of execution

1. Delete the files from §3.
2. Edit `GapsRunner.cpp`, `GapsParameters.{h,cpp}`, `GapsResult.h`,
   `Cogaps.cpp` (§4.1–4.5).
3. Strip the parallelism and `nThreads` parameters (§5, §6).
4. Update `Makevars` / `Makevars.win` / `Makevars.in` (§4.6).
5. Update the C++ tests (§4.7).
6. Update the R layer (including `DistributedCogaps.R`, §4.8–4.9) and regenerate
   `RcppExports.cpp` (§4.10).
7. Rebuild the package and run the C++ tests (`cpp_tests`) and the R tests.

## 9. Acceptance criteria

- [ ] The project builds without errors or warnings about unresolved symbols
      (`ConcurrentAtom*`, `ProposalQueue`, `AsynchronousGibbsSampler`).
- [ ] `grep -rn "Asynchronous\|ProposalQueue\|Concurrent" src/` returns no matches
      (except possibly comments in history).
- [ ] All C++ tests pass; the removed async tests are absent from the build.
- [ ] The R functions `CoGAPS/scCoGAPS/GWCoGAPS` work; the `asynchronousUpdates`/
      `nThreads` parameters are either removed or marked deprecated
      (per the decision in §4.8).
- [ ] **Parity passed** per the §11 protocol — results match the baseline captured
      before removal bit-for-bit (async was already disabled, and the pragmas being
      removed are FP-neutral, so exact equality is expected, not "within tolerance").

## 10. Open questions

1. ~~**R compatibility:**~~ **DECIDED** — deprecated stubs with a `warning` only for
   a non-default value (`nThreads != 1` or `asynchronousUpdates=TRUE`), see §4.8.
2. ~~**OpenMP:**~~ **DECIDED** — variant **B** (complete OpenMP removal, §5) +
   sub-variant **B1** for the public `compiledWithOpenMPSupport()` (stub → `FALSE`,
   §5.4). OpenMP infrastructure is **not needed** for GWCoGAPS/scCoGAPS — their
   parallelism is at the R process level (`BiocParallel`), not C++ threads; inside
   workers single-threading is even desirable (otherwise N×T = oversubscription).
3. ~~**Result format:**~~ **DECIDED** — `averageQueueLengthA/P` are removed from
   `GapsResult` (async diagnostics, downstream does not read them), see §4.5.

**All open questions are closed.** Additionally recorded:
- Branch: work continues in `132-uncertainty-improvements` (not a separate branch).
- `GapsParameters` struct is left unchanged: `maxThreads=1`, `asynchronousUpdates=false`.
- Verification: mandatory parity protocol, see §11.

---

## 11. Parity-check protocol (mandatory)

Goal — prove that the cleanup did not touch the sequential path (insurance against
interacting bugs). Async is already disabled, and the pragmas being removed are
FP-neutral (`omp parallel for` with 1 thread = the same loop order; `omp atomic`
does not change values) — so **exact bit-for-bit equality is expected**.

### 11.1. Capture the baseline BEFORE the edits
On the current `HEAD` (before any changes) build the package and run a matrix of
configurations with a fixed seed, saving the results to RDS:
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
cfg("unc",     uncertainty=GIST.uncertainty)   # uncertainty-matrix path
# distributed (DistributedCogaps.R and the R arguments are affected):
rsc <- scCoGAPS(GIST.matrix, seed=42, nIterations=1000, messages=FALSE,
                nSets=2, BPPARAM=BiocParallel::SerialParam())
saveRDS(list(fl=rsc@featureLoadings, sf=rsc@sampleFactors), "parity_before_scc.rds")
```
> `SerialParam()` for the distributed run — so the result is deterministic and does
> not depend on the number of workers/scheduler.

### 11.2. Capture the same AFTER the edits
Rebuild the package after all changes, run the same calls, save to
`parity_after_*.rds`.

### 11.3. Compare
```r
for (tag in c("dense","sparse","unc","scc")) {
    a <- readRDS(sprintf("parity_before_%s.rds", tag))
    b <- readRDS(sprintf("parity_after_%s.rds",  tag))
    stopifnot(identical(a$fl, b$fl), identical(a$sf, b$sf))
}
```
Criterion — `identical()` (exact match). Any discrepancy = the cleanup touched the
numerical path → investigate before committing to the PR.

### 11.4. Plus standard checks
- `cpp_tests` (Catch) — all pass.
- `R CMD check` / `devtools::test()` — accounting for the edits to
  `test_seed_consistency.R`, `test_top_level.R`, `inst/scripts/debugRuns.R`
  (§4.7-related).
- Check both build paths where things diverged before (see issues 9–10):
  `devtools::install_local()` (release, `-O2`), not only `load_all()`.

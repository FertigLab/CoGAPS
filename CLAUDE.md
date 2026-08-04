# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

CoGAPS (Coordinated Gene Activity in Pattern Sets) is a Bioconductor R package wrapping a C++
Bayesian MCMC matrix factorization algorithm (GAPS). It factors a data matrix `D ≈ A × P`, where
`A` (featureLoadings) is genes × patterns and `P` (sampleFactors) is patterns × samples, and links
the result to gene set statistics.

An earlier assistant-facing description lives in `.sasha-copilot-context/AGENTS.md`; working
agreements with the maintainer are in `.sasha-copilot-context/copilot-rules.md`. Per-issue notes for
in-flight work sit in `.sasha-copilot-context/<issue>/` — including the log of fixed defects
(`132-LLM-assisted-solved-issues.md`), which is the best entry point for why the C++ looks the way
it does. The directory is excluded from the package build via `.Rbuildignore`.

## Setting up a fresh machine

Verified against R 4.6.1 / Bioconductor 3.23.

```r
install.packages("BiocManager")
BiocManager::install(c("devtools", "testthat", "roxygen2", "BiocStyle",
                       "SingleCellExperiment", "fgsea", "gplots", "SeuratObject"))
devtools::install_deps(".", dependencies = TRUE)   # the rest, 21 in DESCRIPTION
```

Three of these are easy to overlook, because nothing fails until it does:

- **`testthat`** is in `LinkingTo`, not just `Suggests` — it ships the Catch2 header
  the C++ test suite compiles against. Without it `src/` does not build at all.
- **`xml2`** is what `tests/testthat/test_cpp.R` uses to turn the Catch report into
  per-case expectations. It is only suggested, so the test silently falls back to a
  single pass/fail check when it is absent.
- **`SeuratObject`** is only needed for `R CMD check`, which refuses to run a
  complete check while a suggested package is missing.

Outside R:

```bash
brew install autoconf autoconf-archive   # or the distro equivalent
```

`autoconf-archive` is required to regenerate `configure` — see the note under
Build & test. It is only needed when `configure.ac` changes.

Do **not** upgrade the generated documentation casually: `DESCRIPTION` pins
`RoxygenNote: 7.3.3`, and running `devtools::document()` under a newer roxygen2
rewrites all 77 `.Rd` files.

## Build & test

All commands are run from the package root.

```r
devtools::load_all()                  # load package
devtools::load_all(recompile = TRUE)  # force C++ recompile
devtools::test()                      # all R tests
testthat::test_file("tests/testthat/test_top_level.R")  # one R test file
devtools::document()                  # regenerate roxygen2 docs + NAMESPACE
Rcpp::compileAttributes()             # after changing // [[Rcpp::export]] signatures
```

```bash
R CMD check --no-manual .   # what CI runs (FertigLab/actions r-build-check)

# after editing configure.ac -- BOTH commands, in this order
aclocal -I /opt/homebrew/share/aclocal   # autoconf-archive macros
autoconf
```

`autoconf` alone is not enough: `configure.ac` uses `AX_COMPILER_VENDOR` and
`AX_COMPILER_VERSION` from autoconf-archive, and it is `aclocal` that makes them
available. Skipping it leaves both macros unexpanded in `configure`, where they
become literal shell commands — `configure` still completes, but prints
`AX_COMPILER_VENDOR: command not found`, leaves `$ax_cv_cxx_compiler_vendor`
empty, and so silently turns `--enable-warnings` into a no-op. That is the state
`configure` is in on `master`; on this branch it is generated correctly, so keep
it that way. `aclocal.m4` is an artefact of this and is not committed.

Note that `--enable-warnings` currently fails the build with `-Werror`, on
`-Wcast-function-type-mismatch` raised inside Rcpp's own `routines.h` under
newer clang. No CoGAPS source file produces a warning.

`DESCRIPTION` pins `RoxygenNote: 7.3.3`. Running `devtools::document()` under a newer roxygen2
rewrites all 77 `.Rd` files and bumps that field — check the resulting diff before committing it.

### C++ unit tests (Catch2, shipped with the `testthat` R package)

These are **not** exported — reach them with `:::`. The file-parser tests read their data paths
from the global environment, so set `gistCsvPath`/`gistTsvPath`/`gistMtxPath`/`gistGctPath` (with
`<<-`) before calling the runner directly.

```r
CoGAPS:::run_catch_unit_tests()                          # all, console output
CoGAPS:::run_catch_unit_tests_by_tag("Test Vector.h")    # by name
CoGAPS:::run_catch_unit_tests_by_tag("[vector]")         # by tag
CoGAPS:::run_catch_unit_tests_by_tag("[vector][green]")  # AND
CoGAPS:::run_catch_unit_tests_by_tag("[vector],[green]") # OR
CoGAPS:::catch_test_case_names()                         # what is compiled in
```

The suite also runs inside `devtools::test()` via `tests/testthat/test_cpp.R`, so a broken C++
test breaks `R CMD check`. That file does **not** use the console form: Catch writes its report
from C++ directly to stdout, where testthat cannot capture it, and one `expect_equal(…, 0L)` would
collapse the whole suite into a single pass/fail. It passes `reporter = "xml"` plus an `output`
file instead, then parses it with `xml2` so each `<TestCase>` becomes its own expectation, named
and located on failure. `output = ""` (the default) still writes to stdout, which is what keeps
the interactive call unchanged.

`catch_test_case_names()` guards against a vacuous pass: with `--disable-cpp-tests`, or on Windows
where `Makevars.win` lists no `cpp_tests` objects, no test case is registered and the runner
returns 0 regardless. The test skips in that case rather than reporting success. See
`src/cpp_tests/README.md` for the longer write-up.

### Build options

Options are declared in `configure.ac`; pass them by setting the matching env var before loading:

```r
Sys.setenv(enable_debug = "yes")      # -g -O0
devtools::load_all(recompile = TRUE)
```

Available toggles: `--enable-debug` (`-g -O0`), `--enable-gaps-debug` (`-DGAPS_DEBUG`),
`--enable-cpp-tests` (on by default; disable to cut compile time), `--enable-checkpoints`
(**off** by default — `-DGAPS_DISABLE_CHECKPOINTS` is set unless enabled), `--enable-warnings`
(`-Wall -Wextra -Werror`), `--enable-simd` (on by default; `sse` disables AVX).

`src/Makevars.win` is maintained by hand (Windows has no `configure`): it hard-codes
`-DGAPS_DISABLE_CHECKPOINTS` and lists no `cpp_tests/*.o`, so the Catch suite is empty on Windows
and `test_cpp.R` passes trivially there. Adding a source file means editing **both**
`configure.ac` and `Makevars.win`.

## Architecture

### Layers

```
R/                      Public API and S4 classes
  CoGAPS.R              CoGAPS(), scCoGAPS(), GWCoGAPS() entry points
  DistributedCogaps.R   subset orchestration, pattern matching, stitching
  SubsetData.R          how data is split into sets (explicit / weighted / uniform)
  class-*.R             S4 class definitions + generics
  methods-*.R           S4 method implementations
  HelperFunctions.R     input validation, dim names, file/RDS handling
  RcppExports.R         auto-generated — never edit by hand

src/                    C++ core
  GapsRunner.cpp/.h     C++ entry: gaps::run() (in-memory or from file)
  Cogaps.cpp            Rcpp glue (cogaps_cpp, cogaps_from_file_cpp)
  GapsParameters/Result/Statistics
  gibbs_sampler/        DenseNormalModel, SparseNormalModel, SingleThreadedGibbsSampler
  atomic/               atomic domain backing the Gibbs sampler
  data_structures/      Matrix, SparseMatrix, HybridMatrix, Vector, SparseIterator, ...
  file_parser/          CSV / TSV / MTX / GCT readers
  math/                 Random, VectorMath, MatrixMath, SIMD
  cpp_tests/            Catch2 tests; test-runner.cpp exposes them to R
```

### Dispatch

`CoGAPS(data, params, nPatterns, ...)` validates inputs, then picks one of three paths:

- `cogaps_cpp` — in-memory matrix (default)
- `cogaps_from_file_cpp` — `data` is a file path
- `distributedCogaps` — `params@distributed` is `"genome-wide"` or `"single-cell"`; splits the data
  into subsets, runs them in parallel via BiocParallel, matches patterns across sets
  (`findConsensusMatrix`/`patternMatch`), then `stitchTogether`s the result

`nPatterns` is required — the `CogapsParams` initializer errors without it, and the `params`
default is `new("CogapsParams", nPatterns = nPatterns)`. Distributed runs want on-disk data
(mtx/tsv/csv/gct); an in-memory matrix warns.

In `stitchTogether`, the **fixed** matrix must be read from
`result[[1]]@metadata$params@fixedPatterns`, not from the worker's `@featureLoadings` /
`@sampleFactors`: when a matrix is fixed its statistics are never accumulated, so those slots are
all zeros.

### Sampler geometry (see `src/README.md` for the full write-up)

Every run drives two samplers over `D ≈ A × P`. Each sampler owns its own copy of `AP` (synced by
`sync()`, built by `extraInitialization()`), its purpose matrix, an uncertainty matrix the size of
`AP`, and a `const` reference to its counterpart. One sampler of the pair is transposed:

- **ASampler** (transposed when `D` is not transposed): `AP` is l×m, `A` is m×k
- **PSampler** (non-transposed): `AP` is m×l, `P` is l×k
- Invariants: `nrows(MyMatrix) == ncols(APMatrix)`, `ncols(MyMatrix) == k`, and after `sync()`,
  `APMatrix == t(other_sampler.APMatrix)`

### Uncertainty model

`S = max(0.1·D, 0.1)` by default — a relative-error model floored so that zeros and `D < 1` get
`S = 0.1`. `DenseNormalModel` materialises it as `mSMatrix`; `SparseNormalModel` never materialises
it (that would defeat sparse storage) and instead recomputes `1/S²` on the fly via the file-local
`invSSq()` helper, with the constant zero-entry uncertainty folded into `mBeta`. The two must agree
— that is what makes `sparseOptimization=TRUE` match a dense run. A user-supplied `uncertainty=`
matrix is used verbatim, and `sparseOptimization` rejects one outright.

**`Vector::pad(val)` fills every allocated element; `padSIMD(val)` fills only the SIMD tail.** Never
use `pad()` on a matrix whose contents matter — that bug silently discarded user uncertainty for
years. SIMD loops read up to `SIMD_INC - 1` elements past the end, so the tail must be non-zero to
avoid `0/0 = NaN`.

### Threading

The asynchronous/OpenMP multi-threaded sampler was **removed** — it broke MCMC detailed balance.
CoGAPS always runs single-threaded. `nThreads` and `asynchronousUpdates` remain in the `CoGAPS()`
signature for backward compatibility, warn when set to a non-default value, and are otherwise
ignored. `compiledWithOpenMPSupport()` is kept, also for compatibility, and always returns `FALSE`.
Do not reintroduce any of these into `allParams` or the C++ parameter struct.

## Conventions

### R

- S4 throughout: new generics go in `class-*.R`, implementations in `methods-*.R`.
- `CogapsResult` extends `LinearEmbeddingMatrix`: `featureLoadings` = Amean, `sampleFactors` = Pmean,
  `loadingStdDev` = Asd, `factorStdDev` = Psd.
- Distributed params (`nSets`, `cut`, `minNS`, `maxNS`) cannot be passed to the `CogapsParams`
  constructor — the initializer errors out. Use `setDistributedParams()` afterwards.
- `scCoGAPS()` / `GWCoGAPS()` are deprecated wrappers for `CoGAPS(..., distributed = ...)`.
- New files must be added to the `Collate:` field in `DESCRIPTION`.

### C++

- Each Catch2 test file needs both `#include <testthat.h>` and `#include "testthat-tweak.h"`
  (the tweak enables tag-based selection). **A test case with no `SECTION` block is never run.**
- Adding any new `.cpp` (source or test) requires adding its `.o` to `GAPS_SOURCE_FILES` in
  `configure.ac`, then running `autoconf`. There is no wildcard build.

### Working with the maintainer

- Conversation is in Russian; everything an outside reader sees — comments, docs, README, commit
  messages — is in English.
- **Commit only when explicitly asked.** Never commit on your own initiative.
- Do not rewrite history (`rebase -i`, force push) once Bioconductor review has started.
- Don't do things "just in case" without asking; ask rather than guess when something is unclear.
- Comments only where a clarification is genuinely needed.

### Merging master into a long-lived branch

`git checkout --ours <file>` takes the branch's version of the **whole file**, silently discarding
master changes that had merged cleanly elsewhere in it. Use `git checkout -m -- <file>` to get the
real three-way merge back and resolve only the conflicting hunks.

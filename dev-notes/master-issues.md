# Defects inherited from master, fixed on this branch

Branch `132-uncertainty-improvements` is about the uncertainty model. It also
carries a handful of fixes that have nothing to do with uncertainty. This file
records why, so that a reviewer does not have to wonder where they came from.

All of them were found while merging `master` into the branch (2026-08) and while
writing tests for exported functions that had never been called from any test —
21 of the 37 exported functions had no coverage at all. None of them originate on
this branch: the branch point is `8a21f281` (merge of PR #131), and the defects
are either older than that or came in from `master` afterwards.

A longer, Russian-language write-up of items 1–5 was prepared separately for the
`master` maintainer; this file is the version that belongs with the code.

| # | What | Origin | Fixed in |
|---|------|--------|----------|
| 1 | `scCoGAPS()` / `GWCoGAPS()` unusable | `24448b6e` on master | `661a0828` |
| 2 | `patternMarkers` lp test asserted nothing | `24448b6e` on master | `4e74b1bf` |
| 3 | `binaryA()` unusable | older than `8a21f281` | `40eb0f58` |
| 4 | `fromCSV()` returned the wrong slot types | older than `8a21f281` | `40eb0f58` |
| 5 | `show(CogapsParams)` errored | older than `8a21f281` | `40eb0f58` |

---

## 1. `scCoGAPS()` / `GWCoGAPS()` failed on every call

`24448b6e` ("require nPatterns always") made `nPatterns` a mandatory argument of
the `CogapsParams` initializer, replacing `.Object@nPatterns <- 7`. `CoGAPS()`
was updated to match — `params = new("CogapsParams", nPatterns = nPatterns)` —
but the two deprecated wrappers kept `params = new("CogapsParams")`, which can no
longer be constructed.

The non-obvious part: passing `nPatterns` did not help either. It goes into `...`
and would reach `CoGAPS()`, but the first statement of each wrapper is
`params@distributed <- ...`, which forces the lazy default and fails first. So
the wrappers worked only with an explicitly built `params` object, and the error
message pointed nowhere useful.

Fixed by giving both wrappers the same signature treatment as `CoGAPS()`.
`tests/testthat/test_deprecated_wrappers.R` now calls them with and without an
explicit `params` — they had no test coverage and no runnable example, which is
why this went unnoticed.

## 2. The `patternMarkers` lp test stopped testing anything

The same commit raised `nPatterns` from 5 to 7 in `test_patternMarkers.R` but
left the `lp` vectors at length 5. `patternMarkers()` then warned "lp length must
equal the number of columns of the Amatrix" and took the invalid-lp path, so the
case meant to cover a well-formed `lp` never exercised it. `expect_no_error()`
does not catch a warning, so the test stayed green while the noise accounted for
2727 of the suite's 2730 warnings.

Fixed by padding the two vectors that are supposed to be valid to length 7. The
middle call keeps its length-4 vector: it asserts the length warning on purpose.

## 3. `binaryA()` failed on every call

```r
binA <- ifelse(calcZ(object) > threshold, 1, 0)
```

`calcZ()` has no default for `whichMatrix`, so every call died with
`argument "whichMatrix" is missing`. Fixed by passing `"featureLoadings"` — the
function is named binary**A** and draws a "Heatmap of Standardized Feature
Matrix".

Still open, left alone deliberately: the name promises a binary matrix, but the
function returns whatever `mtext()` returns. It is a plotting function; changing
its contract is a separate decision.

## 4. `fromCSV()` returned data.frames where matrices were expected

All four matrix slots were read with `read.csv()` and handed to
`new("CogapsResult", ...)`. `CogapsResult` extends `LinearEmbeddingMatrix`, whose
slots are matrices, so a `toCSV()`/`fromCSV()` round-trip produced an object with
the wrong slot types. Values and dimnames did survive; only the type was wrong.
Fixed by wrapping the reads in `as.matrix()`.

## 5. `show()` on `CogapsParams` errored with a checkpoint file set

```r
cat("checkpointInFile          ", checkpointInFile, "\n")
```

Missing `object@`, so R looked for a global variable and printing any params
object with `checkpointInFile` set raised `object 'checkpointInFile' not found`.
Fixed by qualifying it.

This one had been visible all along in `R CMD check` as the "no visible binding
for global variable ‘checkpointInFile’" NOTE — nobody had gone through the NOTEs.

---

## Known gaps left open

Not defects introduced here, and not fixed here either. Recorded so the next
person does not have to rediscover them.

- **`nPatterns` has no default any more.** Dropping it broke `new("CogapsParams")`
  and `scCoGAPS(data)`, both of which worked in 3.27. The old default of 7 was
  arbitrary, but its removal is a backwards-compatibility break worth a conscious
  decision.
- **`configure` on master leaves `AX_COMPILER_VENDOR` / `AX_COMPILER_VERSION`
  unexpanded.** They end up in `configure` as literal shell commands: it still
  runs, but prints `command not found`, leaves `$ax_cv_cxx_compiler_vendor` empty
  and thereby makes `--enable-warnings` a silent no-op. Regenerating needs
  `aclocal` before `autoconf`; see `src/README.md`. The branch's `configure` is
  correct.
- **`--enable-warnings` fails the build under `-Werror`** — on
  `-Wcast-function-type-mismatch` raised inside Rcpp's own `routines.h` with
  newer clang. No CoGAPS source file produces a warning.
- **Windows builds no C++ tests.** `src/Makevars.win` is maintained by hand and
  lists no `cpp_tests` objects, so the Catch suite is empty there.
  `tests/testthat/test_cpp.R` now skips with an explanation instead of passing
  vacuously.
- **`R CMD check` still reports one NOTE** — `std::cout` / `printf` in the
  compiled code (`GapsPrint.h` and Catch2). Clearing it means routing C++ output
  through `Rprintf` throughout.
- **12 exported functions still have no test**, down from 21: `MANOVA`,
  `buildReport`, `compiledWithOpenMPSupport`, `findConsensusMatrix`,
  `getClusteredPatterns`, `getCorrelationToMeanPattern`, `getParam`,
  `getRetinaSubset`, `getUnmatchedPatterns`, `plotPatternMarkers`,
  `plotResiduals`, `setAnnotationWeights`.
- **Weighted subset sampling is still broken** — separate report in
  `annotation-weights-sampling-issue-eng.md`.

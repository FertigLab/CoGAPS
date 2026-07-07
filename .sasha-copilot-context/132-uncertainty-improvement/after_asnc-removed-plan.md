# Post-async-removal: test-coverage audit & suspicious-code plan

**Branch:** `132-uncertainty-improvements`
**Date:** 2026-07-07
**Context:** after the async sampler removal (`async-removal.md`), a review of the
C++ test suite and the surviving code for latent bugs and coverage gaps.

---

## 1. Do all C++ tests compile? — Yes

All 15 test files in the build list (`configure.ac` `GAPS_SOURCE_FILES`) compile
and pass: **46 Catch test cases**. The only file outside the build is
`src/cpp_tests/static_cast_standalone_test.cpp` — a standalone reproducer with its
own `main` (intentional, from the `static_cast<uint64_t>` overflow investigation).

**Caveat:** 6 registered TEST_CASEs are **empty stubs** (`{}`) — they "pass" while
asserting nothing (see §3).

---

## 2. Suspicious code — real latent bugs (ranked)

| ID | Location | Defect | Severity | Trigger |
|----|----------|--------|----------|---------|
| **B1** | `gibbs_sampler/SingleThreadedGibbsSampler.h:265` | `operator>>` reads the `DataModel` with `>>` but then chains `<<` (write) for `mDomain`, `mNumBins`, `mBinLength`, `mNumPatterns`, `mdDomainLength`, `mAlpha`. Checkpoint **restore does not reload the atomic domain / bin geometry** — it writes them into the read archive instead. Compare the correct write path at line 257. | High\* | any checkpoint save→restore |
| **B2** | `math/MatrixMath.cpp:54` (dense), `:72` (sparse) | `nonZeroMean` returns `sum / nNonZeroes` with no guard for `nNonZeroes == 0` → `0/0 = NaN`. Feeds `DenseNormalModel.h:95` / `SparseNormalModel.h:82`: `meanD = nonZeroMean(mDMatrix); mLambda = alpha*sqrt(nPatterns()/meanD)` → poisons `mLambda`, `mMaxGibbsMass`, and every sample. | Med-High | all-zero data matrix or subset |
| **B3** | `math/VectorMath.cpp:7,17,40,50,74` | Dense/hybrid `min`/`max`/`whichMax(Vector/HybridVector)` dereference `v[0]` with no size check. Same empty-container class just fixed for the *sparse* overloads (issue 12); the dense/hybrid ones were left unguarded. `MatrixMath.h:37,49` `min/max(MatrixType)` also call `getCol(0)` with no `nCol()>0` check. | Med | `gaps::max(Vector(0))`, 0-column matrix |
| B4 | `atomic/AtomicDomain.cpp:57,60` | `randomFreePosition` uses `uniform64(1, mDomainLength)` inclusive of `mDomainLength`; a returned `pos == mDomainLength` gives `pos/mBinLength == nElements` → one-past-last row → OOB in `changeMatrix`/`updateAPMatrix`. `move()` (`:197`) correctly uses `rbound - 1`. | Med | `uniform64` returns exactly `mDomainLength` (~1/2^64) |
| B5 | `gibbs_sampler/SingleThreadedGibbsSampler.h:241-243` | `exchange()` computes `mass.value()` **before** the `hasValue()` guard. Benign today only because `OptionalFloat::value()` returns `0.f` when empty and `&&` short-circuits; wrong evaluation order. | Med-Low | `sampleExchange` returns empty `OptionalFloat` |
| B6 | `math/MatrixMath.cpp:19,35`; `MatrixMath.h:72` | `float size = mat.nRow() * mat.nCol();` multiplies two 32-bit `unsigned` in 32-bit arithmetic before widening → overflow for matrices > ~4.3e9 cells → wrong `sparsity`/`mean`. | Med | very large data (e.g. 50k × 100k) |

**\* B1 scope:** checkpoints are **disabled in the default build**
(`-DGAPS_DISABLE_CHECKPOINTS` in `configure.ac:48` and `Makevars.win`;
`checkpointsEnabled() == FALSE`). B1 is therefore **latent** — not reachable by
users of the shipped package — but it is a genuine correctness bug and a one-token
fix (`<<` → `>>`), and a serialization round-trip test would catch it.

Additional robustness concerns noted (lower priority): `exponential()`
(`math/Random.cpp:178`) can yield `+inf` when `uniform()` returns exactly 0;
`mMaxGibbsMass /= mLambda` divides by zero if `alpha == 0`; `mean()`/`sparsity()`
divide by `nRow*nCol` with no empty-matrix guard.

**Verified OK (looked suspicious, but safe):** `SparseIterator<1>` underflow of
`mSparseIndex` on empty is guarded by `atEnd()` before any deref; `AtomicDomain`
erase/move iterator fix-ups (the previously-fixed area) are self-consistent.

---

## 3. Test-coverage gaps

### 3.1 Empty / stub TEST_CASEs (registered, assert nothing) — all in `testSerialization.cpp`

| Line | Tag |
|------|-----|
| 114 | `HybridVector Serialization` |
| 118 | `SparseVector Serialization` |
| 165 | `HybridMatrix Serialization` |
| 169 | `SparseMatrix Serialization` |
| 248 | `GapsParameters Serialization` — checkpoint/restart relies on this; we just kept its fields for format compatibility |
| 252 | `GapsStatistics Serialization` |

Serialization is the checkpoint/resume mechanism; 6 of 12 serializable types have a
registered-but-empty test.

### 3.2 No live coverage at all
- `Random.h`: `truncNormal`, `truncGammaUpper`, `poisson`, `exponential` — the core
  proposal distributions of the Gibbs sampler. (Covered only under `#if 0` in
  `testRandom.cpp:127-204`.)
- `Math.h`: `d_gamma`/`p_gamma`/`q_gamma`/`d_norm`/`p_norm`/`q_norm` — only under
  `#if 0` (`testRandom.cpp:206-214`).
- `AlphaParameters` dense-vs-sparse consistency (the numerical heart of the sampler)
  — only under `#if 0` (`testSparseGibbsSampler.cpp:40-248`).
- `VectorMath`/`MatrixMath`: `whichMax`, `elementSq`, `pmax`, `dot_diff`,
  `sparsity`, `nonZeroMean`, `mean` — no dedicated tests.

### 3.3 Shallow tests
- `testSamplerHighLevel.cpp:53` "Sampler Update" — constructs but never calls
  `update()` (needs `sync()`+`extraInitialization()` first); no assertions.
- `testSparseGibbsSampler.cpp:8` — only checks initial `chiSq()`; never calls
  `update()`. The dense test asserts chiSq decreases; the sparse one does not.

### 3.4 Dead code
- `gibbs_sampler/SparseNegativeBinomialModel.{h,cpp}` — 0-byte files, referenced
  nowhere, not compiled. Delete or mark WIP.

### 3.5 Disabled tests worth reviving (`#if 0`)
- `testRandom.cpp:127-214` — poisson/exponential means + gamma/norm exact values
  (High; needs port from old `gaps::random::` global API to `GapsRng`).
- `testSparseGibbsSampler.cpp:40-248` — dense-vs-sparse `alphaParameters` equality
  (High; needs port to `SingleThreadedGibbsSampler<SparseNormalModel>` API).
- `testSerialization.cpp:256-295` — `AtomicDomain` round-trip (Med-High; needs a
  friend/accessor for `mAtoms`/`mDomainLength`).

---

## 4. Plan

**Phase 1 (now): fix confirmed bugs B1 + B2 + B3 with regression tests.**
- B1: `<<` → `>>` in the sampler deserialization operator; regression = a sampler
  serialize→deserialize round-trip. (Also fills the empty serialization stubs
  conceptually.)
- B2: guard `nonZeroMean` for `nNonZeroes == 0` (return 0); test on an all-zero
  matrix (dense + sparse).
- B3: guard dense/hybrid `min`/`max`/`whichMax` (and `min/max(Matrix)`) for empty
  input; test size-0 vector / 0-column matrix.
- Document each as an entry in `132-LLM-assisted-solved-issues.md` (fixed bugs).

**Phase 2 (later): coverage.**
- Fill the 6 empty serialization TEST_CASEs.
- Revive `testRandom` distribution tests and the `AlphaParameters` dense-vs-sparse
  consistency test against the current API.
- Make `testSamplerHighLevel`/`testSparseGibbsSampler` actually run `update()` and
  assert chiSq decreases.

**Phase 3 (cleanup): decide the fate of `SparseNegativeBinomialModel.{h,cpp}`.**

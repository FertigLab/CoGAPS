# LLM-Assisted Solved Issues — Branch `132-uncertainty-improvements`

Branch tracks GitHub issue #132 ("Uncertainty improvements").
The issues below were diagnosed and solved with LLM assistance.
Issue numbering continues from the manually-fixed set (`132-manually-fixed-issues.md`,
issues 1–8); cross-references such as "issue 9" refer to that shared numbering.

---

## 9. SIMD NaN in `alphaParameters` on Mac Intel (SSE4/AVX)

Full root-cause analysis is in `dev-notes/simd-issue.md`.

The bug manifested only on Mac Intel, where Apple Clang enables SSE4.1 (or AVX)
by default.  On this platform the SIMD loop in `alphaParameters()` increments
its index by 4 (or 8) per step, so the last iteration reads past the last real
data element into the alignment-padding positions that `Vector` allocates but
does not initialise to anything meaningful.

The uncertainty matrix `mSMatrix` is computed by `gaps::pmax(mDMatrix, factor,
mLambda)`, which before the fix filled only the real elements.  The padding
positions were left at zero — the default from the `Vector` constructor.  When
the SIMD loop read one of those padding positions, both the numerator (`pMat`)
and the denominator (`pS * pS`) were zero, giving `0/0 = NaN`.  That NaN
propagated through the horizontal sum of the SIMD register and poisoned the
returned `AlphaParameters`.  A poisoned `AlphaParameters` causes `gibbsMass()`
to return an empty result, so every `sampleBirth()` call was silently rejected.
Consequently PSampler never accumulated atoms: `MyMatrix` stayed all-zero, the
AP product stayed all-zero, and chi-square never decreased no matter how many
iterations were run.

The fix adds `padSIMD(float val)` to `Vector` and `Matrix`, which fills only the
padding positions (beyond `mSize`) with the given value.  Both `gaps::pmax`
overloads call `padSIMD(min_threshold)` before returning.  Padding lanes now
hold `mLambda > 0`, so the ratio there is `0 / mLambda^2 = 0` — no NaN, no
contribution to the accumulator.

Fix commit: `bfd5e09d`

---

## 10. `SIMD_PAD` insufficient on ARM — compiler auto-vectorization heap corruption

### Affected platform

Apple Silicon (M1/M2/M3/M4).  Does **not** affect Intel Mac or Linux.

### Root cause

On ARM there is no SSE/AVX, so `SIMD_INC = 1` and the `SIMD_PAD` macro
allocated only one extra element per `Vector`:

```c
SIMD_PAD(n) = 1 + 1*(n/1) = n + 1
```

The hand-written SIMD loops in `alphaParameters()` and `updateAPMatrix()` were
safe with `SIMD_INC = 1` — they iterate one float at a time, exactly `n` times.
However, when the package was built with optimisation (`-O2` / `-O3`, the R
default in `install_local`), Clang's auto-vectoriser replaced those scalar loops
with 4-wide ARM NEON instructions.  The auto-vectorised code processes floats in
batches of 4, so the last batch for a vector of 25 elements starts at index 24
and reads (or writes) indices 24, 25, 26, 27 — but only 26 elements were
allocated.  Indices 26 and 27 are past the end of the `std::vector` storage.

The writes in `updateAPMatrix` (`pAP.store(ap + i)`) were particularly
destructive: they silently overwrote adjacent heap memory, which could corrupt
the data of a neighbouring `Vector` — for example a column of `mSMatrix`.  When
`chiSq()` subsequently checked `GAPS_ASSERT(mSMatrix(i,j) > 0.f)`, the
corrupted value failed the assert and the sampler terminated with "CoGAPS
terminated".

The bug was invisible with `devtools::load_all()` because that build retains
the existing `src/Makevars` which had `-g -O0` (debug flags from a previous
configure run), suppressing all auto-vectorisation.  `devtools::install_local()`
runs `configure` from scratch, produces a release `Makevars` without `-O0`, and
the bug manifested immediately.

### Relationship to the SSE4/AVX fix (issue 9)

Issue 9 addressed padding values (zeros → `mLambda`).  This issue addresses
padding size: even with correct values in padding, there was simply not enough
padding for the auto-vectoriser to stay within bounds.  Both fixes are required
for correctness across all platforms.

### Fix

Changed `SIMD_PAD` in `Vector.cpp` to always use an effective width of 8,
regardless of `SIMD_INC`:

```c
#define SIMD_PAD_INC 8
#define SIMD_PAD(x) (SIMD_PAD_INC + SIMD_PAD_INC * ((x) / SIMD_PAD_INC))
```

`SIMD_PAD(25)` now allocates 32 elements on every platform, leaving 7 padding
positions (indices 25–31).  The auto-vectoriser on ARM reads/writes at most up
to index 27 — safely within bounds.  The `padSIMD(mLambda)` calls already
present in `gaps::pmax` fill all padding positions with a positive value, so
no NaN can arise there either.

Fix commit: `60f94063`

---

## 11. `AtomicDomain::move()` leaves a stale map iterator in the atom

### Problem

`AtomicDomain::move()` replaced the map entry for an atom (erase old key, insert new key) but never updated `atom->mIterator` to point to the new map node:

```cpp
std::pair<uint64_t, size_t> newpair(newPos, atom->iterator()->second);
mAtomMap.erase(atom->pos());   // atom->mIterator is now dangling
atom->updatePos(newPos);
mAtomMap.insert(newpair);      // atom->mIterator NOT updated — still dangling
```

After `move()`, `atom->mIterator` pointed to a deleted `std::map` node.  If the
same atom was subsequently passed to `erase()`, the line

```cpp
mAtomMap.erase(atom->iterator());   // UB: dangling iterator
```

caused undefined behaviour.  In practice this corrupted the red-black tree,
which meant the next `insert()` placed the atom at a wrong position in the
ordering.  Corrupted neighbour links eventually caused

```
GAPS_ASSERT(a <= b)   // lbound+1 > rbound-1
```

to fire inside `uniform64()` in `SingleThreadedGibbsSampler::move()`, terminating
with "CoGAPS terminated".

The bug triggered reliably with GIST data (1363 × 9, 9 patterns) but not with
the smaller random-matrix test, because only the GIST run produced a `move()`
followed immediately by an `erase()` on the same atom within the first few
thousand iterations.

### Fix

Capture the iterator returned by `mAtomMap.insert()` and store it in the atom:

```cpp
size_t storageIdx = atom->iterator()->second;
mAtomMap.erase(atom->pos());
atom->updatePos(newPos);
atom->setIterator(
    mAtomMap.insert(std::pair<uint64_t, size_t>(newPos, storageIdx)).first);
```

`AtomicDomain` is a friend of `Atom`, so it may call the private `setIterator()`.

A regression test `[atomicdomain][movethenerase]` was added to
`src/cpp_tests/testAtomicDomain.cpp`: it inserts four atoms, moves one, erases
the same atom, and verifies the map contains exactly the expected three keys.

### Changed files

| File | Change |
|---|---|
| `src/atomic/AtomicDomain.cpp` | `move()` calls `atom->setIterator()` after insert |
| `src/cpp_tests/testAtomicDomain.cpp` | regression test `[atomicdomain][movethenerase]` |

---

## 12. `gaps::min/max(SparseVector)` segfault on an empty (all-zero) sparse vector

### Problem

`gaps::min(const SparseVector&)` and `gaps::max(const SparseVector&)` read the
first element through the iterator **before** checking `atEnd()`:

```cpp
float gaps::max(const SparseVector &v)
{
    SparseIterator<1> it(v);
    float mx = get<1>(it);   // reads mData[0] with no bounds check
    while (!it.atEnd()) { ... }
    return mx;
}
```

For an **empty** sparse vector (one with no stored non-zero elements — i.e. an
all-zero row or column), the iterator is immediately `atEnd()` and `get<1>(it)`
dereferences `mData[0]` of an empty backing array → read at address `0x0`
(`SparseVector::getIthElement`), segfault.

This is reachable in normal use, not just tests: `SparseNormalModel`'s constructor
runs `if (gaps::max(mDMatrix) > 50.f)` (`SparseNormalModel.h:86`) on **every**
`sparseOptimization=TRUE` run. Any all-zero gene (row) or sample (column) in the
data — common in sparse single-cell matrices — yields an empty `SparseVector` and
crashes construction before sampling starts.

Related to issue 4 (`min`/`max` initialisation), whose fix did not cover the
empty-`SparseVector` case.

### Diagnosis

The crash surfaced only after `testSamplerHighLevel` was added to the build (it
constructs a `SparseNormalModel` on dummy data whose row 0 / column 0 are all
zero). Root-caused with a standalone AddressSanitizer harness compiled at `-O2`:
the backtrace pointed at `SparseVector::getIthElement` ← `gaps::max(SparseVector)`
← `SparseNormalModel` constructor.

### Fix

Guard the empty case — a sparse vector's absent elements are `0`, and CoGAPS data
is non-negative, so the min/max of an all-zero vector is `0.f`:

```cpp
SparseIterator<1> it(v);
if (it.atEnd()) return 0.f; // empty (all-zero) sparse vector
float mx = get<1>(it);
```

Applied to both `gaps::min` and `gaps::max`. Non-empty behaviour is unchanged.

A regression test `[sparsevector][emptyminmax]` was added to
`src/cpp_tests/testSparseVector.cpp`: it checks `min`/`max` return `0.f` on an
empty `SparseVector` (no crash) and are unaffected for a non-empty one.

### Changed files

| File | Change |
|---|---|
| `src/math/VectorMath.cpp` | `atEnd()` guard in `gaps::min`/`gaps::max` for `SparseVector` |
| `src/cpp_tests/testSparseVector.cpp` | regression test `[sparsevector][emptyminmax]` |

---

## 13. `SingleThreadedGibbsSampler` deserialization used `<<` (write) instead of `>>` (read)

### Problem

The sampler's `operator>>` (`src/gibbs_sampler/SingleThreadedGibbsSampler.h`) read
the `DataModel` base with `>>` but then **chained `<<` (write)** for the sampler's
own members:

```cpp
Archive& operator>>(Archive &ar, SingleThreadedGibbsSampler<DataModel> &s)
{
    operator>>(ar, static_cast<DataModel&>(s)) << s.mDomain << s.mNumBins
        << s.mBinLength << s.mNumPatterns << s.mdDomainLength << s.mAlpha;   // BUG: <<
    return ar;
}
```

So restoring a sampler from a checkpoint did **not** reload the atomic domain or the
bin geometry — those fields were written back into the read archive instead. A run
resumed from a checkpoint started with an empty atomic domain and diverged from an
uninterrupted run.

Found while auditing the code after the async removal. It is a latent bug in the
shipped package because checkpoints are disabled by default
(`-DGAPS_DISABLE_CHECKPOINTS`, `checkpointsEnabled() == FALSE`), but it is a genuine
correctness bug and would bite anyone building with checkpoints enabled.

### Fix

Use `>>` for the whole chain, matching the write path:

```cpp
operator>>(ar, static_cast<DataModel&>(s)) >> s.mDomain >> s.mNumBins
    >> s.mBinLength >> s.mNumPatterns >> s.mdDomainLength >> s.mAlpha;
```

Regression test `[serialization][gibbssampler-roundtrip]` in
`testSerialization.cpp`: build a dense sampler, `sync` + `extraInitialization` +
`update(500)` to populate the domain, serialize, deserialize into a fresh sampler,
and require the restored `nAtoms()` (and `MyMatrix` sum) to match the original.
Before the fix the restored domain stayed empty (`nAtoms() == 0`).

### Changed files

| File | Change |
|---|---|
| `src/gibbs_sampler/SingleThreadedGibbsSampler.h` | `operator>>` chain `<<` → `>>` |
| `src/cpp_tests/testSerialization.cpp` | regression test `[serialization][gibbssampler-roundtrip]` |
---

## 14. `gaps::nonZeroMean` divided by zero (NaN) on an all-zero matrix

### Problem

`gaps::nonZeroMean(const Matrix&)` and `(const SparseMatrix&)`
(`src/math/MatrixMath.cpp`) returned `sum / nNonZeroes` with no guard for
`nNonZeroes == 0`. For an all-zero data matrix (or a subset that contains no
positive entries) this is `0 / 0 = NaN`.

It feeds model initialization directly: `DenseNormalModel.h` / `SparseNormalModel.h`
compute `meanD = nonZeroMean(mDMatrix)` then
`mLambda = alpha * sqrt(nPatterns() / meanD)` and `mMaxGibbsMass /= mLambda`. A NaN
`meanD` poisons `mLambda`, `mMaxGibbsMass`, and every subsequent sample.

### Fix

Guard the empty case (a matrix with no positive entries has non-zero mean `0`):

```cpp
if (nNonZeroes == 0) return 0.f; // all-zero matrix: avoid 0/0 = NaN
return sum / static_cast<float>(nNonZeroes);
```

Applied to both the dense and sparse overloads. Regression test
`[matrix][nonzeromean-empty]` in `testMatrix.cpp` checks `nonZeroMean` of an all-zero
`Matrix` returns `0.f` (a NaN would fail the equality).

### Changed files

| File | Change |
|---|---|
| `src/math/MatrixMath.cpp` | `nNonZeroes == 0` guard in both `nonZeroMean` overloads |
| `src/cpp_tests/testMatrix.cpp` | regression test `[matrix][nonzeromean-empty]` |

---

## 15. Dense/hybrid `gaps::min`/`max`/`whichMax` read element 0 on an empty container

### Problem

The dense and hybrid reductions `gaps::min`/`max(Vector)`,
`gaps::min`/`max(HybridVector)`, `gaps::whichMax(Vector)`
(`src/math/VectorMath.cpp`) initialised the accumulator with `v[0]` with no size
check, and the templated `gaps::min`/`max(MatrixType)` (`src/math/MatrixMath.h`)
called `getCol(0)` with no `nCol() > 0` check. On an empty vector / zero-column
matrix this is an out-of-bounds read.

This is the same empty-container class fixed for the `SparseVector` overloads in
issue 12; the dense/hybrid/matrix siblings were left unguarded (and issue 4, which
changed the accumulator init from `0` to `v[0]`, is what introduced the
empty-vector deref). `v[0]` currently reads SIMD padding rather than crashing on
this platform, but it is a genuine OOB on a vector with no padding.

### Fix

Return `0.f` (index `0` for `whichMax`) for empty input, before any indexing:

```cpp
if (v.size() == 0) return 0.f;      // Vector / HybridVector min/max
if (v.size() == 0) return 0;        // whichMax
if (mat.nCol() == 0) return 0.f;    // min/max(MatrixType)
```

Regression tests: `[vector][emptyminmax]` in `testVector.cpp` (empty `Vector`
min/max/whichMax → 0) and the "zero-column matrix" section of `[matrix][minmax-empty]`
in `testMatrix.cpp`.

### Changed files

| File | Change |
|---|---|
| `src/math/VectorMath.cpp` | empty-size guard in dense/hybrid `min`/`max`/`whichMax` |
| `src/math/MatrixMath.h` | `nCol() == 0` guard in `min`/`max(MatrixType)` |
| `src/cpp_tests/testVector.cpp` | regression test `[vector][emptyminmax]` |
| `src/cpp_tests/testMatrix.cpp` | regression test `[matrix][minmax-empty]` (zero-column section) |

---

## 16. `SparseVector` deserialization dropped all stored values

### Problem

`operator>>(Archive&, SparseVector&)` read the stored (non-zero) values with a loop
bounded by the **destination's** current size:

```cpp
for (unsigned i = 0; i < vec.mData.size(); ++i)   // BUG: destination's count
    ar >> vec.mData[i];
```

`SparseVector` stores only its non-zero entries in `mData`; the count of non-zeros
is not written to the archive. When deserializing into a freshly-constructed
`SparseVector(size)` (empty `mData`), the loop ran zero times, so **none** of the
serialized values were read back — the restored vector was all-zero (and the
archive read pointer was left misaligned for whatever followed).

`SparseMatrix` serialization delegates to this operator per column, so it inherited
the same defect. Both are unused in the current checkpoint path (`SparseNormalModel`
serializes its `HybridMatrix` `mMatrix`, not the `SparseMatrix` `mDMatrix`), which
is why it was never noticed — the same "untested serialization operator" class as
issue 13.

Found while filling the empty `SparseVector`/`SparseMatrix` serialization test
stubs: the round-trip into an empty destination failed. (The pre-existing
`SparseMatrix` stub would have passed even unfixed, because a test that rebuilds the
read target from the same data gives it the right `mData` sizes by accident.)

### Fix

The number of stored values equals the popcount of the bit flags, which are read
first. Derive it and `resize` `mData` before reading (no archive-format change):

```cpp
unsigned nNonZeroes = 0;
for (unsigned i = 0; i < vec.mIndexBitFlags.size(); ++i)
    nNonZeroes += __builtin_popcountll(vec.mIndexBitFlags[i]);
vec.mData.resize(nNonZeroes);
for (unsigned i = 0; i < nNonZeroes; ++i)
    ar >> vec.mData[i];
```

Regression tests `[serialization][sparsevector]` and `[serialization][sparsematrix]`
in `testSerialization.cpp` round-trip into an empty / all-zero destination and
compare the dense reconstruction element-wise.

### Changed files

| File | Change |
|---|---|
| `src/data_structures/SparseVector.cpp` | derive non-zero count from bit flags, `resize` + read `mData` |
| `src/cpp_tests/testSerialization.cpp` | filled `[serialization][sparsevector]` and `[serialization][sparsematrix]` (read into empty target) |

---

## 17. `sparseOptimization` inconsistent with the dense sampler (uncertainty model)

### Problem

`sparseOptimization=TRUE` produced results inconsistent with the dense sampler:
the two paths implement the **same** statistical model, so the sparse version is
meant to be an optimization giving the same answer — but they had diverged.

Both compute the same alpha statistics for the Gibbs update,
`s = Σⱼ P²/S²`, `s_mu = Σⱼ P(D−AP)/S²`, differing only in the uncertainty `S`:

- **Dense** stores `mSMatrix` explicitly and divides by `S²`.
- **Sparse** never stores `S`; it bakes the assumption into the algebra
  (`mZ1` base = the `S=1` sum, per-non-zero corrections replace `S=1` with `S=d`),
  scaled by `mBeta = 100`. Since `1/factor² = 1/0.1² = 100 = mBeta`, the sparse
  model's *effective* uncertainty is `S = factor·D` for observed data and
  `factor` for zeros.

Two independent defects made the effective `S` differ:

1. **Dense zero-floor regression (this branch, from issue 2).** Issue 2 changed
   the dense uncertainty from `pmax(mDMatrix, 0.1f)` (i.e. `max(0.1·D, 0.1)`, floor
   `0.1`) to `pmax(mDMatrix, factor, mLambda)` (floor `mLambda ≈ 0.006`), intending
   `mLambda` as an "atom-size" floor. But `mLambda` is the atom-mass scale
   (used for `mMaxGibbsMass`), not an uncertainty. For zero data this made dense
   `S = mLambda ≈ 0.006` (weight `1/S² ≈ 27000`) versus the sparse `S = 0.1`
   (weight `100`) — a ~230× divergence on every zero entry.
   (The issue 2 note also misdiagnosed master as "S = constant 1 via `pad(1.f)`";
   in fact master's `Vector::pad` loops from `mSize`, so it only touches padding —
   master's real `S` was already `max(0.1·D, 0.1)`.)

2. **Missing floor in sparse (pre-existing, from master).** Sparse used `S = d`
   for *every* non-zero entry with no floor, so `S = factor·d`. For fractional
   data `0 < d < 1` this gives `S < factor`, i.e. a tiny measurement is treated as
   more precise than a zero (weight `1/(factor·d)²` blows up). Master's dense
   floored at `0.1` but master's sparse did not, so dense and sparse already
   mismatched on fractional data even on master.

Symptom on GIST (continuous data): dense-vs-sparse `A·P` reconstruction was
uncorrelated; the alpha statistics diverged by ~230× (caught by the revived
`[sparsegibbs]` consistency test, issue-independent).

### Fix

Use one uncertainty model in both — relative error floored at `factor`:
```
S[i,j] = max(factor · D[i,j], factor)     factor = 0.1
```
so zeros and any `D < 1` get `S = 0.1`, and `D ≥ 1` get `S = 0.1·D`.

- **Dense** (`DenseNormalModel.h`): revert to `gaps::pmax(mDMatrix, factor, factor)`
  (floor = `factor`, not `mLambda`).
- **Sparse** (`SparseNormalModel.cpp`): floor the raw uncertainty at 1
  (`S = max(d, 1)`, effective `factor·max(d,1)`) in all three `alphaParameters`
  functions. `d` plays a dual role (it is both the data value `D` and, previously,
  the uncertainty `S`); the floor lowers only the uncertainty while the data `d`
  is kept in the residual term (`s_mu += v·d·invS2 + v(1−invS2)·AP`,
  `invS2 = 1/max(d,1)²`).

### Verification

- `[sparsegibbs]` was revived (ported from the removed `GibbsSampler<Storage>` API
  to `SingleThreadedGibbsSampler<...NormalModel>` via a test-only `ExposedSampler`
  subclass that surfaces the protected `alphaParameters`) and its data extended to
  **continuous values including `0 < D < 1`**. It asserts sparse and dense
  `alphaParameters` (1D, 2D, symmetry, with-change) match to `TEST_APPROX`; now
  passes (8404 assertions).
- End-to-end on GIST, permutation-invariant metrics confirm equivalence: `A·P`
  reconstruction correlation dense-vs-sparse = `0.999`, equal to the
  dense-vs-dense (different-seed) control; `meanChiSq` for sparse lands within the
  dense-run range. (Raw `featureLoadings` correlation is meaningless here — NMF is
  only identifiable up to pattern permutation, so even two dense runs correlate at
  ≈ `−0.1`.)

Note: this deliberately changes dense/sparse numerics (the pre-fix behaviour was
wrong), so the async-removal parity baseline no longer applies to these paths.

### Changed files

| File | Change |
|---|---|
| `src/gibbs_sampler/DenseNormalModel.h` | uncertainty floor `mLambda` → `factor` (`pmax(mDMatrix, factor, factor)`) |
| `src/gibbs_sampler/SparseNormalModel.cpp` | floor `S = max(d, 1)` in the three `alphaParameters` functions |
| `src/cpp_tests/testSparseGibbsSampler.cpp` | revived `[sparsegibbs]` dense-vs-sparse consistency test; continuous (fractional) data |

## 18. `SparseNormalModel::chiSq()` segfaults when called before `sync()`

### Problem

`chiSq()` computes the fit `Σ (D − A·P)² / S²`, so it needs the *other* factor
matrix (`A` needs `P`, and vice versa). Each sampler receives that pointer,
`mOtherMatrix`, only when `sync()` is called; the constructor leaves it `NULL`.

`SparseNormalModel::chiSq()` dereferences `mOtherMatrix` directly
(`mOtherMatrix->getRow(i)` inside the dot products), so calling `chiSq()` on a
freshly-constructed, not-yet-`sync()`ed sampler is a NULL dereference — an
immediate segfault.

`DenseNormalModel::chiSq()` is **not** affected: it reads its cached, zero-
initialised `mAPMatrix` instead of `mOtherMatrix`, so before `sync()` it simply
returns the "no fit" value (`A·P = 0`). The two models therefore disagreed on
whether `chiSq()`-before-`sync()` is legal — dense tolerated it, sparse crashed.

Not reachable from the production run loop (which always `sync()`s during
initialization), but a real robustness defect for anyone using these classes as a
library, and it surfaced while writing the sampler unit tests.

### Diagnosis

AddressSanitizer backtrace: crash in `HybridMatrix::getRow` ←
`SparseNormalModel::chiSq()`, on the first `chiSq()` call, before any `sync()`.
The cause is the un-set `mOtherMatrix` (NULL) being dereferenced.

### Fix

Guard `SparseNormalModel::chiSq()` against `mOtherMatrix == NULL` and return the
same "no fit" value dense returns. With `A·P = 0`, every dot product is zero: the
`Σ (A·P)²` loop contributes nothing and each stored data value contributes its
`(D − 0)² / D² = 1` term, so the result is `(#non-zero entries) · mBeta`.

For data with all `D ≥ 1` this equals dense's `Σ D²/S² = Σ (D / 0.1D)² = 100`
per entry, i.e. the two models return the identical no-fit chiSq (verified in the
regression test on `D = i + j + 1`, giving `100 · nRow · nCol`).

### Changed files

| File | Change |
|---|---|
| `src/gibbs_sampler/SparseNormalModel.cpp` | `chiSq()` returns the no-fit value when `mOtherMatrix == NULL` instead of dereferencing it |
| `src/cpp_tests/testSparseGibbsSampler.cpp` | regression: `chiSq()` before `sync()` does not crash and equals dense's no-fit chiSq |

## 19. `SparseNormalModel::chiSq()` did not floor the uncertainty (inconsistent with dense and with issue #17)

### Problem

Issue #17 unified the uncertainty model on `S = max(factor·D, factor)`,
`factor = 0.1`, and floored it in all three `alphaParameters` functions — but
**not in `chiSq()`**, which was left computing the residual with the unfloored
`S = D` (`dsq = get<1>(it)²`). So for continuous data with `0 < D < 1` the sparse
`chiSq` used `S = D < 0.1` where the dense `chiSq` (and the sparse
`alphaParameters`) used the floor `S = 0.1`. Result: `meanChiSq` reported for a
`sparseOptimization=TRUE` run diverged from the equivalent dense run on fractional
data. Diagnostic only (goodness-of-fit / annealing readout — the *sampling*
decisions use the already-floored `alphaParameters`), but a real correctness gap
and a landmine if `chiSq` were ever used more widely.

### Root cause / design note

The sparse model deliberately does **not** materialise an `S` matrix (that would
defeat its whole purpose — avoiding dense `nRow×nCol` storage for single-cell
data — and would break the `mZ1`/`mZ2` lookup-table algebra, which depends on the
zero-entry uncertainty being a single constant that factors out into
`mBeta = 1/factor² = 100`). Instead every calculation computes `1/S²` on the fly
from the data value. The floor therefore has to be applied identically at each
call site; issue #17 did three of the four, and `chiSq()` was the straggler.

### Fix

Extract the one uncertainty model into a single file-local helper and call it
from all four sites:
```cpp
static inline float invSSq(float d)      // = 1/max(D,1)^2  (the factor is in mBeta)
{
    float sraw = gaps::max(d, 1.f);
    return 1.f / (sraw * sraw);
}
```
- The three `alphaParameters` functions now call `invSSq(d_val)` instead of the
  inlined `1/max(d,1)²` (no numeric change — just de-duplication).
- `chiSq()`'s non-zero correction becomes the floored residual: the first loop
  adds `A·P²` at the zero weight, and each stored entry corrects it to
  `(D − A·P)² · invSSq(D)` via `D²·invS2 − 2·D·A·P·invS2 + A·P²·(invS2 − 1)`
  (algebraically identical to the old `1 + dot(dot − 2d − dsq·dot)/dsq` when
  `D ≥ 1`, so integer/count data is unchanged).
- The `mOtherMatrix == NULL` no-fit branch (issue #18) likewise became
  `D² · invSSq(D)` so it stays consistent with the floored main path for `D < 1`.

`mZ1`/`mZ2` need no change: they encode the all-entries-are-zero baseline at the
constant `S = factor`, which the per-non-zero corrections adjust — already
consistent with the floor.

### Verification

Extended the `[sparsegibbs]` consistency test (continuous data with `0 < D < 1`,
identical `A`/`P` on a sparse and a dense sampler) to also assert
`sparse.chiSq() == Approx(dense.chiSq())` for both the A- and P-samplers. Fails
before the fix (sparse used `S = D` on the sub-1 entries), passes after. Full cpp
suite green.

### Changed files

| File | Change |
|---|---|
| `src/gibbs_sampler/SparseNormalModel.cpp` | new `invSSq()` helper; `chiSq()` (both branches) and all three `alphaParameters` now floor via it |
| `src/cpp_tests/testSparseGibbsSampler.cpp` | consistency test also compares sparse vs dense `chiSq()` on fractional data |

---

## 20. User-supplied `uncertainty=` was silently discarded (`pad()` overwrote the whole matrix)

### Problem

The `uncertainty=` argument of `CoGAPS()` had no effect whatsoever on a dense run.
Two runs differing only in the uncertainty matrix returned bit-identical results:

```r
r1 <- CoGAPS(D, nPatterns=3, nIterations=200, uncertainty=0.1*as.matrix(D), seed=1)
r2 <- CoGAPS(D, nPatterns=3, nIterations=200, uncertainty=10 *as.matrix(D), seed=1)
getMeanChiSq(r1) == getMeanChiSq(r2)   # TRUE -- 704.3173 in both cases
```

A hundredfold change in the uncertainty left the reported chi-square untouched,
i.e. the sampler was fitting with `S = 1` everywhere regardless of what the caller
passed.

### Root cause

`DenseNormalModel::setUncertainty()` loaded the caller's matrix and then called
`pad()`:

```cpp
mSMatrix = Matrix(unc, transpose, subsetRows, params.dataIndicesSubset);
mSMatrix.pad(1.f); // so that SIMD operations don't divide by zero
```

`Vector::pad(val)` fills **every** allocated element, not just the SIMD padding
lanes (`for (i = 0; i < mData.size(); ++i)`), so the line above replaced the whole
uncertainty matrix with 1.0f. The comment shows the intent was the padding lanes
only; the correct method is `padSIMD()`, which fills `[mSize, mData.size())` and
was added in this branch for the SIMD-NaN fix (issue #9).

The bug predates the branch — it is present in `master` too, on both the default
and the user-supplied path. Issue #17 removed the `pad(1.f)` from the *default*
path (replacing it with `gaps::pmax(mDMatrix, factor, factor)`) and so fixed that
half without the user-supplied half being noticed.

### Fix

```cpp
mSMatrix = Matrix(unc, transpose, subsetRows, params.dataIndicesSubset);
// Only the SIMD padding may be overwritten -- pad() would set *every* element
// to 1.f and so discard the uncertainty the caller passed in. Padding lanes
// get 1.f so that the SIMD loops divide by 1, not by 0.
mSMatrix.padSIMD(1.f);
```

`padSIMD(1.f)` keeps the SIMD-safety property (padding lanes divide by 1, never by
0 -- see issue #9) while leaving the caller's values intact.

### Verification

The same two runs now differ as they must (13338.07 vs 109.872). The chi-square
consistency test gained a case that recomputes the reported `getMeanChiSq()`
against an explicitly passed uncertainty matrix; it fails before the fix
(442 vs 103450) and passes after. That test came from `master`
(`test_chisq.R`), where it was written against the pre-#17 code, and is one of the
things the master merge brought in.

### Changed files

| File | Change |
|---|---|
| `src/gibbs_sampler/DenseNormalModel.h` | `setUncertainty()` uses `padSIMD(1.f)` instead of `pad(1.f)` |
| `tests/testthat/test_chisq.R` | added the explicit-uncertainty round-trip case |

### Note

`SparseNormalModel` is unaffected: `sparseOptimization=TRUE` rejects a
user-supplied uncertainty in `checkInputs()`, so it only ever uses the built-in
model.

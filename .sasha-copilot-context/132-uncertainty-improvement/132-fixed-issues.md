# Manually Fixed Issues — Branch `132-uncertainty-improvements`

Branch tracks GitHub issue #132 ("Uncertainty improvements").
All fixes below were made manually by @favorov.

---

## 1. `pmax` signature: separate `factor` and `min_threshold`

### Problem

The original `gaps::pmax(v, p)` prototype took a single float used simultaneously as
a multiplicative scaling factor **and** as the minimum-value threshold.  The two roles
are conceptually independent and the caller had no way to set them separately.

### Fix

Added a three-argument overload:

```cpp
Vector pmax(const Vector &v, float factor, float min_threshold);
Matrix pmax(const Matrix &m, float factor, float min_threshold);
```

The two-argument versions are kept for back-compatibility (they default
`min_threshold = factor`).  Parameters changed from pass-by-value to
`const &` and the result is now returned as a newly created object
(see §2 below).

### Changed files

| File | Change |
|---|---|
| `src/math/VectorMath.h` | new three-arg overload; `const Vector &` params |
| `src/math/VectorMath.cpp` | implementation of new overload |
| `src/math/MatrixMath.h` | same for Matrix |
| `src/math/MatrixMath.cpp` | same for Matrix |

Commits: `daf2cd94`, `eb767c3b`, `d8590ef`, `bd11860`

---

## 2. Data-driven uncertainty in `DenseNormalModel`

### Problem

The uncertainty matrix `mSMatrix` was initialised with `mSMatrix.pad(1.f)` — a
single constant for every element regardless of the data.  This makes the
chi-square denominator unresponsive to actual data magnitude.

### Fix

Replaced the constant pad with a data-driven computation:

```cpp
float factor = 0.1f;
mSMatrix = gaps::pmax(mDMatrix, factor, mLambda);
// mSMatrix[i,j] = max(mDMatrix[i,j] * factor, mLambda)
```

`mLambda` is set before `mSMatrix` is computed and represents the expected atom
size, so the minimum uncertainty is now data-scale-aware.

### Changed files

| File | Change |
|---|---|
| `src/gibbs_sampler/DenseNormalModel.h` | replace `pad(1.f)` with `gaps::pmax(mDMatrix, factor, mLambda)` |

Commit: `3efffbc1`

---

## 3. `elementSq` / `pmax` mutation bug

### Problem

`Vector elementSq(Vector v)` and `Vector pmax(Vector v, float p)` accepted vectors
by value, mutated the local copy, and returned it.  The intent was to return a new
object, but the pass-by-value idiom is fragile and was already causing confusion
with callers that expected the original to be unchanged.

### Fix

Changed signatures to `const Vector &` and constructed the result explicitly
before returning.

Commits: `d8590ef`, `bd11860`

---

## 4. `gaps::min` / `gaps::max` initialised with `0` instead of first element

### Problem

```cpp
float mn = 0.f;   // wrong: should be v[0]
for (...) mn = (v[i] < mn) ? v[i] : mn;
```

For any vector whose every element is positive, `gaps::min()` returned `0` instead
of the true minimum.  Symmetrically, `gaps::max()` returned `0` for all-negative
vectors.  This corrupted any downstream computation that depends on the true
data range.

Affected overloads: `min(Vector)`, `min(HybridVector)`, `min(SparseVector)`,
`max(Vector)`, `max(HybridVector)`.

### Fix

Initialise accumulators with the first element (`v[0]` / `get<1>(it)`) and start
the loop at index 1.

### Changed files

| File | Change |
|---|---|
| `src/math/VectorMath.cpp` | all five overloads corrected |
| `src/math/MatrixMath.h` | `min`/`max` for Matrix corrected |

Commits: `db47de2f`, `5741cfc8`

---

## 5. `Vector::pad` started at `mSize` instead of `0`

### Problem

```cpp
void Vector::pad(float val)
{
    for (unsigned i = mSize; i < mData.size(); ++i)  // wrong start
        mData[i] = val;
}
```

`pad` was supposed to overwrite **all** allocated elements with `val`.  Starting at
`mSize` left the real elements (`[0, mSize)`) unchanged, so calls like
`mSMatrix.pad(1.f)` silently left the matrix full of zeros.

### Fix

Changed loop start to `0`.

### Changed files

| File | Change |
|---|---|
| `src/data_structures/Vector.cpp` | loop start `mSize → 0` |

Commit: `88fc63e9`

---

## 6. `AtomicDomain`: removed `MutableMap`; fixed `erase` and `move`

### Problem

`MutableMap` was used to update the position key of an atom in-place via
`updateKey()`.  This operation is inherently unsafe (it bypasses std::map's
ordering invariant).  Additionally, `erase()` had two bugs:

1. When the last atom in `mAtoms` was swapped into the erased slot, the neighbour
   atoms' stored indices were not updated, leaving dangling left/right index references.
2. `dst.mIndex = index` assignment was redundant (the assert already verified
   equality); however the map entry was not being updated to reflect the
   post-swap storage index.

### Fix

- Replaced `mAtomMap.updateKey(...)` in `move()` with an explicit erase + insert:
  ```cpp
  std::pair<uint64_t, size_t> newpair(newPos, atom->iterator()->second);
  mAtomMap.erase(atom->pos());
  atom->updatePos(newPos);
  mAtomMap.insert(newpair);
  ```
- After copying the last atom into the erased slot, updated the left/right
  neighbour references:
  ```cpp
  if (dst.hasLeft())  mAtoms[dst.leftIndex()].setRightIndex(index);
  if (dst.hasRight()) mAtoms[dst.rightIndex()].setLeftIndex(index);
  ```
- Changed `GAPS_ASSERT(size() > 0)` in `erase()` to
  `GAPS_ASSERT_MSG(size() > 0, "empty AtomicDomain tries to erase an atom")`.
- Removed the now-unused `GAPS_ASSERT(size() > 0)` from `randomFreePosition()`
  (called when domain may legitimately be empty).

### Changed files

| File | Change |
|---|---|
| `src/atomic/AtomicDomain.cpp` | `erase`, `move` rewritten; asserts updated |
| `src/atomic/AtomicDomain.h` | removed `MutableMap`; public accessor added |

Commits: `9dadfc36`, `89100aaf`, `3ebef56c`, `e3f1b601`, `9dadfc36`

---

## 7. All atom index types converted to `size_t`

### Problem

Atom indices were inconsistently typed: `int` in some places, `unsigned` in others,
and `size_t` implied by `std::vector` — making comparisons and casts unsafe and
triggering signed/unsigned warnings.

### Fix

Converted all atom storage indices to `size_t` throughout `Atom.h`,
`AtomicDomain.h`, `AtomicDomain.cpp`.  `random32` → `random64` where an
unsigned 64-bit random position is needed.

Commit: `fdfcd37e`

---

## 8. `static_cast<uint64_t>` of near-max `double` on old Mac Intel Clang

### Problem

`SingleThreadedGibbsSampler` stored domain length as `double mDomainLength` and
converted it back to `uint64_t` via `static_cast`.  Old Apple Clang (x86-64)
converts `double` values very close to `UINT64_MAX` incorrectly — the cast
overflows to `0`.  This silently corrupted the atomic domain length calculation.

Documented with a standalone reproducer: `src/cpp_tests/static_cast_standalone_test.cpp`.

### Fix

Added `uint64_t AtomicDomain::DomainLength() const` accessor that returns
`mDomainLength` directly as `uint64_t`.  Callers now use `mDomain.DomainLength()`
instead of the double round-trip.  The `double mdDomainLength` field in
`SingleThreadedGibbsSampler` was renamed (`mdDomainLength`) to make it distinct
from the safe integer accessor.

### Changed files

| File | Change |
|---|---|
| `src/atomic/AtomicDomain.h` | added `DomainLength()` inline accessor |
| `src/atomic/AtomicDomain.cpp` | minor cleanup |
| `src/gibbs_sampler/SingleThreadedGibbsSampler.h` | use `mDomain.DomainLength()` |
| `src/cpp_tests/static_cast_standalone_test.cpp` | standalone reproducer (not part of test suite) |

Commits: `7f469a30`, `ac8a2e1a`

---

## 9. SIMD NaN in `alphaParameters` on Mac Intel (SSE4/AVX)

Full root-cause analysis is in `.sasha-copilot-context/simd-issue.md`.

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


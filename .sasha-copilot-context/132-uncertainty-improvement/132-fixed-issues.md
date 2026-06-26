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

Already documented in detail in `.sasha-copilot-context/simd-issue.md`.

**Short summary:** `gaps::pmax` left SIMD padding positions at zero.  When
`alphaParameters()` read padding lanes, it computed `0/0 = NaN`, which
poisoned `gibbsMass()` and caused `PSampler` to never accumulate atoms.
Fixed by `padSIMD(min_threshold)` called at the end of both `pmax` overloads.

Fix commit: `bfd5e09d`

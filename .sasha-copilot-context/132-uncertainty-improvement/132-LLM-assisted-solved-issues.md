# LLM-Assisted Solved Issues — Branch `132-uncertainty-improvements`

Branch tracks GitHub issue #132 ("Uncertainty improvements").
The issues below were diagnosed and solved with LLM assistance.
Issue numbering continues from the manually-fixed set (`132-manually-fixed-issues.md`,
issues 1–8); cross-references such as "issue 9" refer to that shared numbering.

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

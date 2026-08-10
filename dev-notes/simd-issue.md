# Bug: SIMD padding zeros cause NaN in alphaParameters, silencing all sampleBirth() calls on Mac Intel

## Affected platforms

Fails on **Mac Intel** (Apple Clang with SSE4.1 or AVX enabled by default).
Passes on Apple M1/M2/M3/M4 (ARM, scalar fallback), MinGW32 (SIMD explicitly
disabled), and Intel/Linux (R toolchain typically does not define `__SSE4_1__`,
so the scalar fallback is used there too).

Failing test: `CoGAPS:::run_catch_unit_tests_by_tag("[densesinglesampler]")`

## Root cause

`DenseNormalModel::alphaParameters()` (and the two overloads
`alphaParameters(r1,c1,r2,c2)` and `alphaParametersWithChange()`) runs a SIMD
loop bounded by `mDMatrix.nRow()`:

```cpp
unsigned size = mDMatrix.nRow();   // e.g. 25
// ...
for (gaps::simd::Index i(0); i < size; ++i)   // increments by SIMD_INC (4 or 8)
{
    pMat.load(mat + i);
    pS.load(S + i);
    gaps::simd::PackedFloat ratio(pMat / (pS * pS));
    partialS     += pMat * ratio;
    partialS_mu  += ratio * (pD - pAP);
}
```

Because `gaps::simd::Index` increments by `SIMD_INC` (4 for SSE4, 8 for AVX),
the last iteration reads up to `SIMD_INC - 1` elements **beyond** the last
real element, i.e. into the SIMD-alignment padding that `Vector` allocates via
`SIMD_PAD(n)`.

`mSMatrix` is the uncertainty matrix.  It is computed in the
`DenseNormalModel` constructor as:

```cpp
mSMatrix = gaps::pmax(mDMatrix, factor, mLambda);
```

`gaps::pmax` (both the `Vector` and `Matrix` overloads) filled only the
`mSize` real elements and left the SIMD padding positions at **0** (the value
set by the `Vector` constructor).

Therefore, when the SIMD loop reads a padding position:

```
pMat[k]  = 0  (padding of the other matrix — also zero)
pS[k]    = 0  (padding of mSMatrix — zero, not filled by pmax)

ratio[k] = pMat[k] / (pS[k] * pS[k]) = 0 / 0 = NaN
partialS[k] += pMat[k] * ratio[k]    = 0 * NaN = NaN
```

NaN propagates through `PackedFloat::scalar()` (which sums all SIMD lanes)
and poisons the returned `AlphaParameters`.

A poisoned `AlphaParameters` causes `gibbsMass()` to return an empty
`OptionalFloat`, so every `sampleBirth()` call is rejected.  As a consequence,
**the P matrix (PSampler) never accumulates any atoms**: `MyMatrix` stays
all-zero, the AP product stays all-zero, and chi-square never decreases.

### Why Mac Intel specifically

On Mac Intel, Apple Clang defines `__SSE4_1__` (or `__AVX__`) by default,
activating the SIMD code path (`SIMD_INC = 4` or `8`).  On the other
platforms:

| Platform | SIMD path | `SIMD_INC` | Bug triggered |
|---|---|---|---|
| Apple M1–M4 (ARM) | no SSE/AVX → scalar | 1 | No |
| MinGW32 | explicitly disabled | 1 | No |
| Intel/Linux (R GCC toolchain) | no `__SSE4_1__` → scalar | 1 | No |
| **Mac Intel (Apple Clang)** | SSE4.1 or AVX | **4 or 8** | **Yes** |

With `SIMD_INC = 1` the loop iterates exactly `size` times and never touches
padding memory.

### Historical note

The original code used `mSMatrix.pad(1.f)`, which set **all** allocated
elements — including the SIMD padding — to `1.0f`.  This was safe because
padding lanes gave `ratio = 0 / 1 = 0`.  When `pad(1.f)` was replaced by
`gaps::pmax(...)` to compute data-driven uncertainty, the SIMD-safety property
was inadvertently lost.

## Fix

Added `padSIMD(float val)` methods to `Vector` and `Matrix` that fill only
the padding positions (`mSize .. mData.size() - 1`) with the given value.
Both `gaps::pmax` overloads now call `padSIMD(min_threshold)` before
returning, so that every SIMD padding element is set to `mLambda > 0`.

With the fix, a padding lane contributes:

```
ratio[k] = 0 / mLambda^2 = 0
partialS[k] += 0 * 0 = 0
```

No NaN, no contribution to the sum — correct behaviour on all platforms.

### Changed files

| File | Change |
|---|---|
| `src/data_structures/Vector.h` | declare `void padSIMD(float val)` |
| `src/data_structures/Vector.cpp` | implement `padSIMD`: fills indices `[mSize, mData.size())` |
| `src/data_structures/Matrix.h` | declare `void padSIMD(float val)` |
| `src/data_structures/Matrix.cpp` | implement `padSIMD`: calls `padSIMD` on each column |
| `src/math/VectorMath.cpp` | `gaps::pmax(Vector,…)` calls `res.padSIMD(min_thr)` |
| `src/math/MatrixMath.cpp` | `gaps::pmax(Matrix,…)` calls `rmat.padSIMD(min_threshold)` |

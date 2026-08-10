# The Uncertainty Model in CoGAPS — Dense vs. Sparse

Reference for how measurement uncertainty (`S`) enters the CoGAPS Gibbs sampler,
and how the two data models — `DenseNormalModel` and `SparseNormalModel` —
compute the same statistics from two very different representations.

Branch `132-uncertainty-improvements`. Written after issues #17, #18, #19 unified
the two models onto one uncertainty formula.

---

## 1. The statistical model

CoGAPS factors a data matrix `D` (genes × samples, or its transpose) as a product
of two non-negative matrices, `D ≈ A · P`. Each observed entry `D[i,j]` is treated
as a Gaussian measurement of the true value `(A·P)[i,j]` with its own standard
deviation `S[i,j]`:

```
D[i,j] ~ Normal( (A·P)[i,j] , S[i,j]^2 )
```

The negative log-likelihood (up to constants) is the chi-square:

```
chiSq = Σ_{i,j} ( D[i,j] − (A·P)[i,j] )^2 / S[i,j]^2
```

`S` is the **uncertainty matrix**. It is an *input assumption*, not something the
sampler learns. Everything in this document is about how `S` is defined and how it
is threaded through the two implementations.

---

## 2. The uncertainty formula (single source of truth)

Both models use the **identical** relative-error model, floored at a constant:

```
S[i,j] = max( factor · D[i,j] , factor )        with   factor = 0.1
```

i.e. a 10 % relative error, but never smaller than `0.1` in absolute terms:

| data value `D`      | uncertainty `S`   | interpretation                         |
|---------------------|-------------------|----------------------------------------|
| `D ≥ 1`             | `0.1 · D`         | relative 10 % error                    |
| `0 < D < 1`         | `0.1` (the floor) | small values are not "super-precise"   |
| `D = 0` (a zero)    | `0.1` (the floor) | zeros share one constant uncertainty   |

The floor is what makes the two models agree and is the crux of issues #17/#19:
without it, a tiny measurement `D = 0.01` would get `S = 0.001` and be trusted
~10 000× more than a zero, which is nonsense. With the floor, every `D < 1`
(zeros included) gets the same `S = 0.1`.

> **Historical note.** Before issue #17 the two models disagreed: the dense floor
> had drifted to `mLambda ≈ 0.006` (issue #2, wrong — `mLambda` is an atom-size
> scale, not an uncertainty), and the sparse model never floored small non-zeros
> at all (a defect inherited from `master`). Issues #17 (alphaParameters) and #19
> (chiSq) put both models back on `S = max(0.1·D, 0.1)`.

---

## 3. Parameter glossary

These names appear throughout the two model classes. Only the first two are about
uncertainty; the rest are listed because they are easy to mistake for it.

### `factor` = 0.1  — the uncertainty coefficient
The relative error **and** the absolute floor: `S = max(factor·D, factor)`.
Defined literally in `DenseNormalModel.h` (`float factor = 0.1f;`); in the sparse
model it is implicit (see `mBeta` below and the `max(d,1)` floor).

### `mBeta` = 1/factor² = 100  — the sparse precision constant *(sparse only)*
`mBeta` is the **precision** (`1/S²`) of a floored entry, i.e. of any entry with
`S = factor`:

```
mBeta = 1 / factor^2 = 1 / 0.1^2 = 100
```

The sparse model factors this constant out of every calculation and multiplies it
back in once at the end. That is *why* the sparse code can avoid ever materialising
an `S` matrix: for the majority of entries (zeros and `D < 1`) the weight `1/S²`
is the *same* number, `mBeta`, so it need not be stored per entry. Only the
non-zero, `D ≥ 1` entries deviate, and they are corrected individually
(section 6). The dense model has no `mBeta`: it divides by the stored `S` directly.

### `invSSq(d)` = 1 / max(d,1)²  — the per-entry weight, factor removed *(sparse only)*
A file-local helper in `SparseNormalModel.cpp`, the single point where the sparse
floor lives. It returns `1/S²` in units where `factor` has been pulled out into
`mBeta`. The full weight of an entry is therefore `mBeta · invSSq(d)`:

```
1/S^2 = 1 / (factor · max(d,1))^2 = (1/factor^2) · (1/max(d,1)^2) = mBeta · invSSq(d)
```

For `d ≥ 1`: `invSSq = 1/d²`. For `d < 1` (and `d = 0`): `invSSq = 1` (the floor,
so the full weight is just `mBeta`).

### `mLambda`  — **NOT uncertainty** (both models)
`mLambda = alpha · sqrt(nPatterns / meanD)` is the scale of the **sparsity /
atom-size prior** (the exponential prior on atom masses in the atomic domain). It
is used in `gibbsMass(..., mLambda)` to shift the proposal mean by `−mLambda`, and
to scale `mMaxGibbsMass`. It has nothing to do with `S`. Conflating it with the
uncertainty floor was exactly the issue #2 mistake.

### `mAnnealingTemp`  — **NOT uncertainty** (both models)
Equilibration temperature; multiplies the `AlphaParameters` (`s`, `s_mu`) during
annealing. Orthogonal to `S`.

### `AlphaParameters { s, s_mu }`  — the Gaussian conditional of one atom
When the sampler proposes changing a single matrix element by a mass `δ`, the
conditional posterior of `δ` is Gaussian, summarised by two numbers:

- **`s`**  — the **precision** (`1/variance`) of that Gaussian,
- **`s_mu`** — the precision-weighted mean (`s · mean`).

`gibbsMass()` turns them into a truncated-normal draw with
`mean = s_mu / s` (or `(s_mu − lambda)/s` with the prior) and `sd = 1/sqrt(s)`.
Both are sums over the affected data column, weighted by `1/S²` — this is the one
place uncertainty enters the *sampling* (as opposed to the *reporting* `chiSq`):

```
s    = Σ_i  mat[i]^2                / S[i]^2
s_mu = Σ_i  mat[i] · (D[i] − AP[i]) / S[i]^2
```

where `mat` is the relevant column of the *other* factor matrix (`P` when updating
`A`, and vice-versa). Derivation: substituting `AP → AP + δ·mat` into `chiSq`
gives a quadratic `chiSq(δ) = const − 2δ·s_mu + δ²·s`, whose coefficients are the
two sums above.

---

## 4. The two representations at a glance

| aspect                       | `DenseNormalModel`                         | `SparseNormalModel`                                   |
|------------------------------|--------------------------------------------|------------------------------------------------------|
| stores `S`?                  | **yes** — full `mSMatrix` (nRow×nCol)      | **no** — never materialised                          |
| how `1/S²` is obtained       | divide by the stored `S`                   | `mBeta · invSSq(d)`, computed on the fly             |
| `factor` lives in            | `factor = 0.1f` → `pmax(D, factor, factor)`| `mBeta = 100` + `max(d,1)` floor in `invSSq`         |
| data storage                 | dense `D`, dense `AP` (cached)             | sparse `D`, `A·P` never cached                        |
| pre-computed helpers         | `mAPMatrix` (cached A·P)                   | `mZ1`, `mZ2` lookup tables (rebuilt in `sync()`)     |
| custom user uncertainty      | **supported** (`setUncertainty`)          | **ignored** (nop — always the default model)         |
| memory footprint             | O(nRow·nCol) — fine for bulk data          | O(non-zeros) — required for single-cell scale        |

The two are kept consistent only by matching math and by the `[sparsegibbs]`
cross-check test — *not* by shared code. Making them share one `factor` constant
is a possible future refactor (see section 8).

---

## 5. Dense model — explicit `S`, direct division

`DenseNormalModel` is the straightforward implementation. It builds the uncertainty
matrix once in the constructor and then just divides by it.

**Construction** (`DenseNormalModel.h`, constructor):
```cpp
float factor = 0.1f;
mSMatrix = gaps::pmax(mDMatrix, factor, factor);   // S = max(factor*D, factor)
```
`pmax(M, mult, floor)` returns `max(mult*M[i,j], floor)`; here both are `factor`,
giving `max(0.1*D, 0.1)`. SIMD tail lanes are padded with `factor` (a positive
value) so later divisions never hit `0/0`.

**chiSq** (`DenseNormalModel.cpp:55`):
```cpp
chisq += ((D[i,j] − AP[i,j]) / S[i,j])^2;     // summed over all i,j
```
Direct transcription of the model. Uses the cached `mAPMatrix`.

**alphaParameters** (`DenseNormalModel.cpp:161`, SIMD):
```cpp
ratio       = mat / (S * S);          // = mat / S^2
partialS    += mat * ratio;           // Σ mat^2 / S^2      -> s
partialS_mu += ratio * (D − AP);      // Σ mat(D−AP) / S^2  -> s_mu
```
- The **2-pattern** overload (`r1==r2`, columns `c1,c2`) replaces `mat` with
  `mat1 − mat2` — used by exchange / two-atom moves.
- **`alphaParametersWithChange(row,col,ch)`** uses residual `D − (AP + ch·mat)`:
  it evaluates the statistics *as if* the element had already been changed by
  `ch`, without mutating `mAPMatrix`.

**Custom uncertainty** (`setUncertainty`): the user may supply their own `S`
matrix; it replaces `mSMatrix` (padded with `1.f` for SIMD safety). This is the
one capability the sparse model does not have.

---

## 6. Sparse model — baseline + correction, `S` never stored

`SparseNormalModel` exists for data (single-cell) where an `nRow×nCol` dense matrix
does not fit in memory. It stores `D` sparsely and **never** forms `S`, `AP`, or
their dense products. Instead it uses the structure of the uncertainty model:

> Every zero entry — and there are many — has the **same** uncertainty
> `S = factor`, hence the same weight `1/S² = mBeta`. So compute everything *as if
> all entries were zero* (a cheap dense-algebra baseline), then **correct** only
> the stored non-zero entries.

### The lookup tables (`mZ1`, `mZ2`), rebuilt in `sync()`
Let `P` denote the *other* matrix (`mOtherMatrix`). Then:
```
mZ1[k]    = Σ_i P[i,k]^2            // squared column norms of P
mZ2[k,l]  = Σ_i P[i,k] · P[i,l]     // Gram matrix P^T P
```
These are the "everything is a zero at `S = factor`" baselines:
- `mZ1[col]` is the baseline `s` (`Σ mat²`, weight 1, i.e. before ×mBeta).
- `mZ2` yields the baseline `s_mu` via `−dot(A_row, mZ2.col(col))`, because
  `Σ_i mat[i]·AP[i] = Σ_k A[row,k]·mZ2[col,k]` (the zeros contribute `−mat·AP`,
  since their `D = 0`).

They are regenerated whenever the other matrix changes — that is what `sync()`
does. (This is why the sampler must interleave `A.update; P.sync(A); P.update;
A.sync(P)`: the sparse tables would otherwise be stale.)

### The per-non-zero correction (`invSSq`)
Iterating the sparse non-zero structure, each stored `D` entry is switched from the
baseline weight `1` to its true floored weight `invSSq(D)`:

**alphaParameters** (`SparseNormalModel.cpp:193`):
```cpp
float invS2 = invSSq(d_val);                 // = 1/max(d,1)^2
s    += v_val * v_val * (invS2 − 1.f);        // fix Σ mat^2 weight: 1 -> invS2
s_mu += v_val*d_val*invS2 + v_val*(1−invS2) * (A_row · P_row);   // fix residual
...
return AlphaParameters(s, s_mu) * mBeta;       // restore the factor
```
The `(invS2 − 1)` / `(1 − invS2)` terms are precisely "remove the baseline
contribution, add the true one". For `d ≥ 1`, `invS2 = 1/d²`; for `d < 1`,
`invS2 = 1`, so the correction is zero — a floored entry already matches the
baseline weight (both are `S = factor`). Multiplying by `mBeta` at the very end
reinstates `factor`.

**chiSq** (`SparseNormalModel.cpp:53`) has the same shape:
```cpp
// baseline: sum A*P^2 over ALL entries (as if every entry were a zero at S=factor)
for all i,j:  chisq += (A·P)[i,j]^2;
// correction: for each stored non-zero D, turn A*P^2 into the floored residual
chisq += d*d*invS2 − 2*d*AP*invS2 + AP*AP*(invS2 − 1);   // = (D−AP)^2*invS2 − AP^2
return chisq * mBeta;
```
Adding the correction to the baseline `AP²` yields `(D − AP)² · invSSq(D)` for each
non-zero, exactly the floored per-entry chi-square. For `d ≥ 1` this is
algebraically identical to the pre-#19 formula `1 + dot(dot − 2d − dsq·dot)/dsq`,
so integer/count data is bit-for-bit unchanged; only `0 < D < 1` entries move.

**No-fit / pre-sync branch** (issue #18): before `sync()`, `mOtherMatrix` is
`NULL` and `A·P = 0`. The baseline loop is zero, so each stored `D` contributes
`(D−0)² · invSSq(D) = D²·invSSq(D)`, `×mBeta`. This matches the dense
`Σ D²/S²` and, for `D ≥ 1`, equals `mBeta` per entry (e.g. `100·nRow·nCol` on
fully-observed integer data). Guarding this branch is what fixed the pre-sync
segfault.

### Why the sparse model cannot just "store an `S` matrix"
Two independent reasons, both fatal:
1. **Memory.** A dense `S` matrix is exactly the footprint the sparse model exists
   to avoid (single-cell matrices are mostly zeros).
2. **Algebra.** The `mZ1`/`mZ2` baseline trick works *only* because the zero-entry
   weight is a single constant that factors out as `mBeta`. A per-entry `S` on the
   zeros would break the factorisation, and the "compute over all, correct the
   non-zeros" shortcut collapses.

So in the sparse model the uncertainty is *computed inline everywhere* (via
`invSSq` + `mBeta`), never stored.

---

## 7. Equivalence and verification

The two models produce the same numbers (up to float ordering) on the same
`A`, `P`, and `D`:

- **Unit test** `[sparsegibbs]` (`testSparseGibbsSampler.cpp`): sets identical
  `A`/`P` on a sparse and a dense sampler over continuous data that includes
  `0 < D < 1`, and asserts `alphaParameters` (1D, 2D, symmetry, with-change) and
  `chiSq` agree to `Approx(0.1 %)`.
- **End-to-end** on GIST: `A·P` reconstruction correlation dense-vs-sparse `= 0.999`
  (equal to a dense-vs-dense different-seed control), and `meanChiSq` for a sparse
  run lands within the dense-run range.
  - Caveat: raw `featureLoadings` correlation is meaningless as a metric — NMF is
    only identifiable up to a permutation of patterns, so even two dense runs
    correlate at ≈ −0.1. Use `A·P` reconstruction or `meanChiSq`.

---

## 8. Known duplication (not yet refactored)

The constant `factor = 0.1` currently lives in **two** places that must stay in
lock-step:

- dense: `float factor = 0.1f;` (→ `pmax(D, factor, factor)`), and
- sparse: `mBeta = 100` (= `1/factor²`) plus the `max(d,1)` floor inside `invSSq`.

Agreement is enforced only by comments and the `[sparsegibbs]` test — the same
fragility that produced issues #17/#19. A cheap, numerics-preserving refactor would
introduce one shared `constexpr float UNCERTAINTY_FACTOR = 0.1f` and derive both
`pmax`'s floor and `mBeta`/`invSSq` from it, so a future change to the model
propagates to both branches automatically. Deferred (does not intersect the
current Phase 3/4 work).

---

## 9. Code map

| what | file:line |
|---|---|
| uncertainty formula (dense build) | `src/gibbs_sampler/DenseNormalModel.h:96,110` |
| dense chiSq | `src/gibbs_sampler/DenseNormalModel.cpp:55` |
| dense alphaParameters (1D / 2D / with-change) | `src/gibbs_sampler/DenseNormalModel.cpp:161,185,216` |
| dense custom uncertainty | `src/gibbs_sampler/DenseNormalModel.h:113` (`setUncertainty`) |
| `invSSq` (single sparse floor) | `src/gibbs_sampler/SparseNormalModel.cpp:25` (doc block from :17) |
| `mBeta` init (= 100) | `src/gibbs_sampler/SparseNormalModel.h:77` |
| sparse chiSq (baseline + correction + no-fit) | `src/gibbs_sampler/SparseNormalModel.cpp:53` |
| sparse alphaParameters (1D / with-change / 2D) | `src/gibbs_sampler/SparseNormalModel.cpp:193,237,283` |
| `mZ1`/`mZ2` lookup tables | `src/gibbs_sampler/SparseNormalModel.cpp:337` (`generateLookupTables`) |
| `s`/`s_mu` → truncated-normal draw | `src/gibbs_sampler/AlphaParameters.cpp:27,38` (`gibbsMass`) |
| cross-check test | `src/cpp_tests/testSparseGibbsSampler.cpp` (`[sparsegibbs]`) |

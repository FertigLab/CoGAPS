# Annotation-weighted subset sampling — problem report and options

Status: **open design issue**, deliberately left out of branch
`132-uncertainty-improvements` (which is about the uncertainty model). This
report documents the problem so the decision and evidence are not lost.

Surfaced while tidying the R `testthat` suite (Phase 4): the test
`test_subset_data.R` "subsetting data with annotation weights" is entirely
commented out, with a TODO. Investigating why it is commented out revealed a
genuine correctness gap in the feature it was meant to cover.

---

## 1. What the feature is

For **distributed** CoGAPS (`GWCoGAPS` / `scCoGAPS`, i.e. `distributed =
"genome-wide"` or `"single-cell"`), the data is partitioned into `nSets`
subsets that are factored in parallel and then stitched back together.

`setAnnotationWeights(params, annotation, weights)` lets the user bias that
partitioning. Every gene (genome-wide) or sample (single-cell) carries a
category label from `annotation`; every category carries a number in `weights`.
Subsets are then drawn so that heavily-weighted categories are over-represented,
instead of the default uniform partition.

```r
weight <- c(A = 1, B = 2, C = 3)
anno   <- sample(names(weight), nrow(data), replace = TRUE)
params <- setAnnotationWeights(CogapsParams(), anno, weight)
res    <- CoGAPS(data, params, distributed = "genome-wide", ...)
```

The two values are stored in the `CogapsParams` slots `samplingAnnotation` and
`samplingWeight` (`R/class-CogapsParams.R`), set together by
`setAnnotationWeights` (`R/methods-CogapsParams.R`).

---

## 2. How it works today

`createSets()` (`R/SubsetData.R`) dispatches to `sampleWithAnnotationWeights()`
whenever `samplingAnnotation` is set:

```r
sampleWithAnnotationWeights <- function(allParams, setSize)
{
    weight <- allParams$gaps@samplingWeight
    weight <- weight[order(names(weight))]
    groups <- sort(unique(allParams$gaps@samplingAnnotation))

    lapply(1:allParams$gaps@nSets, function(i)
    {
        # 1. draw setSize category tickets, proportional to weight
        groupCount <- sample(groups, size = setSize, replace = TRUE, prob = weight)
        # 2. for each category, draw that many member indices WITH REPLACEMENT
        sort(unlist(sapply(groups, function(g)
        {
            groupNdx <- which(allParams$gaps@samplingAnnotation == g)
            sample(groupNdx, size = sum(groupCount == g), replace = TRUE)
        })))
    })
}
```

Both draws use `replace = TRUE`. That is the root of the problem: step 2 can pick
the same gene several times within one subset, and nothing links subsets, so a
gene may land in several subsets or in none.

---

## 3. The defects (empirically confirmed)

Reproduction — GIST (1363 genes), three categories A/B/C with weights 1/2/3,
`distributed = "genome-wide"`, `nSets = 4`:

| observation | value | should be |
|---|---|---|
| duplicates **within** each subset | ~75 of 340 (**~22 %**) | 0 (a gene at most once per subset) |
| unique genes covered across all subsets | **597 of 1363** (~44 %) | all genes represented |
| `nrow(res@featureLoadings)` | **1360** (= 4 × 340, i.e. counts duplicates) | 1363 (one row per gene) |

So a weighted run (a) puts the same gene in a subset multiple times, (b) leaves
~56 % of genes unsampled, and (c) emits duplicated genes as **separate rows** in
`featureLoadings` instead of collapsing them. The output matrix is therefore
mis-shaped (rows ≠ genes) and its rows are not unique.

This matches the original author's TODO, verbatim from the commented test:

> address how weighted sampling works with duplicates, do we need to allow
> passing a value for setSize in this case? we should collapse down using the
> mean; prevent multiple copies from being in the same set

## 4. Why the existing test is commented out

The commented test asserts the *fixed* behaviour:

```r
expect_equal(nrow(result@featureLoadings), nrow(GIST.matrix))          # 1363
expect_equal(sum(sapply(sets, length)), nrow(result@featureLoadings))  # full coverage
```

Against the current implementation the first assertion is false (`1360 ≠ 1363`),
so the test cannot pass without changing the feature. It was commented out rather
than fixed — leaving a silent gap: the `sampleWithAnnotationWeights` path has
**zero** test coverage.

---

## 5. Options

### A. Delete the stub, file the design gap as a tracked issue  *(recommended for this branch)*
Remove the commented-out test and open an issue capturing §3 (the quantified
defect) and §6 (a fix sketch). Rationale: the test describes behaviour the
feature does not yet provide; making it pass is a feature change out of scope for
the uncertainty branch. A commented-out test is worse than none — it reads as
"covered" while testing nothing.
- Pro: keeps branch 132 focused; the gap is recorded, not lost.
- Con: the sampling path stays untested until the feature is fixed.

### B. Add a characterization test of the *current* behaviour
`sampleWithAnnotationWeights` has no tests at all. Add one that asserts only what
is true today — it runs without error, returns `nSets` subsets, and the
high-weight category is over-represented relative to a uniform draw — with an
explicit comment (and a filed issue) that duplicate handling is a known gap.
- Pro: real regression coverage of an untested code path, cheaply.
- Con: risks blessing buggy behaviour; must be clearly labelled as
  characterization, not correctness.

### C. Fix the feature, then test it  *(out of scope here)*
Implement the intended semantics and test the fixed behaviour. This is a feature
change to the distributed subsetting and result assembly, not test tidying, and
belongs in its own branch/PR — not in the uncertainty branch.

**Recommendation:** **A** for branch 132 (delete + issue). If we want to avoid
losing coverage entirely, **B** is an acceptable middle ground. **C** is a
separate piece of work.

---

## 6. Sketch of a real fix (for the issue / option C)

Not implemented here; recorded so the issue is actionable.
1. **Dedupe within a subset.** Draw category tickets as now, but draw member
   indices **without replacement** per category (`replace = FALSE`), capping the
   count at the category size. Decide the policy when a category is exhausted
   (spill to other categories, or shrink the subset).
2. **Collapse duplicates by mean.** If a gene still ends up represented more than
   once (across the stitched result), average its rows rather than emitting
   copies — as the TODO says ("collapse down using the mean").
3. **Coverage / `setSize`.** Clarify the contract: is weighted sampling meant to
   cover every gene (a weighted *partition*) or to draw a weighted *sample*
   (coverage < 100 % by design)? The current output shape (rows ≠ genes) implies
   the former was intended. Possibly expose `setSize` so the caller controls it.
4. Re-enable the test to assert the chosen contract (unique rows, expected shape,
   category proportions within tolerance).

---

## 7. Code map

| what | location |
|---|---|
| params slots `samplingAnnotation` / `samplingWeight` | `R/class-CogapsParams.R` (slot docs ~27–29) |
| `setAnnotationWeights` generic / method | `R/class-CogapsParams.R`; `R/methods-CogapsParams.R` |
| weighted draw (the defect) | `R/SubsetData.R` — `sampleWithAnnotationWeights()` |
| dispatch into it | `R/SubsetData.R` — `createSets()` (the `samplingAnnotation` branch) |
| commented-out test (the stub) | `tests/testthat/test_subset_data.R` — "subsetting data with annotation weights" |

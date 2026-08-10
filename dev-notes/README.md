# CoGAPS developer notes — index

Nothing in this directory is part of the package: `dev-notes/` is excluded by
`.Rbuildignore`, so it reaches neither the tarball nor `R CMD check` nor
Bioconductor. It is kept in the repository so that the reasoning behind the code
travels with it.

One file describes the way we work rather than any particular piece of work:
[`rus/.agent-rules-rus.md`](rus/.agent-rules-rus.md), the working agreements with
the maintainer. It is deliberately in Russian and lives in `rus/` for that reason;
everything an outside reader sees stays in English. The assistant-facing
description of the project — how it is built, tested and laid out — is
`../CLAUDE.md`. Everything else below belongs to branch
`132-uncertainty-improvements`, which tracks GitHub issue #132 ("Uncertainty
improvements").

## Defect journals

Fixed defects share one issue numbering across the first two files: 1–8 were fixed
manually, 9 onwards with LLM assistance. A cross-reference such as "issue 8" always
refers to that shared numbering.

| File | What it is |
|---|---|
| [`132-manually-fixed-issues.md`](132-manually-fixed-issues.md) | Issues 1–8, fixed by hand by @favorov: the `pmax` signature split, data-driven uncertainty in `DenseNormalModel`, the `static_cast<uint64_t>` overflow, and so on. |
| [`132-LLM-assisted-solved-issues.md`](132-LLM-assisted-solved-issues.md) | Issues 9–20, diagnosed and fixed with LLM assistance — SIMD padding NaN, `SIMD_PAD` on ARM, the stale iterator in `AtomicDomain::move()`, checkpoint deserialization, empty-container guards, the two models brought onto one uncertainty formula (17–19), and `Vector::pad()` silently discarding a user-supplied `uncertainty=` matrix (20). The longest file here and the best entry point for why the C++ looks the way it does. |
| [`master-issues.md`](master-issues.md) | Defects inherited from `master` and fixed here even though they have nothing to do with uncertainty — found while merging `master` in (2026-08) and while covering exported functions that no test had ever called. Records where they came from, so a reviewer does not have to wonder. |

## Removal of the asynchronous sampler

Read in this order; each file was written before the next step was taken.

| File | What it is |
|---|---|
| [`remove-async-plan-eng.md`](remove-async-plan-eng.md) | The spec: why the OpenMP sampler breaks MCMC detailed balance, and the full list of edits its removal requires. |
| [`async-removal.md`](async-removal.md) | The retrospective report of executing that spec. |
| [`after_async-removed-plan.md`](after_async-removed-plan.md) | The audit that followed: does the C++ suite still compile, which test cases are empty stubs, and a ranked list of latent bugs in the surviving code. Issues 9 onwards grow out of this list. |

## Reference write-ups

| File | What it is |
|---|---|
| [`uncertainty-model-eng.md`](uncertainty-model-eng.md) | How measurement uncertainty `S` enters the sampler, and how `DenseNormalModel` and `SparseNormalModel` compute the same statistics from two very different representations. Written once issues 17–19 had unified them onto one formula. |
| [`simd-issue.md`](simd-issue.md) | The standalone root-cause analysis behind issue 9: SIMD padding zeros produce `NaN` in `alphaParameters()`, which silences every `sampleBirth()` call on Mac Intel. |

## Open questions

| File | What it is |
|---|---|
| [`annotation-weights-sampling-issue-eng.md`](annotation-weights-sampling-issue-eng.md) | Annotation-weighted subset sampling for distributed CoGAPS is incorrect. Deliberately **not** fixed on this branch; the report keeps the evidence and the options for whoever takes it on. |

## Reproducer

| Directory | What it is |
|---|---|
| [`static-cast-uint64-reproducer/`](static-cast-uint64-reproducer/) | A standalone program for issue 8: `static_cast<uint64_t>` of a `double` near `UINT64_MAX` overflows to `0` on old Apple Clang. Has its own `main()`, is not listed in `configure.ac`, and is never compiled by the build — it lives here rather than in `src/cpp_tests/`, where it looked like a unit test. |

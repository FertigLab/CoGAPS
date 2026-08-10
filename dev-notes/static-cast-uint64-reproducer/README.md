# `static_cast<uint64_t>` of a near-max `double` — standalone reproducer

This directory holds the reproducer written while investigating issue 8
(`132-manually-fixed-issues.md` §8). It is **not** part of the package: it has its
own `main()`, is not listed in `configure.ac`, and is never compiled by the build.
It used to sit in `src/cpp_tests/`, where it looked like a unit test and shipped
inside the package tarball; it lives here instead so that `src/` contains only
code that is actually built.

## What it demonstrates

`SingleThreadedGibbsSampler` used to keep the atomic domain length as a `double`
and convert it back with `static_cast<uint64_t>`. On old Apple Clang (x86-64),
casting a `double` whose value is very close to `UINT64_MAX` does not yield the
neighbouring integer — it overflows to `0`. The atomic domain length was then
silently wrong.

The program prints `UINT64_MAX` and `UINT64_MAX - 10`, converts both to `double`
and back, and shows the round-trip result. On an affected toolchain the values
come back as `0` instead of the expected magnitude; note that `UINT64_MAX - 10`
is not exactly representable as a `double` in the first place, so the round trip
cannot be lossless even where it does not overflow.

## How to run it

It is deliberately self-contained — no R, no package headers:

```
c++ static_cast_standalone_test.cpp -std=gnu++17 -O0
./a.out
```

On an unaffected toolchain (checked on Apple arm64) both values round-trip to
`18446744073709551615`, i.e. `UINT64_MAX - 10` rounds *up* to `UINT64_MAX`:

```
uint64_t value 1 (max) : 18446744073709551615
uint64_t value 2 (max-10) : 18446744073709551605
value 1 as double : 1.84467e+19
value 2 as double : 1.84467e+19
value 1 as uint64_t casted back from double : 18446744073709551615
value 2 as uint64_t casted back from double : 18446744073709551615
```

That is the *benign* outcome — lossy, but the magnitude survives. The bug is the
last two lines coming back as `0`. So running this on arm64 or on Linux/GCC will
not reproduce it; an old Apple Clang on x86-64 is needed.

## The fix in the package

`AtomicDomain` gained a `uint64_t DomainLength() const` accessor that returns the
length as an integer directly, and the callers use `mDomain.DomainLength()`
instead of the `double` round trip. The remaining `double` field was renamed to
`mdDomainLength` so that the unsafe value is visibly distinct from the safe
accessor; it is still needed for the birth/death probability arithmetic, which is
genuinely floating point (`SingleThreadedGibbsSampler.h`).

See `132-manually-fixed-issues.md` §8 for the full write-up and the commits.

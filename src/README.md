# Some notes about the structures CoGAPS uses

## A(P)Sampler
In the GapsRunner, all the sampling events are done by two samplers, one for a decomposition matrix, 
$D=AP$ ASampler and PSampler. The type of the sampler objects depends on the data model it uses (Dense or Sparse).
The data is organised in the sampler as folows.

Let's say that there are $l$ rows (usually, genes) in $D$ matrix, $m$ (samples) columns in $D$ and the decomposition runs with $k$ patterns. So, $D$ and $AP$ are both of $m \times l$ size.

Each sampler carries its own copy of $AP$, they are copied by the `sync()` method and they are calculated by the `extraInitialization()`. Each sampler has it purpose matrix, and the uncertanoty matrix. The latter is og the same size as the AP. Each sampler has a `const sampler &` reference to its *vis-a-vis*.

A sampler can be transposed or not. The run involves a transposed and a nontransposes sampler. If the input D matrix is not transposed by the call, the left (A) sampler is transposed, the right (P) is not.

A transposed (left, A, if the D is non-transposed) carries a transposed AP matrix of $l \times m$ and the puspose (A) matrix of $m \times k$ size.

A non-transposed (right, P, if the D is non-transposed) carries a nontransposed AP matrix of $m \times l$ and the puspose (P) matrix of $l \times k$ size.

So, for any correct sampler, nrows(MyMatrix)==ncols(APMatrix); ncols(MyMatrix) is the pattern munber, and APMatrix has the same dimesions as transposed other\_sampler.APMatrix. After `sync()` APMatrix==tr(other\_sampler.APMatrix).

## Debug, etc congigure options
All the options for config are given in ../configure.ac

After changing, regenerate `configure`. Plain `autoconf` is **not** enough:
`configure.ac` uses `AX_COMPILER_VENDOR` and `AX_COMPILER_VERSION` from
autoconf-archive, and those are pulled in by `aclocal`, not by `autoconf`. Run
both, from the package root:

```
aclocal -I /opt/homebrew/share/aclocal   # path to the autoconf-archive macros
autoconf
```

If you skip `aclocal`, the two `AX_*` macros are left unexpanded and end up in
`configure` as literal shell commands. It still "works" — you just get

```
./configure: line 2940: AX_COMPILER_VENDOR: command not found
building on  compiler version
```

and `$ax_cv_cxx_compiler_vendor` stays empty, which silently disables the
vendor-dependent branch: `--enable-warnings` then adds no flags at all. The
`configure` on `master` has exactly this problem; the one on this branch was
regenerated correctly and does not.

`aclocal.m4` is a build artefact of that procedure and is not committed.

To pass and option to configure when running devtools::load_all or similar, set the environment varianle, for example, to run ./configure --enable-debug, run the following in R:

Sys.setenv(enable_debug="yes")
devtools::load_all(recompile = TRUE) 

We suppose we are in the root of the package.
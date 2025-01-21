# Some notes about the structures CoGAPS uses

## A(P)Sampler
In the GapsRunner, all the sampling events are done by two samplers, one for a decomposition matrix, 
$D=AP$ ASampler and PSampler. The type of the sampler objects depends on the data model it uses (Dense or Sparse).
The data is organised in the sampler as folows.

Let's say that there are $l$ rows (usually, genes) in $D$ matrix, $m$ (samples) columns in $D$ and the decomposition runs with $k$ patterns. So, $D$ and $AP$ are both of $m \times l$ size.

Each sampler carries its own copy of $AP$, they are copied by the `sync()` method and they are calculated by the `extraInitialization()`. Each sampler has it purpose matrix, and the uncertanoty matrix. The latter is og the same size as the AP. Each sampler has a `const sampler &` reference to its *vis-a-vis*.

A sapler can be transposed or not. The run ivolves a transposed and a nontransposes sampler. If the input D matrix is not transposed by the call, the left (A) sampler is transposed, the right (P) is not.

A transposed (left, A, if the D is nontransposed) carries a transposed AP matrix of $l \times m$ and the puspose (A) matrix of $m \times k$ size.

A nontransposed (right, P, if the D is nontransposed) carries a nontransposed AP matrix of $m \times l$ and the puspose (P) matrix of $l \times k$ size.

So, for any correct sampler, nrows(MyMatrix)==ncols(APMatrix); ncols(MyMatrix) is the pattern munber, and APMatrix has the same dimesions as transposed other\_sampler.APMatrix. After `sync()` APMatrix==tr(other\_sampler.APMatrix).
 

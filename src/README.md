# Some notes about the structures CoGAPS uses

## A(P)Sampler
In the GapsRunner, all the sampling events are done by two samplers, one for a decomposition matrix, 
$D=AP$ ASampler and PSampler. The type of the sampler objects depends on the data model it uses (Dense or Sparse).
The data is organised in the sampler as folows.

Let's say that there are $l$ rows (usually, genes) in $D$ matrix, $m$ (samples) columns in $D$ and the decomposition runs with $k$ patterns. So, $D$ is $m \times l$

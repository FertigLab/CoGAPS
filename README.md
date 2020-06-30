# CoGAPS (Version: 3.7.0)

[![Bioc](https://bioconductor.org/images/logo_bioconductor.gif)](https://bioconductor.org/packages/CoGAPS)
[![downloads](https://bioconductor.org/shields/downloads/release/CoGAPS.svg)](http://bioconductor.org/packages/stats/bioc/CoGAPS/)
[![Build Status](https://travis-ci.org/FertigLab/CoGAPS.svg?branch=master)](https://travis-ci.org/FertigLab/CoGAPS)

Coordinated Gene Activity in Pattern Sets (CoGAPS) implements a Bayesian MCMC matrix factorization algorithm, GAPS, and links it to gene set statistic methods to infer biological process activity. It can be used to perform sparse matrix factorization on any data, and when this data represents biomolecules, to do gene set analysis.

# Installing CoGAPS

*CoGAPS* is a bioconductor R package ([link](https://bioconductor.org/packages/CoGAPS)) and so the release version can be installed
as follows:

```
install.packages("BiocManager")
BiocManager::install("CoGAPS")
```

The most up-to-date version of *CoGAPS* can be installed directly from the 
*FertigLab* Github Repository:

```
BiocManager::install("FertigLab/CoGAPS")
```

# Using CoGAPS

Follow the vignette [here](https://www.bioconductor.org/packages/devel/bioc/vignettes/CoGAPS/inst/doc/CoGAPS.html)

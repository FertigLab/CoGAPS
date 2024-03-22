<img src="https://user-images.githubusercontent.com/25310425/169565420-56958b50-29a2-4032-afb3-08447577d074.png" width="150">

[![R build status](https://github.com/FertigLab/CoGAPS/workflows/r-build-check/badge.svg)](https://github.com/FertigLab/CoGAPS/actions?workflow=r-build-check)

# CoGAPS

[![Bioc](https://bioconductor.org/images/logo_bioconductor.gif)](https://bioconductor.org/packages/CoGAPS)
[![downloads](https://bioconductor.org/shields/downloads/release/CoGAPS.svg)](http://bioconductor.org/packages/stats/bioc/CoGAPS/)

Coordinated Gene Activity in Pattern Sets (CoGAPS) implements a Bayesian MCMC matrix factorization algorithm, GAPS, and links it to gene set statistic methods to infer biological process activity. It can be used to perform sparse matrix factorization on any data, and when this data represents biomolecules, to do gene set analysis.

# Installing CoGAPS

Via Bioconductor:

```
install.packages("BiocManager")
BiocManager::install("FertigLab/CoGAPS")
```

The most up-to-date version of *CoGAPS* can be installed directly from the
*FertigLab* Github Repository:

```
devtools::install_github("FertigLab/CoGAPS")
```

# Using CoGAPS

Follow the vignette [here](https://github.com/FertigLab/CoGAPS/blob/master/vignettes/CoGAPS.Rmd) and available as static html [here](https://rpubs.com/jeanettejohnson/1018399)

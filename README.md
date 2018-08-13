# CoGAPS Version: 3.3.10

[![Bioc](https://bioconductor.org/images/logo_bioconductor.gif)](https://bioconductor.org/packages/CoGAPS)
[![downloads](https://bioconductor.org/shields/downloads/CancerInSilico.svg)](https://bioconductor.org/packages/CoGAPS)

Coordinated Gene Activity in Pattern Sets (CoGAPS) implements a Bayesian MCMC matrix factorization algorithm, GAPS, and links it to gene set statistic methods to infer biological process activity. It can be used to perform sparse matrix factorization on any data, and when this data represents biomolecules, to do gene set analysis.

# Installing CoGAPS

*CoGAPS* is a bioconductor R package and so the release version can be installed
as follows:

```
source("https://bioconductor.org/biocLite.R")
biocLite("CoGAPS")
```

The most up-to-date version of *CoGAPS* can be installed directly from the 
*FertigLab* Github Repository:

```
## Method 1 using biocLite
biocLite("FertigLab/CoGAPS", dependencies = TRUE, build_vignettes = TRUE)

## Method 2 using devtools package
devtools::install_github("FertigLab/CoGAPS")
```

There is also an option to install the development version of *CoGAPS*, 
while this version has the latest experimental features, it is not guaranteed
to be stable.

```
## Method 1 using biocLite
biocLite("FertigLab/CoGAPS", ref="develop", dependencies = TRUE, build_vignettes = TRUE)

## Method 2 using devtools package
devtools::install_github("FertigLab/CoGAPS", ref="develop")
```

# Using CoGAPS

Follow the vignette here: http://htmlpreview.github.io/?https://github.com/FertigLab/CoGAPS/blob/develop/vignettes/CoGAPS.html
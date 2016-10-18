---
title: "GWCoGAPS and PatternMarkers Vignette"
author: "Genevieve L. Stein-O'Brien"
date: \today
output: BiocStyle::pdf_document
bibliography: AppNote.bib
vignette: >
  %\VignetteIndexEntry{Overview of CNPBayes package}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc} 
---

# Introduction

NMF algorithms associate gene expression changes with biological processes (e.g., time-course dynamics or disease subtypes). Compared with univariate gene associations, the relative gene weights of NMF solutions do not clearly identify gene biomarkers for these processes. Therefore, we developed a novel PatternMarkers statistic to extract genes uniquely representative of NMF patterns for enhanced visualization and subsequent biological validation. Finding unbiased gene markers with PatternMarkers requires whole-genome data. However, NMF algorithms typically do not converge for the tens of thousands of genes in genome-wide profiling. Therefore, we also developed GWCoGAPS to simultaneously cut runtime and ensure convergence for whole genome NMF with a sparse, MCMC NMF algorithm CoGAPS. The software includes a shiny application PatternMatcher to compare NMF patterns. 

# Executing GWCoGAPS

In this chapter, we describe how to run both the GWCoGAPS algorithm and a manual
pipeline for genome wide NMF analysis.

## GWCoGAPS

### Methods
The GWCoGAPS algorithm is run by calling the GWCoGAPS function in the CoGAPS
R package as follows:


```r
library(CoGAPS) 
GWCoGAPS(D, S, nFactor, nSets, nCores, saveBySetResults, fname,
    PatternsMatchFN = patternMatch4Parallel, Cut, minNS, ...)
```

**Input Arguments**
The inputs that must be set each time are only the nSets, nFactor, and data and standard deviation matrices, with all other inputs having default values. Additional inputs to the gapsRun function, as outlined previously, can be used to taylor the analysis based on the expected dimensionality of the data. GWCoGAPS specific arguments are as follows:

\begin{description}
\item[D]{data matrix}
\item[S]{uncertainty matrix (std devs for chi-squared of Log Likelihood)}
\item[nFactor]{number of patterns (basis vectors, metagenes), which must be
greater than or equal to the number of rows of FP}
\item[nSets]{number of sets for parallelization}
\item[nCores]{number of cores for parallelization. If left to the default NA, nCores = nSets. }
\item[saveBySetResults]{logical indicating whether to save by intermediary by set results. Default is FALSE.}
\item[fname]{character string used to label file output. Default is "GWCoGAPS.out"}
\item[PatternsMatchFN]{function to use for pattern matching across sets. Default is patternMatch4Parallel.}
\item[Cut]{number of branches at which to cut dendrogram used in patternMatch4Parallel}
\item[minNS]{minimum of individual set contributions a cluster must contain used in patternMatch4Parallel}
\end{description}

\par Once the GWCoGAPS algorithm has been run, the inferred patterns and corresponding amplitudes can processed in the same way as output from gapsRun, gapsMapsRun, or CoGAPS.

### Example: Simulated data

\par In this example, we use the same simulated data in SimpSim (SimpSim.D), as previously described, with three known patterns (SimpSim.P) and corresponding amplitude (SimpSim.A) with specified activity in two gene sets (GSets).


```r
library('CoGAPS')
data('SimpSim')
GWCoGAPS(SimpSim.D, SimpSim.S, nFactor=3, nCores=NA, nSets=2,
            fname="test" ,PatternsMatchFN = postCoGAPSPatternMatch,
            sampleSnapshots = "TRUE", numSnapshots = 3)
plotGAPS(AP.Fixed$A, AP.Fixed$P, 'ModSimFigs')
```


## Manual Pipeline

The functions that compose the core of the GWCoGAPS algorithm are also provided in the CoGAPS R package and can be used for a manual pipeline of the same analysis or to create a custom pipeline. Additionally, a web based shiny application is included for visualizing and editing the pattern matching process.

### PatternMatcher Shiny App

**Input Arguments**
The inputs that must be set each time is the PBySet. Additionally, the order and sample.color inputs are only available when calling patternMatcher in an interactive R session.

\begin{description}
\item[PBySet]{list of matched set solutions for the Pmatrix from an NMF algorithm}
\item[out]{optional name for saving output}
\item[order]{optional vector indicating order of samples for plotting. Default is NULL.}
\item[sample.color]{optional vector of colors of same lenght as colnames. Default is NULL.}
\end{description}

![A screenshot of the PatternMatcher Shiny App plotting the by-set results for pattern 1 of the simulated data](ShinyApp.pdf)

### Example: Simulated data}

\par This example will manually generate the same result as calling the GWCoGAPS function as given in the example in the previous section.


```r
data('SimpSim')
D<-SimpSim.D
S<-SimpSim.S
nFactor=3
nSets<-2
nCores=NA
numSnapshots=3
saveBySetResults=TRUE
fname="test"
minNS=NA

# break the data into sets
genesInSets<-createGWCoGAPSSets(data=D, nSets=nSets, keep=FALSE)

#generate seeds for parallelization
nut<-generateSeeds(chains=nSets, seed=-1)

# establish the number of cores that you are able to use
if(is.na(nCores)){nCores<-nSets}
registerDoParallel(cores=nCores)

# run CoGAPS for each set
AP <- foreach(i=1:nSets) %dopar% {
  D <- as.matrix(D[genesInSets[[i]],])
  S <- as.matrix(S[genesInSets[[i]],])
  gapsRun(D=D, S=S, nFactor=nFactor,seed=nut[i],numSnapshots=numSnapshots)
}

BySet<-reOrderBySet(AP=AP,nFactor=nFactor,nSets=nSets)
matchedPs<-patternMatch4Parallel(Ptot=BySet$P,nP=nFactor,nS=nSets,cnt=nFactor,minNS=minNS,bySet=TRUE)

# use shiny for pattern matching
selectPBySet<-PatternMatcher(PBySet=matchedPs[["PBySet"]])
library(reshape2)
selectPBySet<-dcast(selectPBySet,  BySet ~ Samples)
rownames(selectPBySet)<-selectPBySet$BySet
selectPBySet<-as.matrix(selectPBySet[,-1])
matchedPs<-patternMatch4Parallel(Ptot=selectPBySet,nP=nFactor,nS=nSets,cnt=nFactor,minNS=minNS,bySet=FALSE)

#generate seeds for parallelization
nut<-generateSeeds(chains=nSets, seed=-1)
#final number of factors
nFactorFinal<-dim(matchedPs)[1]

# run fixed CoGAPS
Fixed <- foreach(i=1:nSets) %dopar% {
  D <- as.matrix(D[genesInSets[[i]],])
  S <- as.matrix(S[genesInSets[[i]],])
  AP <- gapsMapRun(D, S, FP=matchedPs, nFactor=nFactorFinal, fixedMatrix = "P",seed=nut[i],numSnapshots=numSnapshots)
}

#extract A and Asds
As4fixPs <- postFixed(AP.fixed=Fixed,setPs=matchedPs)
```


# PatternMarkers

In this chapter, we describe the PatternMarkers statistic.

## PatternMarkers

The PatternMarkers statistic finds the genes most uniquely associated with a given pattern or linear combination of patterns by computing

\begin{equation}
\sqrt{\left(\bf{A}_{i}-l_{p}\right)^{\textit{t}} \left(\bf{A}_{i}-l_{p}\right)}
\label{eq:PatternMarkers}
\end{equation}

where  are the elements of the $\bf{A}$ matrix for the $\textit{i}^{th}$ gene scaled to have a maximum of one and \textit{l} is the $\textit{p}^{th}$ user specified norm. The genes are then ranked. In the case where \textit{l} is the identity vector, the ranking is run separately for each of the K patterns. Unique sets are generated by thresholding using the first gene to have a lower ranking, i.e. better fit to, another patterns.

### Methods

The PatternMarkers statistic is run by calling thepatternMarkers function in the CoGAPS R package as follows:


```r
 patternMarkers(Amatrix = AP$Amean, scaledPmatrix = FALSE, Pmatrix = NA,
  threshold = "cut", lp = NA, full = FALSE, ...)
```

**Input Arguments**
\begin{description}
\item[Amatrix]{A matrix of genes by weights resulting from CoGAPS or other NMF decomposition}
\item[scaledPmatrix]{logical indicating whether the corresponding pattern matrix was fixed to have max 1 during decomposition}
\item[Pmatrix]{the corresponding Pmatrix (patterns X samples) for the provided Amatrix (genes x patterns). This must be supplied if scaledPmatrix is FALSE.}
\item[threshold]{the type of threshold to be used. The default "cut" will thresholding by the first gene to have a lower ranking, i.e. better fit to, a pattern. Alternatively, threshold="all" will return all of the genes in rank order for each pattern.}
\item[lp]{a vector of weights for each pattern to be used for finding markers. If NA markers for each pattern of the A matrix will be used.}
\item[full]{logical indicating whether to return the ranks of each gene for each pattern}
\end{description}

Once the PatternMarkers statistic has been run, a heatmap of each markers expression level can be displayed using the plotPatternMarkers function as follows:


```r
plotPatternMarkers(data = NA, patternMarkers = PatternMarkers,
      patternPalette = NA, sampleNames = NA, samplePalette = NA,
      colDenogram = TRUE, heatmapCol = "bluered", scale = "row", ...)
```
 
**Input Arguments**
\begin{description}
\item[data]{the dataset from which the patterns where generated}
\item[patternMarkers]{the list of genes generated from the patternMarkers function}
\item[patternPalette]{a vector indicating what color should be used for each pattern}
\item[sampleNames]{names of the samples to use for labeling }
\item[samplePalette]{a vector indicating what color should be used for each sample}
\end{description}

### Example: Simulated data


```r
# GWCoGAPS output
PatternMarkers<-patternMarkers(Amatrix=AP.fixed$A,scaledPmatrix=TRUE,threshold="cut")
plotPatternMarkers(data=D,patternMarkers=PatternMarkers,patternPalette=c("grey","navy","orange"))
# manual pipeline output
PatternMarkers<-patternMarkers(Amatrix=As4fixPs$A,scaledPmatrix=TRUE,threshold="cut")
plotPatternMarkers(data=D,patternMarkers=PatternMarkers,patternPalette=c("grey","navy","orange"))
```

Figure \ref{fig:PM1} shows the results from running plotPatternMarkers on the PatternMarkers generated from the the GWCoGAPS (left) and manual pipeline (right) results generated from the simulated data as previously illustrated.
\begin{figure}[h]
\begin{center}
\includegraphics[width=0.4\linewidth]{GWCoGAPSPMs}
\includegraphics[width=0.4\linewidth]{ManPipePMs}
\end{center}
\caption{Heatmap of PatternMarkers expression levels for simulated data.}
\label{fig:PM1}
\end{figure}

# Feedback

Please send feedback to Genevieve L. Stein-O'Brien \texttt{gsteino1@jhmi.edu}.

If you want to send a bug report, please first try to reproduce the error.  The code is stochastic and we have seen many transient errors arising in the boost libraries which rarely repeat.  Send the data and please describe what you think should have happened, and what did happen.

\bibliographystyle{plain}
\bibliography{AppNote}



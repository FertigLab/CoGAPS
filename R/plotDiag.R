#heatmaps and line plots

#Input command for R, save cogaps results to variable "results"
#Plots(results$AMean, results$PMean, results$AStd, results$PStd, results$ChiSq, results$AtomsAEquil, results$AtomsASamp, results$AtomsPEquil, results$AtomsPSamp)

#'\code{plotDiag} plots a series of diagnostic plots
#'
#'@param gapsRes list returned by gapsRun, gapsMapRun, or CoGAPS
#'@export


plotDiag <-function(gapsRes)  {


  AMean <- gapsRes$Amean
  PMean <- gapsRes$Pmean
  ASD <- gapsRes$Asd
  PSD <- gapsRes$Psd
  ChiSq <- gapsRes$chiSqValues
  AtomsAEquil <- gapsRes$atomsAEquil
  AtomsASamp <- gapsRes$atomsASamp
  AtomsPEquil <- gapsRes$atomsPEquil
  AtomsPSamp <- gapsRes$atomsPSamp



  par(ask=TRUE)


  nbreak=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,5,10,15,20,25)
  mbreak=c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,5,10,15,20,25)
  jbreak=c(0.0001,0.0005,0.0008,0.001,0.005,0.008,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,5,10,15,20,25)



  heatmap.2(AMean, breaks=nbreak, Rowv = FALSE, Colv = FALSE,dendrogram="none",
            trace="none", density.info="none",main="Heatmap of A Matrix")

  hist(AMean, breaks=50, main="Histogram of A Matrix")


  heatmap.2(PMean, breaks=mbreak, Rowv = FALSE, Colv = FALSE,dendrogram="none",
            trace="none", density.info="none",main="Heatmap of P Matrix")

  hist(PMean, main="Histogram of P Matrix")


  heatmap.2(ASD, Rowv = FALSE,breaks=jbreak, Colv = FALSE,dendrogram="none",
            trace="none", density.info="none",main="A Standard Deviation Matrix")


  heatmap.2(PSD, Rowv = FALSE, breaks=jbreak, Colv = FALSE,dendrogram="none",
            trace="none", density.info="none",main="P Standard Deviation Matrix")


  plot(ChiSq, main="Chi Squared Values")


  par(mfrow=c(2, 2))
  plot(AtomsAEquil, main="Atoms A Equilibrium")
  plot(AtomsASamp, main="Atoms A Sample")
  plot(AtomsPEquil, main="Atoms P Equilibrium")
  plot(AtomsPSamp, main="Atoms P Sample")


}

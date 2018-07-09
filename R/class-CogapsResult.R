<<<<<<< HEAD
#' CogapsResult
#' @export
#' @example 
#' # Return output from CogapsResult
#' CogapsResult <- CoGAPS(GIST.D, CogapsParams)
#' @description Contains all output from Cogaps run
setClass("CogapsResult", slots=c(
    Amean = "matrix",
    Asd = "matrix",
    Pmean = "matrix",
    Psd = "matrix",
    seed = "numeric",
    meanChiSq = "numeric",
    diagnostics = "list"
))

#' Constructor for CogapsResult
#' @return initialized CogapsResult object
#' @importFrom methods callNextMethod
setMethod("initialize", "CogapsResult",
function(.Object, ...)
{
    .Object <- callNextMethod(.Object, ...)
    .Object
})

setMethod("show", signature("CogapsResult"),
function(object)
{
    nGenes <- nrow(object@Amean)
    nPatterns <- ncol(object@Amean)
    nSamples <- ncol(object@Pmean)

    print(paste("CogapsResult object with", nGenes, "genes and", nSamples, "samples"))
    print(paste(nPatterns, "patterns were learned"))
})

setMethod('plot', signature('CoGAPSResult'),
function(object)
{
    colors <- rainbow(nrow(object@Pmean))
    xlimits <- c(0, ncol(object@Pmean) + 1)
    ylimits <- c(0, (max(object@Pmean) * 1.05))

    plot(NULL, xlim=xlimits, ylim=ylimits, ylab="Relative Amplitude")

    for (i in 1:nrow(object@Pmean))
    {
        lines(x=1:ncol(object@Pmean), y=object@Pmean[i,], col=colors[i])
        points(x=1:ncol(object@Pmean), y=object@Pmean[i,], col=colors[i], pch=i)
    }

    legend("bottom", paste("Pattern", 1:nrow(object@Pmean), sep = ""),
    pch = 1:nrow(object@Pmean), lty=1, cex=0.8, col=colors, bty="y", ncol=5)
})

setMethod("MergeResultsWithSCE", signature("CogapsResult", "SingleCellExperiment"),
function(result, SCE)
{
    SCE@reducedDims <- SimpleList(Amean=result@Amean, Pmean=result@Pmean)
    return(SCE)
})

=======
#' CogapsResult
#' @export
#'
#' @description Contains all output from Cogaps run
setClass("CogapsResult", slots=c(
    Amean = "matrix",
    Asd = "matrix",
    Pmean = "matrix",
    Psd = "matrix",
    seed = "numeric",
    meanChiSq = "numeric",
    diagnostics = "list"
))

#' Constructor for CogapsResult
#' @return initialized CogapsResult object
#' @importFrom methods callNextMethod
setMethod("initialize", "CogapsResult",
function(.Object, ...)
{
    .Object <- callNextMethod(.Object, ...)
    .Object
})

setMethod("show", signature("CogapsResult"),
function(object)
{
    nGenes <- nrow(object@Amean)
    nPatterns <- ncol(object@Amean)
    nSamples <- ncol(object@Pmean)

    print(paste("CogapsResult object with", nGenes, "genes and", nSamples, "samples"))
    print(paste(nPatterns, "patterns were learned"))
})

setMethod('plot', signature(x='CogapsResult'),
function(x)
{
    colors <- rainbow(nrow(object@Pmean))
    xlimits <- c(0, ncol(object@Pmean) + 1)
    ylimits <- c(0, (max(object@Pmean) * 1.05))

    plot(NULL, xlim=xlimits, ylim=ylimits, ylab="Relative Amplitude")

    for (i in 1:nrow(object@Pmean))
    {
        lines(x=1:ncol(object@Pmean), y=object@Pmean[i,], col=colors[i])
        points(x=1:ncol(object@Pmean), y=object@Pmean[i,], col=colors[i], pch=i)
    }

    legend("bottom", paste("Pattern", 1:nrow(object@Pmean), sep = ""),
    pch = 1:nrow(object@Pmean), lty=1, cex=0.8, col=colors, bty="y", ncol=5)
})

setMethod('plot.heatmap', signature(x='CogapsResult'),
function(x)
{
  heatmap(Amean, Rowv=NA, Colv=NA)
  if (outputPDF != "")
  {
    dev.off()
  }    
  
})

# Plots patterns matrix as a line plot
setMethod('plot.smoothPatterns', signature(x='CogapsResult'),
          function(x)
          {
            # set optional NULL variables
            if (all(is.null(x)))
            {
              x <- 1:ncol(P)
            }
            
            # specify the break points for the plot
            if (all(is.null(breaks)))
            {
              breaks <- c(0,(ncol(P)+1))
            }
            
            # make the breaks in a uniform format
            if (length(breaks)==1)
            {
              breaks <- as.numeric(unique(unlist(strsplit(sub("\\(","",sub("\\]","",
                                                                           levels(cut(x,breaks)))),split=","))))
            }
            
            # check that the style of breaks matches the number of groups
            if (length(breakStyle) == 1)
            {
              breakStyle <- rep(breakStyle, length(breaks) - 1)
            }
            else
            {
              if (length(breakStyle) != length(breaks) -1)
              {
                if (length(breakStyle) == length(breaks))
                {
                  breaks <- c(min(x)-1.*abs(min(x)),breaks)
                }
                else
                {
                  stop('number of plot boundaries must match number of breaks')
                }
              }
            }
            
            # check that dimensions agree
            if (ncol(P) != length(x))
            {
              stop('length of x coordinates must match number of samples')
            }
            
            # reorder samples according to the group in which they obtain their maximum
            if (orderP)
            {
              PMax <- apply(P,1,max)
              xMax <- seq(from=ncol(P)+1, length.out=nrow(P))
              xMax <- xMax[order(PMax,decreasing=TRUE)]
              xTmp <- x
              PTmp <- P
              for (iP in order(PMax,decreasing=TRUE))
              {
                if (length(xTmp) < 1)
                {
                  break
                }
                xMax[iP] <- xTmp[which.max(PTmp[iP,])]
                PTmp <- PTmp[,which(xTmp!=xMax[iP])]
                xTmp <- xTmp[which(xTmp!=xMax[iP])]
              }
              POrder <- order(xMax)
              # reorder the pattern matrix
              P <- P[POrder,]
            }
            
            # find which group each sample belongs to
            sampleGroups <- cut(x,breaks[c(TRUE,breakStyle)])
            
            # specify the structure of these plots based upon the breaks
            if (!add)
            {
              split.screen(c(nrow(P), length(which(breakStyle))))
            }
            else
            {
              split.screen(c(nrow(P), length(which(breakStyle))), erase=FALSE)
            }
            
            scr <- 1
            
            for (iP in 1:nrow(P))
            {
              for (k in levels(sampleGroups))
              {
                screen(scr)
                xBorders <- as.numeric(unlist(strsplit(sub("\\]","",sub("\\(","",k)),split=",")))
                softBreaks <- breaks[which(breaks > xBorders[1] & breaks < xBorders[2])]
                plotBreaksTmp <- sort(c(xBorders,softBreaks))
                
                idxTmp <- which(x >= plotBreaksTmp[1] & x<plotBreaksTmp[2])
                
                if (k==levels(sampleGroups)[1])
                {
                  ylab <- paste('Pattern',iP,sep='')
                }
                else
                {
                  ylab <- ''
                }
                
                if (add)
                {
                  plot(x[idxTmp], loess(P[iP,idxTmp]~x[idxTmp])$fit, col=lineCol, type='l', axes=FALSE, xlab='', ylab='',lwd=3,
                       ylim=c(0,max(P[iP,])),xlim=xBorders,...)
                }
                else
                {
                  plot(x[idxTmp], loess(P[iP,idxTmp]~x[idxTmp])$fit, col=lineCol, lty=2, lwd=3, type='l', xlab=k, ylab=ylab,
                       ylim=c(0,max(P[iP,])),xlim=xBorders,...)
                }
                
                if (length(softBreaks) > 0)
                {
                  for (b in 1:length(softBreaks))
                  {
                    if (!add)
                    {
                      abline(v=softBreaks[b],lty=2,col=lineCol)
                    }
                    idxTmp <- which(x >= plotBreaksTmp[b+1] & x < plotBreaksTmp[b+2])
                    lines(x[idxTmp], loess(P[iP,idxTmp]~x[idxTmp])$fit, col=lineCol, type='l', lty=2, lwd=3)
                  }
                }
                
                if (plotPTS)
                {
                  points(x[which(sampleGroups==k)], P[iP, which(sampleGroups==k)],col=pointCol,pch=19, ...)
                }
                scr <- scr + 1
              }
            }
            close.screen(all.screens=TRUE)
          })

setMethod('plot.diagnostic', signature(x='CogapsResult'),
          function(x)
          {
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
            })


#' @export
setGeneric("MergeResultsWithSCE", function(result, SCE)
    {standardGeneric("MergeResultsWithSCE")})

#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("MergeResultsWithSCE", signature("CogapsResult", "SingleCellExperiment"),
function(result, SCE)
{
    SCE@reducedDims <- SimpleList(Amean=result@Amean, Pmean=result@Pmean)
    return(SCE)
})

>>>>>>> 2cdb0256a0675cee0896dab037fe9fe505b0fffb

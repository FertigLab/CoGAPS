#' Plot Decomposed A and P Matrices
#'
#' @details plots the output A and P matrices as a
#' heatmap and line plot respectively
#' @param A the mean A matrix
#' @param P the mean P matrix
#' @param outputPDF optional root name for PDF output, if
#'  not specified, output goes to screen
#' @export
plotGAPS <- function(A, P, outputPDF=NULL)
{
    if (is.null(outputPDF))
        dev.new()
    else
        pdf(file=paste(outputPDF, "-Patterns", ".pdf", sep=""))

    arrayIdx <- 1:ncol(P)
    maxP <- max(P)
    nPatt <- dim(P)[1]
    matplot(arrayIdx, t(P), type='l', lwd=3, xlim=c(1,ncol(P)), ylim=c(0,maxP),
        xlab='Samples',ylab='Relative Strength',col=rainbow(nPatt))
    title(main='Inferred Patterns')
    legend("topright", paste("Pattern", 1:nPatt, sep = ""), pch = 1:nPatt,
        lty=1,cex=0.8,col=rainbow(nPatt),bty="y",ncol=5)

    if (is.null(outputPDF))
        dev.new()
    else
        pdf(file=paste(outputPDF, "-Amplitude", ".pdf", sep=""))
 
    heatmap(A, Rowv=NA, Colv=NA)

    if (is.null(outputPDF))
        dev.off()
}

#' Plot the P Matrix
#'
#' @details plots the P matrix in a line plot with error bars
#' @param Pmean matrix of mean values of P
#' @param Psd matrix of standard deviation values of P
#' @export
plotP <- function(Pmean, Psd)
{
    nFactor <- nrow(Pmean)
    nObs <- ncol(Pmean)
    colors <- rainbow(nFactor)
    ylimits <- c(0, (max(Pmean + Psd) * 1.05))

    plotCI(x=1:nObs, y=Pmean[1,], col=colors[1], uiw=Psd[1,],
        ylim=ylimits, type='l', ylab="Relative Amplitude")

    for (i in 2:nFactor)
    {
        points(1:nObs, Pmean[i,], col=colors[i], pch=i)
        plotCI(1:nObs, Pmean[i,], col=colors[i], uiw=Psd[i,], add=TRUE)
    }

    legend("bottom", paste("Pattern", 1:nFactor, sep = ""),
        pch = 1:nFactor, lty=1, cex=0.8, col=colors, bty="y", ncol=5)
}

#' Diagnostic Plots
#'
#' @details plots a series of diagnostic plots
#' @param gapsRes list returned by CoGAPS
#' @export
plotDiag <-function(gapsRes)
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

    nbreak <- c(seq(0.1, 1.0, 0.1), seq(5, 25, 5))
    mbreak <- c(seq(0.01, 0.1, 0.01), nbreak)
    jbreak <- c(0.0001, 0.0005, 0.0008, 0.001, 0.005, 0.008, mbreak)

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

# plotAtoms: function to plot Rows of P with error bars
# History: v1.0 - MFI 8/2014

# Inputs: output of gapsRun, gapsMapRun, or CoGAPS

# Output: screen plot atoms

#'\code{plotAtoms} a simple plot of the number of atoms
#'from one of the vectors returned with atom numbers
#'
#'@param gapsRes the list resulting from applying GAPS
#'@param type the atoms to plot, values are "sampA", "sampP" ,
#'"equilA", or "equilP" to plot sampling or equilibration teop
#'atom numbers
#'@export
plotAtoms<-function(gapsRes, type='sampA')
{
    if (type == 'sampA')       atoms <- gapsRes$atomsASamp
    else if (type == 'sampP')  atoms <- gapsRes$atomsPSamp
    else if (type == 'equilA') atoms <- gapsRes$atomsAEquil
    else                       atoms <- gapsRes$atomsPEquil

    plot(atoms, xlab='Sample Number', ylab='Number of Atoms',
        main='Number of Atoms During MCMC Sampling')
}



# plotSmoothPatterns: function to plot pattern matrix
# History: v1.0 EJF original CoGAPS

# Inputs: P - P matrix
#         x - optional variables
#         breaks, breakStyle - plotting methods
#         orderP - whether to order patterns
#         plotPTS - whether to plot points on lines
#         pointCol - color of points
#         lineCol - color of line

# Output: plot of P matrix rows as lines

#'\code{plotSmoothPatterns} plots the output A and P matrices as a
#' heatmap and line plot respectively
#'
#'@param P the mean A matrix
#'@param x optional variables
#'@param breaks breaks in plots
#'@param breakStyle style of breaks
#'@param orderP whether to order patterns
#'@param plotPTS whether to plot points on lines
#'@param pointCol color of points
#'@param lineCol color of line
#'@param add logical specifying if bars should be added to an already existing plot; defaults to `FALSE'.
#'@param ... arguments to be passed to/from other methods.  For the default method these can include further arguments (such as `axes', `asp' and `main') and graphical parameters (see `par') which are passed to `plot.window()', `title()' and `axis'.
#'@export

plotSmoothPatterns <- function(P, x=NULL, breaks=NULL, breakStyle=T, orderP=!all(is.null(x)), plotPTS=F, pointCol='black',
                               lineCol='grey', add=F, ...)
{

    # set optional NULL variables
    if (all(is.null(x))) {
        x <- 1:ncol(P)
    }

    # specify the break points for the plot
    if (all(is.null(breaks))) {
        breaks <- c(0,(ncol(P)+1))
    }

    # make the breaks in a uniform format
    if(length(breaks)==1) {
        breaks <- as.numeric(unique(unlist(strsplit(sub("\\(","",sub("\\]","",levels(cut(x,breaks)))),split=","))))
    }

    # check that the style of breaks matches the number of groups
    if (length(breakStyle) == 1) {
        breakStyle <- rep(breakStyle, length(breaks) - 1)
    } else {
        if (length(breakStyle) != length(breaks) -1) {
            if (length(breakStyle) == length(breaks)) {
                breaks <- c(min(x)-1.*abs(min(x)),breaks)
            } else {
                stop('CoGAPS: plotSmoothPatterns: number of plot boundaries must match number of breaks in the plot')
            }
        }
    }


    # check that dimensions agree
    if (ncol(P) != length(x)) {
        stop('CoGAPS: plotSmoothPatterns: length of x coordinates must match number of samples in the columns of the P matrix')
    }

    # If desired, reorder samples according to the group in which they obtain their maximum
    if (orderP) {
        PMax <- apply(P,1,max)
        xMax <- seq(from=ncol(P)+1,length.out=nrow(P))[order(PMax,decreasing=T)]
        xTmp <- x
        PTmp <- P
        for (iP in order(PMax,decreasing=T)) {
            if (length(xTmp) < 1) {
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
    sampleGroups <- cut(x,breaks[c(T,breakStyle)])

    # specify the structure of these plots based upon the breaks
    if (!add) {
        split.screen(c(nrow(P), length(which(breakStyle))))
    } else {
        split.screen(c(nrow(P), length(which(breakStyle))), erase=F)
    }

    scr <- 1

    for (iP in 1:nrow(P)) {
        for (k in levels(sampleGroups)) {

            screen(scr)

            xBorders <- as.numeric(unlist(strsplit(sub("\\]","",sub("\\(","",k)),split=",")))

            softBreaks <- breaks[which(breaks > xBorders[1] & breaks < xBorders[2])]

            plotBreaksTmp <- sort(c(xBorders,softBreaks))

            idxTmp <- which(x >= plotBreaksTmp[1] & x<plotBreaksTmp[2])

            if (k==levels(sampleGroups)[1]) {
                ylab <- paste('Pattern',iP,sep='')
            } else {
                ylab <- ''
            }

            if (add) {
                plot(x[idxTmp], loess(P[iP,idxTmp]~x[idxTmp])$fit, col=lineCol, type='l', axes=F, xlab='', ylab='',lwd=3,
                     ylim=c(0,max(P[iP,])),xlim=xBorders,...)
            } else {
                plot(x[idxTmp], loess(P[iP,idxTmp]~x[idxTmp])$fit, col=lineCol, lty=2, lwd=3, type='l', xlab=k, ylab=ylab,
                     ylim=c(0,max(P[iP,])),xlim=xBorders,...)
            }

            if (length(softBreaks) > 0) {
                for (b in 1:length(softBreaks)) {
                    if (!add) {
                        abline(v=softBreaks[b],lty=2,col=lineCol)

                    }
                    idxTmp <- which(x >= plotBreaksTmp[b+1] & x < plotBreaksTmp[b+2])
                    lines(x[idxTmp], loess(P[iP,idxTmp]~x[idxTmp])$fit, col=lineCol, type='l', lty=2, lwd=3)

                }
            }

            if (plotPTS) {
                points(x[which(sampleGroups==k)], P[iP, which(sampleGroups==k)],col=pointCol,pch=19, ...)
            }

            scr <- scr + 1

        }
    }
    close.screen(all.screens=T)
}

#' plotPatternMarkers
#'
#' @param data the dataset from which the patterns where generated
#' @param patternMarkers the list of genes generated from the patternMarkers function
#' @param patternPalette a vector indicating what color should be used for each pattern
#' @param sampleNames names of the samples to use for labeling 
#' @param samplePalette  a vector indicating what color should be used for each sample
#' @param colDenogram logical indicating whether to display sample denogram
#' @param heatmapCol pallelet giving color scheme for heatmap
#' @param scale character indicating if the values should be centered and scaled in either 
#' the row direction or the column direction, or none. The default is "row". 
#' @param ... additional graphical parameters to be passed to \code{heatmap.2}
#' @export
#' @return heatmap of the \code{data} values for the \code{patternMarkers}
#' @seealso  \code{\link{heatmap.2}}
#' @examples \dontrun{
#' plotPatternMarkers(data=p,patternMarkers=PatternMarkers,patternPalette=NA,sampleNames=pd$sample,
#' samplePalette=pd$color,colDenogram=TRUE,heatmapCol="bluered", scale='row')
#'}

plotPatternMarkers <- function(
    data=NA,# the dataset from which the patterns where generated
    patternMarkers=NA,# the list of genes generated from the patternMarkers function
    patternPalette=NA,# a vector indicating what color should be used for each pattern
    sampleNames=NA,#names of the samples to use for labeling
    samplePalette=NULL,# a vector indicating what color should be used for each sample
    colDenogram=TRUE,# logical indicating whether to
    heatmapCol="bluered",
    scale='row',
    ...){

if(is.null(samplePalette)){samplePalette<-rep("black",dim(data)[2])}

### coloring of genes
patternCols <- rep(NA, length(unlist(patternMarkers)))
names(patternCols) <- unlist(patternMarkers)
for (i in 1:length(patternMarkers)) {
    patternCols[patternMarkers[[i]]] <- patternPalette[i]
}

heatmap.2(as.matrix(data[unlist(patternMarkers),],...),
          symkey=F,  symm=F, 
          scale=scale,
          col=heatmapCol,
          distfun=function(x) as.dist((1-cor(t(x)))/2),
          hclustfun=function(x) hclust(x,method="average"),
          Rowv=F,Colv=colDenogram,trace='none',
          RowSideColors = as.character(patternCols[unlist(patternMarkers)]),
          labCol= sampleNames,
          cexCol=.8,
          ColSideColors = as.character(samplePalette),
          rowsep = cumsum(sapply(patternMarkers,length)))
}

#' Binary Heatmap for Standardized A Matrix
#'
#' @details creates a binarized heatmap of the A matrix
#'  in which the value is 1 if the value in Amean is greater than
#'  threshold * Asd and 0 otherwise
#' @param Amean the mean estimate for the A matrix
#' @param Asd the standard deviations on Amean
#' @param threshold the number of standard deviations above zero
#'  that an element of Amean must be to get a value of 1
#' @export
binaryA <- function(Amean, Asd, threshold=3)
{
    BinA_Map <- ifelse(Amean/Asd > threshold, 1, 0)
    colnames(BinA_Map) <- colnames(BinA_Map, do.NULL=F, prefix="Pattern ")
    rownames(BinA_Map) <- rep(" ", nrow(BinA_Map))

    heatmap.2(BinA_Map, Rowv=FALSE, Colv=FALSE, dendrogram="none",
        scale="none", col = brewer.pal(3, "Blues"), trace="none",
        density.info="none", cexCol=1.3, srtCol=45,
        lmat=rbind(c(0,3), c(2,1), c(0,4)),
        lwid=c(1,10), lhei=c(1,4,1.2),
        main="Heatmap of Standardized A Matrix")
    message(paste("(Threshold = ", threshold, ")"))
}

#residuals: function to create a heatmap for the residuals
# History: v1.0 - JS, KC 7/2014

# Inputs: AMean_Mat - mean of A matrix
#         PMean_Mat - mean of P matrix
#         D - input data matrix
#         S - input standard deviation matrix

# Output: screen plot of residuals

#'\code{residuals} calculate residuals and produce heatmap
#'
#'@param AMean_Mat matrix of mean values for A from GAPS
#'@param PMean_Mat matrix of mean values for P from GAPS
#'@param D original data matrix run through GAPS
#'@param S original standard deviation matrix run through GAPS
#'@export


residuals=function(AMean_Mat, PMean_Mat, D, S)
{

    M_Mean <- AMean_Mat%*%PMean_Mat
    Resid_M_Mean<-as.matrix((D - M_Mean)/S)
    colnames(Resid_M_Mean) <- colnames(D)
    rownames(Resid_M_Mean) <- rownames(D)

    scaledRdYlBu <- colorRampPalette(brewer.pal(9,"RdYlBu"))(100)
    heatmap.2(Resid_M_Mean, Rowv = FALSE, Colv = FALSE,dendrogram="none",
        scale="none",col = scaledRdYlBu, trace="none",density.info="none",
        cexCol=1.33,srtCol=45,lmat=rbind(c(0, 3),c(2,1),c(0,4) ),
        lwid=c(1,10),lhei=c(1, 4, 1.2 ), main="Heatmap of Residuals")
}

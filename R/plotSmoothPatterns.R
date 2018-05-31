#' Plot Smooth Patterns
#'
#' @details plots the output A and P matrices as a heatmap and a
#' line plot respectively
#' @param P the mean A matrix
#' @param x optional variables
#' @param breaks breaks in plots
#' @param breakStyle style of breaks
#' @param orderP whether to order patterns
#' @param plotPTS whether to plot points on lines
#' @param pointCol color of points
#' @param lineCol color of line
#' @param add logical specifying if bars should be added to an already existing
#' plot; defaults to `FALSE'.
#' @param ... arguments to be passed to/from other methods.  For the default
#' method these can include further arguments (such as `axes', `asp' and
#' `main') and graphical parameters (see `par') which are passed to
#' `plot.window()', `title()' and `axis'.
#' @return plot
plotSmoothPatterns <- function(P, x=NULL, breaks=NULL, breakStyle=TRUE,
orderP=!all(is.null(x)), plotPTS=FALSE, pointCol='black', lineCol='grey',
add=FALSE, ...)
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
}

# plotGAPS: function to plot the decomposed A and P matrices
# History: EJF - original CoGAPS 

# Inputs: A - A matrix
#         P - P matrix
#         outputPDF - string for PDF output, if null uses screen

# Output: PDF file if indicated, otherwise none

#'\code{plotGAPS} plots the output A and P matrices as a
#' heatmap and line plot respectively
#'
#'@param A the mean A matrix
#'@param P the mean P matrix
#'@param outputPDF optional root name for PDF output, if
#'not specified, output goes to screen
#'@export

plotGAPS <- function(A, P, outputPDF="") {
	if (outputPDF != "") {
	  pdf(file=paste(outputPDF,"-Patterns",".pdf",sep=""))
	} else {
        	  dev.new()
	}

	arrayIdx <- 1:ncol(P)
    maxP <- max(P)
    nPatt <- dim(P)[1]
	matplot(arrayIdx, t(P), type='l', lwd=3, xlim=c(1,ncol(P)), ylim=c(0,maxP),
        xlab='Samples',ylab='Relative Strength',col=rainbow(nPatt))
	title(main='Inferred Patterns')
    legend("topright", paste("Pattern", 1:nPatt, sep = ""), pch = 1:nPatt,
        lty=1,cex=0.8,col=rainbow(nPatt),bty="y",ncol=5)
  
	if (outputPDF == "") {
	  dev.new()
	} else {
	  dev.off()
	  pdf(file=paste(outputPDF,"-Amplitude",".pdf",sep=""))
	}

	heatmap(A, Rowv=NA, Colv=NA)
	
	if (outputPDF != "") {
		dev.off()
	}
  
}

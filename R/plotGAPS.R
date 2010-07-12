# plot the decomposed A and P matrices
plotGAPS <- function(A, P, outputPDF="") {
	if (outputPDF != "") {
	  pdf(file=paste(outputPDF,"_Patterns",".pdf",sep=""))
	} else {
	  dev.new()
	}

	arrayIdx <- 1:ncol(P)
	matplot(arrayIdx, t(P), type='l', lwd=10, xlim=c(1,ncol(P)), ylim=c(0,1))
	title(main='Inferred patterns')
  
	if (outputPDF == "") {
	  dev.new()
	} else {
	  dev.off()
	  pdf(file=paste(outputPDF,"_Amplitude",".pdf",sep=""))
	}

	heatmap(A, Rowv=NA, Colv=NA) 
	
	if (outputPDF != "") {
		dev.off()
	}
  
}

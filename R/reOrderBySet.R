#' reOrderBySet
#'
#' @description <restructures output of gapsRun into a list containing each sets solution for Amean, Pmean, and Asd>
#' @param AP output of gapsRun in parallel
#' @param nFactor number of patterns
#' @param nSets number of sets
#' @param match which matrix to use for downstream matching. default is P
#' @return a list containing the \code{nSets} sets solution for Amean under "A", Pmean under "P", and Asd under "Asd"
reOrderBySet <- function(AP, nFactor, nSets, match="P")
{
	if(match=="P")
	{
		P<-do.call(rbind,lapply(AP, function(x) t(x@sampleFactors)))
		rownames(P)<-paste(rep(1:nSets,each=nFactor),rep(1:nFactor,nSets),sep=".")
		A<-lapply(AP, function(x) x@featureLoadings)
		Asd<-lapply(AP, function(x) x@featureStdDev)
		names(A)=names(Asd)<-paste(rep("Set",nSets),rep(1:nSets),sep="")
		return(list("A"=A,"Asd"=Asd,"P"=P))
	}

	if(match=="A")
	{
		A<-do.call(cbind,lapply(AP, function(x) x@featureLoadings))
		colnames(A)<-paste(rep(1:nSets,each=nFactor),rep(1:nFactor,nSets),sep=".")
		P<-lapply(AP, function(x) t(x@sampleFactors))
		Asd<-lapply(AP, function(x) x@featureStdDev)
		names(P)=names(Asd)<-paste(rep("Set",nSets),rep(1:nSets),sep="")
		return(list("A"=A,"Asd"=Asd,"P"=P))
	} 
}

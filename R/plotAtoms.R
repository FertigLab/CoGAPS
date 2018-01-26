#' Plot Number of Atoms
#'
#' @details a simple plot of the number of atoms
#' from one of the vectors returned with atom numbers
#' @param gapsRes the list resulting from applying GAPS
#' @param type the atoms to plot, values are "sampA", "sampP" ,
#' "equilA", or "equilP" to plot sampling or equilibration teop
#' atom numbers
#' @return plot
#' @examples
#' # Load the sample data from CoGAPS
#' data(SimpSim)
#' # Run plotAtoms
#' plotAtoms(results,type="sampA")
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

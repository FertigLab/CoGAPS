#Calculates significant genes in each pattern according to certain threshold
#Returns the significant gene names as well as well as the means of these matrices and number of genes in each

gapsInterPattern <- function(Amean, Asd, sdThreshold = 3)
{
    #number of rows and cols of Asd
    numGenes = length(Asd[,1]);
    numCols = length(Asd[1,]);

    #Vector holding the number of each significant gene in each column
    sigGeneNums = data.frame();

    #Temp number of sig genes in the col
    sigCount = 0;

    #Number to catch the amount of columns without significant genes
    numEmptyCols = 0;

    #Keep an array of the significant gene counts
    significantGeneNums = c(0);

    #Names of the genes from the data matrix
    geneNames = names(Asd[,1]);

    #Names of the genes that are significant from the data matrix
    sigGeneNames = "0";

    #The numerator of our statistic
    dimensionStatNumerator = 0;

    #The denominator of our statistic
    dimensionStatDenominator = 0;

    #The value of our statistic
    dimensionStatistic = 0;

    #A matrix holding the values of our statistics
    dimensionStatisticMatrix = matrix(nrow = numCols, ncol = numCols);

    #The mean of the statistic matrix
    sbar = 0;

    #A list to return the significant genes in each col and the statistic matrix
    results = list(list());

    #Scan in the significant genes from each column of Asd
    #The columns of sigGeneNums hold the significant genes for each col of Asd
    for(i in 1:numCols)
    {
        sigCount = 0;
        for(j in 1:numGenes)
        {
            if((Amean[j,i] - (sdThreshold*Asd[j,i])) > 0)
            {
                sigCount = sigCount + 1;
                sigGeneNums[sigCount, i] = j;
            }
        }

        if(sigCount == 0)
        {
            sigGeneNums[1, i] = 0;
            numEmptyCols = numEmptyCols + 1;
        }

        #Save the number of sigGenes
        significantGeneNums[i] = sigCount;

        #Get the names and store them
        if(sigCount != 0)
        {
            for(k in 1:sigCount)
            {
                sigGeneNames[k] = geneNames[sigGeneNums[k, i]];
            }
            results[[1]][[i]] = sigGeneNames;
            sigGeneNames = "0";
        }
        else
        {
            results[[1]][[i]] = "None";
        }
    }

    if(any(significantGeneNums == 0))
    {
        zeroSigCols = which(significantGeneNums == 0);
        print("Warning: No Significant Genes in Pattern(s): ");

        for(z in 1:length(zeroSigCols))
        {
            print(zeroSigCols[z]);
        }
    }

    #Now that we have the significant genes want to see if these genes are significant in other columns
    #Fill across the row then down the column

    #This compares the significant genes in the specified by 'j' to the same genes in Asd in the column specified by 'i'
    for(i in 1:numCols)
    {
        for(j in 1:numCols)
        {
            #Grab the number of significant genes from the interested column
            sigCount = sum(sigGeneNums[,j] > 0, na.rm = TRUE);

            if(sigCount != 0)
            {
                dimensionStatDenominator = sigCount;
                dimensionStatNumerator = 0;

                #loop through the number of significant genes and compare to these genes in the 'ith' col of Asd.
                #Find the significance of these genes, Amean - threshold * Asd
                for(k in 1:sigCount)
                {
                    if((Amean[sigGeneNums[k,j],i] - (sdThreshold*Asd[sigGeneNums[k,j],i])) > 0)
                    {
                        dimensionStatNumerator = dimensionStatNumerator + 1;
                    }
                }

                dimensionStatistic = dimensionStatNumerator/dimensionStatDenominator;

                dimensionStatisticMatrix[i, j] = dimensionStatistic;
            }
            else
            {
                dimensionStatisticMatrix[i, j] = 0;
            }
        }
    }

    #Find mean of the matrices (excluding the diagonal elements)
    #we can subtract numCols as the diagonal elements are 1
    sbar = ((sum(dimensionStatisticMatrix) - (numCols - numEmptyCols))/(length(dimensionStatisticMatrix) - (numCols - numEmptyCols)));

    results[[2]] = significantGeneNums;
    results[[3]] = t(dimensionStatisticMatrix);
    results[[4]] = sbar;

    names(results) = c("SignificantGeneNames", "SignificantGeneTotals", "SeparationMatrix", "InterPatternValue");

    return(results);
}

#Calculates significant genes in each pattern according to certain threshold
#Returns the significant gene names as well as well as the correlation matrices between these genes and the means of these matrices

gapsIntraPattern <- function(Amean, Asd, DMatrix, sdThreshold = 3)
{
    #number of rows and cols of Asd
    numGenes = length(Asd[,1]);
    numCols = length(Asd[1,]);

    #number of samples in DMatrix
    numSamp = ncol(DMatrix);

    #Vector holding the number of each significant gene in each column
    sigGeneNums = data.frame();

    #Temp number of sig genes in the col
    sigCount = 0;

    #Keep an array of the significant gene counts
    significantGeneNums = c(0);

    #A matrix to hold the significant genes in D for the current pattern
    #The matrix just acts as a subset of D, just eliminates non relevant rows
    tempSubsetD = matrix();

    #A matrix holding the values of our correlation coefficients between genes for the current column
    tempGeneCorrMatrix = matrix();

    #A list to hold all the correlation matrices
    geneCorrMatrices = list();

    #A list to hold all the means
    geneCorrMatrMeans = list();

    #The mean of all the correlation matrices
    cbar = 0;

    #A list to return the means and the matrices
    results = list();

    #Scan in the significant genes from each column of Asd
    #The columns of sigGeneNums hold the significant genes for each col of Asd
    for(i in 1:numCols)
    {
        sigCount = 0;
        for(j in 1:numGenes)
        {
            if((Amean[j,i] - (sdThreshold*Asd[j,i])) > 0)
            {
                sigCount = sigCount + 1;
                sigGeneNums[sigCount, i] = j;
            }
        }

        if(sigCount == 0)
        {
            sigGeneNums[1, i] = 0;
        }

        #Save the number of sigGenes
        significantGeneNums[i] = sigCount;
    }

    #If a pattern has no significant genes this is clearly an error so return such
    if(any(significantGeneNums == 0))
    {
        zeroSigCols = which(significantGeneNums == 0);
        warning("Warning: No Significant Genes in Pattern(s): ");

        for(z in 1:length(zeroSigCols))
        {
            message(zeroSigCols[z]);
        }
    }


    #Now that we have the significant genes want to grab these from our original D matrix
    #and find the sigGene x sigGene correlation matrix and find its mean

    for(j in 1:numCols)
    {
        #Grab the number of significant genes from the interested column
        sigCount = sum(sigGeneNums[,j] > 0, na.rm = TRUE);

        if(sigCount != 0)
        {

            #loop through the number of significant genes and pull out the rows of D that represent these genes.
            #Then find the correlation between them with the built in R corr function
            tempSubsetD = matrix(nrow = sigCount, ncol = numSamp);
            for(k in 1:sigCount)
            {
                #Subset D based on significant Genes
                #need to transpose as it reads this in as column vector otherwise
                tempSubsetD[k,] = t(DMatrix[sigGeneNums[k,j], ]);
            }

            #Find the correlation between these genes in D
            #Need to transpose as it calculates correlations between the columns
            tempGeneCorrMatrix = cor(t(tempSubsetD));

            #Find the mean of this matrix
            tempGeneCorrMatrMean = mean(tempGeneCorrMatrix);

        }
        else
        {
            tempGeneCorrMatrix = 0;
            tempGeneCorrMatrMean = 0;
        }

        #Save these in the overall list
        geneCorrMatrices[[j]] = tempGeneCorrMatrix;
        geneCorrMatrMeans[[j]] = tempGeneCorrMatrMean;

    }

    #Return as an overall list of lists
    # We return Corr Matrices themselves, their means, and the means of the means (cbar)
    results[[1]] = geneCorrMatrices;
    results[[2]] = geneCorrMatrMeans;

    #Return as an overall list of lists
    for(i in 1:numCols)
    {
        cbar = cbar + results[[2]][[i]];
    }

    cbar = cbar/numCols;
    results[[3]] = cbar;

    names(results) = c("CorrelationMatrices", "CorrelationMatrixMeans", "IntraPatternValue");

    return(results);
}


#' patternMarkers
#'
#' @param Amatrix A matrix of genes by weights resulting from CoGAPS or other NMF decomposition
#' @param scaledPmatrix logical indicating whether the corresponding pattern matrix was fixed to have max 1 during decomposition
#' @param Pmatrix the corresponding Pmatrix (patterns X samples) for the provided Amatrix (genes x patterns). This must be supplied if scaledPmatrix is FALSE.
#' @param threshold # the type of threshold to be used. The default "all" will distribute genes into pattern with the lowest ranking. The "cut" thresholding by the first gene to have a lower ranking, i.e. better fit to, a pattern.
#' @param lp a vector of weights for each pattern to be used for finding markers. If NA markers for each pattern of the A matrix will be used.
#' @param full logical indicating whether to return the ranks of each gene for each pattern
#'
#' @return By default a non-overlapping list of genes associated with each \code{lp}. If \code{full=TRUE} a data.frame of
#' genes rankings with a column for each \code{lp} will also be returned.
#' @export
#'
#' @examples \dontrun{
#' patternMarkers(Amatrix=AP$Amean,scaledPmatrix=FALSE,Pmatrix=NA,threshold="All",full=TRUE)
#' }
#'
patternMarkers <- function(
    Amatrix=NA, #A matrix of genes by weights resulting from CoGAPS or other NMF decomposition
    scaledPmatrix=FALSE, # logical indicating whether the corresponding pattern matrix was fixed to have max 1 during decomposition
    Pmatrix=NA, #the corresponding Pmatrix (patterns X samples) for the provided Amatrix (genes x patterns). This must be supplied if scaledPmatrix is FALSE.
    threshold="all", # the type of threshold to be used. The default "All" will distribute genes into pattern with the highest ranking. 
   # The "cut" thresholding by the first gene to have a lower ranking, i.e. better fit to, a pattern. 
    lp=NA, # a vector of weights for each pattern to be used for finding markers. If NA markers for each pattern of the A matrix will be used.
    full=FALSE #logical indicating whether to return the ranks of each gene for each pattern.
){


if(scaledPmatrix==FALSE){
    if(!is.na(Pmatrix)){
      pscale <- apply(Pmatrix,1,max)   # rescale p's to have max 1
      Amatrix <- sweep(Amatrix, 2, pscale, FUN="*")   # rescale A in accordance with p's having max 1
  }
    else(warning("P values must be provided if not already scaled"))
  }
# find the A with the highest magnitude
Arowmax <- t(apply(Amatrix, 1, function(x) x/max(x)))
pmax<-apply(Amatrix, 1, max)
# determine which genes are most associated with each pattern
ssranks<-matrix(NA, nrow=nrow(Amatrix), ncol=ncol(Amatrix),dimnames=dimnames(Amatrix))#list()
ssgenes<-matrix(NA, nrow=nrow(Amatrix), ncol=ncol(Amatrix),dimnames=NULL)
nP=dim(Amatrix)[2]
if(!is.na(lp)){
    if(length(lp)!=dim(Amatrix)[2]){
        warning("lp length must equal the number of columns of the Amatrix")
    }
        sstat <- apply(Arowmax, 1, function(x) sqrt(t(x-lp)%*%(x-lp)))
        ssranks[order(sstat),i] <- 1:length(sstat)
        ssgenes[,i]<-names(sort(sstat,decreasing=FALSE))
} else {
    for(i in 1:nP){
        lp <- rep(0,dim(Amatrix)[2])
        lp[i] <- 1
        sstat <- apply(Arowmax, 1, function(x) sqrt(t(x-lp)%*%(x-lp)))
        ssranks[order(sstat),i] <- 1:length(sstat)
        ssgenes[,i]<-names(sort(sstat,decreasing=FALSE))
    }
}
if(threshold=="cut"){
        geneThresh <- sapply(1:nP,function(x) min(which(ssranks[ssgenes[,x],x] > apply(ssranks[ssgenes[,x],],1,min))))
        ssgenes.th <- sapply(1:nP,function(x) ssgenes[1:geneThresh[x],x])
        #geneThresh <- apply(sweep(ssranks,1,t(apply(ssranks, 1, min)),"-"),2,function(x) which(x==0))
        #ssgenes.th <- lapply(geneThresh,names)
} else if(threshold=="all"){
        pIndx<-apply(ssranks,1,which.min)
        gBYp <- lapply(sort(unique(pIndx)),function(x) names(pIndx[pIndx==x]))
        ssgenes.th <- lapply(1:nP, function(x) ssgenes[which(ssgenes[,x] %in% gBYp[[x]]),x])
} else {stop("Threshold arguement not viable option")}

if(full){return(list("PatternMarkers"=ssgenes.th,"PatternRanks"=ssranks))
} else{return("PatternMarkers"=ssgenes.th)}
}


#' patternMatch4Parallel
#'
#' @param Ptot a matrix containing the total by set estimates of Pmean output from \code{reOrderBySet}
#' @param nSets number of parallel sets used to generate \code{Ptot}
#' @param cnt  number of branches at which to cut dendrogram
#' @param minNS minimum of individual set contributions a cluster must contain
#' @param cluster.method the agglomeration method to be used for clustering
#' @param ignore.NA logical indicating whether or not to ignore NAs from potential over dimensionalization. Default is FALSE.
#' @param bySet logical indicating whether to return list of matched set solutions from \code{Ptot}
#' @param ... additional parameters for \code{agnes}
#'
#' @return a matrix of concensus patterns by samples. If \code{bySet=TRUE} then a list of the set contributions to each
#' concensus pattern is also returned.
#' @export
#' @seealso \code{\link{agnes}}
#'
#'

patternMatch4Parallel <- function(Ptot,
  nSets, #number of sets
  cnt, # number of branches at which to cut dendrogram
  minNS, # minimum of sets a cluster must contain
  cluster.method="complete",
  ignore.NA=FALSE, # ignore NAs from potential over dimensionalization
  bySet=FALSE,
  ...){

#### read in CoGAPS results
cdir <- getwd()
#if(!is.null(path)){setwd(path)}
if(!is.null(minNS)){minNS=nSets/2}

if(ignore.NA==FALSE){if(anyNA(Ptot)){warning('Non-sparse matrixes produced. Reducing the number of patterns asked for and rerun.')}}
if(ignore.NA==TRUE){Ptot<-Ptot[complete.cases(Ptot),]}

####################################################################
# corr dist
corr.dist=cor(t(Ptot))
corr.dist=1-corr.dist
# cluster
#library(cluster)
clust=agnes(x=corr.dist,diss=T,method=cluster.method)
cut=cutree(as.hclust(clust),k=cnt)
#save.image(file=paste("CoGAPS.",nP,"P.",nS,"Set.CorrClustCut",cnt,".RData"))
####################################################################
#drop n<4 and get weighted Avg
cls=sort(unique(cut))
cMNs=matrix(nrow=cnt,ncol=dim(Ptot)[2])
rownames(cMNs)=cls
colnames(cMNs)=colnames(Ptot)

RtoMeanPattern <- list()
PByClust <- list()
for(i in cls){
   if(is.null(dim(Ptot[cut == i, ]))==TRUE){
       cMNs[i,] <- Ptot[cut == i, ]
       RtoMeanPattern[[i]] <- rep(1,length(Ptot[cut == i, ]))
       PByClust[[i]] <- t(as.matrix(Ptot[cut == i, ]))
   }
  else{
  cMNs[i,]=colMeans(Ptot[cut==i,])
  PByClust[[i]] <- Ptot[cut==i,]
  nIN=sum(cut==i)
  RtoMeanPattern[[i]] <- sapply(1:nIN,function(j) {round(cor(x=Ptot[cut==i,][j,],y=cMNs[i,]),3)})
  }
}

#drop n<minNS 
PByClustDrop <- list()
RtoMPDrop <- list()
for(i in cls){
  if(is.null(dim(PByClust[[i]]))==TRUE){next}
  if(dim(PByClust[[i]])[1]<minNS){next}
  else{
    #indx <- which(RtoMeanPattern[[i]]>.7,arr.ind = TRUE)
    PByClustDrop <- append(PByClustDrop,list(PByClust[[i]]))
    RtoMPDrop <- append(RtoMPDrop,list(RtoMeanPattern[[i]]))
  }
}


### split by corr  (build in drop if below minNS)
PByCDS <- list()
RtoMPDS <- list()
for(j in 1:length(PByClustDrop)){
  if(is.null(dim(PByClustDrop[[j]]))==TRUE){
      next
      }
  if(dim(PByClustDrop[[j]])[1]<minNS+nSets){
    PByCDS <- append(PByCDS,PByClustDrop[j])
    RtoMPDS <- append(RtoMPDS,RtoMPDrop[j])
  }
  if(dim(PByClustDrop[[j]])[1]>=minNS+nSets){
    corr.distPBCD=cor(t(PByClustDrop[[j]]))
    corr.distPBCD=1-corr.distPBCD
    clustPBCD=agnes(x=corr.distPBCD,diss=T,method="complete")
    cutPBCD=cutree(as.hclust(clustPBCD),k=2)
    g1 <- PByClustDrop[[j]][cutPBCD==1,]
    PByCDS <- append(PByCDS,list(g1))
    RtoMPDS <- append(RtoMPDS,list(sapply(1:dim(g1)[1],function(z) round(cor(x=g1[z,],y=colMeans(PByClustDrop[[j]][cutPBCD==1,])),3))))
    g2 <- PByClustDrop[[j]][cutPBCD==2,]
    if (is.null(dim(g2)[1])==FALSE){
        PByCDS <- append(PByCDS,list(g2))
        RtoMPDS <- append(RtoMPDS,list(sapply(1:dim(g2)[1],function(z) round(cor(x=g2[z,],y=colMeans(PByClustDrop[[j]][cutPBCD==2,])),3))))
    }
  }
#print(j)
#print(str(PByCDS))
}

#weighted.mean(PByClustDrop[[1]],RtoMPDrop[[1]])
PByCDSWavg<- t(sapply(1:length(PByCDS),function(z) apply(PByCDS[[z]],2,function(x) weighted.mean(x,(RtoMPDS[[z]])^3))))
rownames(PByCDSWavg) <- lapply(1:length(PByCDS),function(x) paste("Pattern",x))

#save
#save(PByCDSWavg,file=paste("PAatternSummary.UnScaled.",fname,".CoGAPS.",nP,"P.",nS,"Set.CorrClustCut",cnt,".RData",sep=""))

#scale ps
Pmax <- apply(PByCDSWavg,1,max)
PByCDSWavgScaled <- t(sapply(1:dim(PByCDSWavg)[1],function(x) PByCDSWavg[x,]/Pmax[x]))
rownames(PByCDSWavgScaled) <- rownames(PByCDSWavg)

if(bySet){
# return by set and final
PBySet<-PByCDS
return(list("consenusPatterns"=PByCDSWavgScaled,"PBySet"=PBySet))
} else {return(PByCDSWavgScaled)}

}


####################################################################
####################################################################


#' PatternMatcher Shiny Ap
#'
#' @param PBySet list of matched set solutions for the Pmatrix from an NMF algorithm
#' @param out optional name for saving output
#' @param order optional vector indicating order of samples for plotting. Default is NULL.
#' @param sample.color optional vector of colors of same length as colnames. Default is NULL.
#'
#' @return either an index of selected sets' contributions or the editted \code{PBySet} object
#' @export
#'
#' @examples \dontrun{
#' patternMatcher(PBySet,out,order,sample.color)
#' }
#'
#'
patternMatcher<-function(PBySet=NULL,out=NULL,order=NULL, sample.color=NULL) {

runApp(list(
  ui = pageWithSidebar(
    # Application title
    headerPanel('NMF Pattern Matching'),
    # Side pannel with controls
    sidebarPanel(
      # to upload file
      fileInput('file1',
                'Choose .Rda File',
                accept=c('.RData','.Rda','R data object','.rda')
      ),
      #
      uiOutput("pickplot"),
      uiOutput("checkbs"),
      downloadButton('downloadData', 'Download'),
      actionButton("end", "Return")
    ),
    # Main panel with plots
    mainPanel(
      plotOutput('plot1')
    )
  ),

  server = function(input, output, session) {
    #load in the data
    df<-reactive({
      if(!is.null(PBySet)){
        df<-PBySet
        return(df)
      }
      inFile <- input$file1 # get the path to the input file on the server
      if (is.null(inFile)){return(NULL)}
      load(inFile$datapath) #load it
      df <- get(load(inFile$datapath))# get the name of the object that was loaded
      return(df)
    })

    # get data name
    datName<-reactive({
      if(!is.null(out)){
        datName<-paste(out,'.SelectedPatterns.Rda',sep="")
        return(datName)
      }
      inFile <- input$file1
      if (is.null(inFile) & is.null(out)){
        datName<-"SelectedPatterns.Rda"
        return(datName)
      }
      if (is.null(inFile)){return(NULL)}
      fn<-strsplit(inFile$name,"[.]")[[1]][1]
      datName<-paste(fn,'.SelectedPatterns.Rda',sep="")
      return(datName)
    })


    mdf=reactive({# use that to give options for subsetting, some formatting may need to be removed
      dfx=df()
      if (is.null(dfx)){return(NULL)}
      mdf=melt(dfx,stringsAsFactors=FALSE) # melt the elements of the list
      colnames(mdf)<-c("BySet","Samples","value","Patterns")
      mdf$BySet<-as.character(mdf$BySet) # change them to characters
      mdf$Samples<-as.character(mdf$Samples)
      mdf$value=as.numeric(mdf$value) #make sure value is numeric for plotting
      str(mdf)
      return(mdf)
    })


    output$pickplot <- renderUI({# menu to select which matrix to look at/edit
      if (is.null(df())){return(NULL)}
      mdf2=mdf()
      selectInput("whichplot", "Select the Pattern to Plot",choices=unique(mdf2$Patterns))
    })


    output$checkbs <- renderUI({# make the checkboxes for each row of each matrix
      if (is.null(df())){return(NULL)}
      mdf2=mdf()
      lapply(unique(mdf2$Patterns), function(i) {
        subss <- unique(mdf2$BySet[mdf2$Patterns==i]) # find the rows (after it has been melted)
        tmp=sprintf("input.whichplot ==  %g", i) # create the javascript code to make this a conditional panel
        conditionalPanel(
          condition = tmp,
          checkboxGroupInput(paste("subs",i,sep=""), i, choices=subss, selected=subss) # the actual checkboxes for each, subs1, subs2, subsn
        )
      })
    })


    output$plot1 <- renderPlot({#plot the data, subset to the desired columns
      # if there has not been an uploaded matrix yet, don't even try to make a plot
      if (is.null(df())){return(NULL)}
      if (is.null(input$whichplot)){return(NULL)}
      par(mar = c(5.1, 4.1, 0, 1))
      mdf2=mdf() # grab the melted data frame to use the ggplot2 plot
      x=input$whichplot # which matrix to show
      x=as.numeric(x)
      tmp=input[[paste("subs",x,sep="")]] # get the rows that have been selected
      mdf2x=mdf2[which(mdf2$BySet%in%tmp),]
      if (!is.null(order) & !is.null(sample.color)){
       ggplot(mdf2x, aes(x=Samples, y=value, col=BySet,group=BySet))+
        geom_line() + scale_x_discrete(limits=order) +
        theme(axis.text.x = element_text(angle=45,family="Helvetica-Narrow", hjust = 1,colour = sample.color))
      } else if(!is.null(sample.color) & is.null(order)) {
      ggplot(mdf2x, aes(x=Samples, y=value, col=BySet,group=BySet))+
          geom_line() +
          theme(axis.text.x = element_text(angle=45,family="Helvetica-Narrow", hjust = 1,colour = sample.color))
      } else if(!is.null(order) & is.null(sample.color) ) {
        ggplot(mdf2x, aes(x=Samples, y=value, col=BySet,group=BySet))+
          geom_line() + scale_x_discrete(limits=order) +
          theme(axis.text.x = element_text(angle=45,family="Helvetica-Narrow", hjust = 1))
      } else {
        ggplot(mdf2x, aes(x=Samples, y=value, col=BySet,group=BySet))+
          geom_line() +
          theme(axis.text.x = element_text(angle=45,family="Helvetica-Narrow", hjust = 1))
      }
      #pplot
      #browser()
    })

    # create and download the final result file
    output$downloadData <- downloadHandler(
      filename = datName(), # set the file name
      content = function(file) {
        PatternsSelect <- lapply(1:length(mdf()), function(i) {input[[paste("subs",i,sep="")]]})
        save(PatternsSelect, file=file) # generate the object to save
      }
    )
    #stop app and return to R
    observeEvent(input$end, {
      mdf2=mdf()
      PatternsSelect <- sapply(1:length(df()), function(i) {input[[paste("subs",i,sep="")]]})
      selectPBySet <- mdf2[which(mdf2$BySet%in%PatternsSelect),]
      stopApp(returnValue = selectPBySet)
    })



  }

)
)

}

# reorderByPatternMatch: function to identify corresponding rows
#                        between two P matrices
# History: EJF - original CoGAPS

# Inputs: P - matrix to be matched
#         matchTo - matrix to match P to

# Output: P with reordered rows to match matchTo

#'\code{reorderByPatternMatch} plots the output A and P matrices as a
#' heatmap and line plot respectively
#'
#'@param P matrix to be matched
#'@param matchTo matrix to match P to
#'@export

reorderByPatternMatch <- function(P, matchTo) {

    # check that P and the matchTo matrix have the same dimensions for valid matching
    if (nrow(matchTo) != nrow(P) | ncol(matchTo) != ncol(P)) {
        stop('CoGAPS: reorderByPatternMatch: dimensions of P and matchTo must agree')
    }

    # ensuring that rownames match for simplicty of matching process
    row.names(matchTo) <- row.names(P)

    # compute the correlation between each entry
    corP <- cor(t(matchTo),t(P))

    # initalize the new matrix
    pMatch <- rep(0, nrow(P))
    names(pMatch) <- row.names(P)

    # match patterns in order of correlation
    for (p in 1:(nrow(P)-1)) {
        ptemp <- which(corP==max(corP),arr.ind=T)
        pMatch[row.names(corP)[ptemp[1]]] <- colnames(corP)[ptemp[2]]
        if (length(corP) > 1) {
            corP <- corP[-ptemp[1],-ptemp[2]]
        }
    }
    pMatch[which(pMatch==0)] <- setdiff(names(pMatch), pMatch)

    # return matched patterns
    return(P[pMatch,])

}

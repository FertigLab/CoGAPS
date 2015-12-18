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
	
	#The mean of the current gene correlation matrix
	tempCorrMatrMean = 0;
	
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

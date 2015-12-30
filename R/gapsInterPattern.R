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

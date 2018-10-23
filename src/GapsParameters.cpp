#include "GapsParameters.h"

Archive& operator<<(Archive &ar, const GapsParameters &p)
{
    ar << p.seed << p.nGenes << p.nSamples << p.nPatterns << p.nIterations
        << p.alphaA << p.alphaP << p.maxGibbsMassA << p.maxGibbsMassP
        << p.singleCell << p.useSparseOptimization;
    return ar;
}

Archive& operator>>(Archive &ar, GapsParameters &p)
{
    ar >> p.seed >> p.nGenes >> p.nSamples >> p.nPatterns >> p.nIterations
        >> p.alphaA >> p.alphaP >> p.maxGibbsMassA >> p.maxGibbsMassP
        >> p.singleCell >> p.useSparseOptimization;
    return ar;
}
    
void GapsParameters::calculateDataDimensions(const std::string &file)
{
    FileParser fp(file);
    
    nGenes = transposeData ? fp.nCol() : fp.nRow();
    nSamples = transposeData ? fp.nRow() : fp.nCol();

    if (subsetData && subsetGenes)
    {
        nGenes = dataIndicesSubset.size();
    }
    
    if (subsetData && !subsetGenes)
    {
        nSamples = dataIndicesSubset.size();
    }
}

void GapsParameters::calculateDataDimensions(const Matrix &mat)
{
    nGenes = transposeData ? mat.nCol() : mat.nRow();
    nSamples = transposeData ? mat.nRow() : mat.nCol();

    if (subsetData && subsetGenes)
    {
        nGenes = dataIndicesSubset.size();
    }
    
    if (subsetData && !subsetGenes)
    {
        nSamples = dataIndicesSubset.size();
    }
}

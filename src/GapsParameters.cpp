#include "GapsParameters.h"

void GapsParameters::peekCheckpoint(const std::string &file)
{
    Archive ar(file, ARCHIVE_READ);
    ar >> nPatterns >> seed >> nIterations >> whichFixedMatrix;
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

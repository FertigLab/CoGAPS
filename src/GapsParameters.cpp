#include "GapsParameters.h"
#include "utils/Archive.h"
#include "utils/GapsPrint.h"
#include "file_parser/FileParser.h"

void GapsParameters::print() const
{
    gaps_printf("\n---- C++ Parameters ----\n\n");
    gaps_printf("transposeData: %s\n", transposeData ? "TRUE" : "FALSE");
    gaps_printf("nGenes: %d\n", nGenes);
    gaps_printf("nSamples: %d\n", nSamples);
    gaps_printf("nPatterns: %d\n", nPatterns);
    gaps_printf("nIterations: %d\n", nIterations);
    gaps_printf("seed: %d\n", seed);
    gaps_printf("\n");
    gaps_printf("maxThreads: %d\n", maxThreads);
    gaps_printf("printMessages: %s\n", printMessages ? "TRUE" : "FALSE");
    gaps_printf("outputFrequency: %d\n", outputFrequency);
    gaps_printf("snapshotFrequency: %d\n", snapshotFrequency);
    gaps_printf("\n");
    gaps_printf("useSparseOptimization: %s\n", useSparseOptimization ? "TRUE" : "FALSE");
    gaps_printf("asynchronousUpdates: %s\n", asynchronousUpdates ? "TRUE" : "FALSE");
    gaps_printf("takePumpSamples: %s\n", takePumpSamples ? "TRUE" : "FALSE");
    gaps_printf("\n");
    gaps_printf("runningDistributed: %s\n", runningDistributed ? "TRUE" : "FALSE");
    gaps_printf("printThreadUsage: %s\n", printThreadUsage ? "TRUE" : "FALSE");
    gaps_printf("workerID: %d\n", workerID);
    gaps_printf("\n");
    gaps_printf("alphaA: %f\n", alphaA);
    gaps_printf("alphaP: %f\n", alphaP);
    gaps_printf("maxGibbsMassA: %f\n", maxGibbsMassA);
    gaps_printf("maxGibbsMassP: %f\n", maxGibbsMassP);
    gaps_printf("\n");
    gaps_printf("useCheckPoint: %s\n", useCheckPoint ? "TRUE" : "FALSE");
    gaps_printf("checkpointInterval: %d\n", checkpointInterval);
    gaps_printf("checkpointFile: %s\n", checkpointFile.c_str());
    gaps_printf("checkpointOutFile: %s\n", checkpointOutFile.c_str());
    gaps_printf("\n");
    gaps_printf("subsetData: %s\n", subsetData ? "TRUE" : "FALSE");
    gaps_printf("subsetGenes: %s\n", subsetGenes ? "TRUE" : "FALSE");
    gaps_printf("dataIndicesSubset.size(): %lu\n", dataIndicesSubset.size());
    gaps_printf("\n");
    gaps_printf("useFixedPatterns: %s\n", useFixedPatterns ? "TRUE" : "FALSE");
    gaps_printf("whichMatrixFixed: %c\n", whichMatrixFixed);
    gaps_printf("fixedPatterns.nRow(): %d\n", fixedPatterns.nRow());
    gaps_printf("fixedPatterns.nCol(): %d\n", fixedPatterns.nCol());
    gaps_printf("\n------------------------\n\n");
}

GapsParameters::GapsParameters(){}

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

Archive& operator<<(Archive &ar, const GapsParameters &p)
{
    ar << p.seed << p.nGenes << p.nSamples << p.nPatterns << p.nIterations
        << p.alphaA << p.alphaP << p.maxGibbsMassA << p.maxGibbsMassP
        << p.useSparseOptimization << p.checkpointInterval;
    return ar;
}

Archive& operator>>(Archive &ar, GapsParameters &p)
{
    ar >> p.seed >> p.nGenes >> p.nSamples >> p.nPatterns >> p.nIterations
        >> p.alphaA >> p.alphaP >> p.maxGibbsMassA >> p.maxGibbsMassP
        >> p.useSparseOptimization >> p.checkpointInterval;
    return ar;
}

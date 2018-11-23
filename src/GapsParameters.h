#ifndef __COGAPS_GAPS_PARAMETERS_H__
#define __COGAPS_GAPS_PARAMETERS_H__

#include "data_structures/Matrix.h"
#include "file_parser/FileParser.h"

#include <string>
#include <vector>

struct GapsParameters
{
public:

    template <class DataType>
    explicit GapsParameters(const DataType &data, bool t_transposeData=false,
        bool t_subsetData=false, bool t_subsetGenes=false,
        const std::vector<unsigned> t_dataIndicesSubset=std::vector<unsigned>());

    Matrix fixedMatrix;

    std::vector<unsigned> dataIndicesSubset;

    std::string checkpointFile;
    std::string checkpointOutFile;

    uint32_t seed;

    unsigned nGenes;
    unsigned nSamples;
    unsigned nPatterns;
    unsigned nIterations;
    unsigned maxThreads;
    unsigned outputFrequency;
    unsigned checkpointInterval;

    float alphaA;
    float alphaP;
    float maxGibbsMassA;
    float maxGibbsMassP;

    bool useFixedMatrix;
    bool subsetData;
    bool useCheckPoint;
    bool transposeData;
    bool singleCell;
    bool printMessages;
    bool subsetGenes;
    bool printThreadUsage;
    bool useSparseOptimization;

    char whichFixedMatrix;
    unsigned workerID;
    bool runningDistributed;

private:

    void calculateDataDimensions(const std::string &file);
    void calculateDataDimensions(const Matrix &mat);
};

Archive& operator<<(Archive &ar, const GapsParameters &p);
Archive& operator>>(Archive &ar, GapsParameters &p);

template <class DataType>
GapsParameters::GapsParameters(const DataType &data, bool t_transposeData,
bool t_subsetData, bool t_subsetGenes,
const std::vector<unsigned> t_dataIndicesSubset)
    :
fixedMatrix(Matrix()),
dataIndicesSubset(t_dataIndicesSubset),
checkpointFile(std::string()),
checkpointOutFile("gaps_checkpoint.out"),
seed(0),
nGenes(0),
nSamples(0),
nPatterns(3),
nIterations(1000),
maxThreads(1),
outputFrequency(500),
checkpointInterval(250),
alphaA(0.01f),
alphaP(0.01f),
maxGibbsMassA(100.f),
maxGibbsMassP(100.f),
useFixedMatrix(false),
subsetData(t_subsetData),
useCheckPoint(false),
transposeData(t_transposeData),
singleCell(false),
printMessages(true),
subsetGenes(t_subsetGenes),
printThreadUsage(true),
useSparseOptimization(false),
whichFixedMatrix('N'),
workerID(1),
runningDistributed(false)
{
    calculateDataDimensions(data);
}

#endif
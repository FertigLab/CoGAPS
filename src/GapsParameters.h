#ifndef __COGAPS_GAPS_PARAMETERS_H__
#define __COGAPS_GAPS_PARAMETERS_H__

struct GapsParameters
{
    Matrix fixedMatrix;

    std::vector<unsigned> dataIndicesSubset;

    std::string checkpointFile;
    std::string checkpointOutFile;

    uint32_t seed;

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

    char whichFixedMatrix;

    GapsParameters() :
        checkpointOutFile("gaps_checkpoint.out"),
        seed(0),
        nIterations(1000),
        maxThreads(1),
        outputFrequency(500),
        checkpointInterval(250),
        alphaA(0.01f),
        alphaP(0.01f),
        maxGibbsMassA(100.f),
        maxGibbsMassP(100.f),
        useFixedMatrix(false),
        subsetData(false),
        useCheckpoint(false),
        transposeData(false),
        singleCell(false),
        printMessages(true),
        subsetGenes(false),
        printThreadUsage(true),
        whichFixedMatrix('N')
    {}

    void peekCheckpoint(const std::string &file)
    {
        Archive ar(file, ARCHIVE_READ);
        ar >> nPatterns >> seed >> nIterations >> whichFixedMatrix;
    }
};

#endif
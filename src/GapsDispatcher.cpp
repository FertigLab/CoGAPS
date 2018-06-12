#include "GapsDispatcher.h"

#include <string>

void GapsDispatcher::runOneCycle(unsigned k)
{

}

Rcpp::List GapsDispatcher::run()
{
    unsigned iterPerCycle = mMaxIterations / 10;
    unsigned whichCycle = 1;
    while (whichCycle <= 10)
    {
        runOneCycle(iterPerCycle);
        RowMatrix pMaster = syncPatterns();
        broadcastPatterns(pMaster, static_cast<float>(whichCycle) / 10.f);
        whichCycle++;
    }
    
    for (unsigned i = 0; i < mRunners.size(); ++i)
    {
        mRunners[i]->enableSampling();
        mRunners[i]->fixPatternMatrix();
    }

    runOneCycle(mMaxIterations);

    // stitch together matrices
}

void GapsDispatcher::useDefaultUncertainty()
{

}

void GapsDispatcher::setUncertainty(const std::string &pathToMatrix)
{

}

void GapsDispatcher::setUncertainty(const RowMatrix &S)
{

}


void GapsDispatcher::loadData(const RowMatrix &D)
{

}

void GapsDispatcher::loadData(const std::string &pathToData);   
{

}

    

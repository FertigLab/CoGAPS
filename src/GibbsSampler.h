#ifndef __COGAPS_GIBBS_SAMPLER_H__
#define __COGAPS_GIBBS_SAMPLER_H__

#include "GapsAssert.h"
#include "Archive.h"
#include "Matrix.h"
#include "Random.h"
#include "Algorithms.h"
#include "ProposalQueue.h"
#include "AtomicDomain.h"

#include <Rcpp.h>

// forward declarations needed for friend classes
class AmplitudeGibbsSampler;
class PatternGibbsSampler;
class GapsStatistics;

/************************** GIBBS SAMPLER INTERFACE **************************/

template <class T, class MatA, class MatB>
class GibbsSampler
{
private:

    friend T; // prevent incorrect inheritance - only T can construct

    GibbsSampler(const Rcpp::NumericMatrix &D, const Rcpp::NumericMatrix &S,
        unsigned nrow, unsigned ncol, float alpha);

protected:

    MatA mMatrix;
    MatB* mOtherMatrix;
    MatB mDMatrix;
    MatB mSMatrix;
    MatB mAPMatrix;

    ProposalQueue mQueue;
    AtomicDomain mDomain;

    float mLambda;
    float mMaxGibbsMass;
    float mAnnealingTemp;
    
    unsigned mNumRows;
    unsigned mNumCols;
    uint64_t mBinSize;

    T* impl();

    void processProposal(const AtomicProposal &prop);
    void birth(uint64_t pos);
    void death(uint64_t pos, float mass);
    void move(uint64_t src, float mass, uint64_t dest);
    void exchange(uint64_t p1, float mass1, uint64_t p2, float mass2);
    float gibbsMass(unsigned row, unsigned col);

public:

    void update(unsigned nSteps);
    void setAnnealingTemp(float temp);
    
    float chi2() const;
    float nAtoms() const;

    // serialization
    //friend Archive& operator<<(Archive &ar, GibbsSampler &sampler);
    //friend Archive& operator>>(Archive &ar, GibbsSampler &sampler);
};

class AmplitudeGibbsSampler : public GibbsSampler<AmplitudeGibbsSampler, ColMatrix, RowMatrix>
{
private:

    friend GibbsSampler;
    friend PatternGibbsSampler;

    unsigned getRow(uint64_t pos) const;
    unsigned getCol(uint64_t pos) const;
    bool canUseGibbs(unsigned row, unsigned col) const;
    bool canUseGibbs(unsigned r1, unsigned c1, unsigned r2, unsigned c2) const;
    void updateAPMatrix(unsigned row, unsigned col, float delta);

public:

    AmplitudeGibbsSampler(const Rcpp::NumericMatrix &D,
        const Rcpp::NumericMatrix &S, unsigned nFactor, float alpha=0.f,
        float maxGibbsMass=0.f);

    void sync(PatternGibbsSampler &sampler);
};

class PatternGibbsSampler : public GibbsSampler<PatternGibbsSampler, RowMatrix, ColMatrix>
{
private:

    friend GibbsSampler;
    friend AmplitudeGibbsSampler;

    unsigned getRow(uint64_t pos) const;
    unsigned getCol(uint64_t pos) const;
    bool canUseGibbs(unsigned row, unsigned col) const;
    bool canUseGibbs(unsigned r1, unsigned c1, unsigned r2, unsigned c2) const;
    void updateAPMatrix(unsigned row, unsigned col, float delta);

public:

    PatternGibbsSampler(const Rcpp::NumericMatrix &D,
        const Rcpp::NumericMatrix &S, unsigned nFactor, float alpha=0.f,
        float maxGibbsMass=0.f);

    void sync(AmplitudeGibbsSampler &sampler);
};

/******************* IMPLEMENTATION OF TEMPLATED FUNCTIONS *******************/

template <class T, class MatA, class MatB>
GibbsSampler<T, MatA, MatB>::GibbsSampler(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, unsigned nrow, unsigned ncol, float alpha)
    :
mMatrix(nrow, ncol), mOtherMatrix(NULL), mDMatrix(D), mSMatrix(S),
mAPMatrix(D.nrow(), D.ncol()), mQueue(nrow * ncol, alpha),
mAnnealingTemp(0.f), mNumRows(nrow), mNumCols(ncol)
{
    mBinSize = std::numeric_limits<uint64_t>::max() / static_cast<uint64_t>(mNumRows * mNumCols);
    uint64_t remain = std::numeric_limits<uint64_t>::max() % (mNumRows * mNumCols);
    mQueue.setDomainSize(std::numeric_limits<uint64_t>::max() - remain);
}

template <class T, class MatA, class MatB>
T* GibbsSampler<T, MatA, MatB>::impl()
{
    return static_cast<T*>(this);
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::update(unsigned nSteps)
{
    unsigned n = 0;
    while (n < nSteps)
    {
        /*
        GAPS_ASSERT(nSteps - (queue.size() + n) >= 0);
        mQueue.populate(mDomain, nSteps - (mQueue.size() + n))

        // would making this a mulitple of nCores be better?
        unsigned nJobs = mQueue.size();
        for (unsigned i = 0; i < nJobs; ++i) // can be run in parallel
        {
            processProposal(mQueue[i]);
        }
        mQueue.clear();
        n += nJobs;
        GAPS_ASSERT(n <= nSteps);
        */
        mQueue.populate(mDomain, 1);
        GAPS_ASSERT(mQueue.size() == 1);
        processProposal(mQueue[0]);
        mQueue.clear(1);
        GAPS_ASSERT(mQueue.size() == 0);
        n++;
    }
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::processProposal(const AtomicProposal &prop)
{
    GAPS_ASSERT(prop.type == 'B' || prop.type == 'D' || prop.type == 'M' || prop.type == 'E');
    switch (prop.type)
    {
        case 'B': birth(prop.pos1); break;
        case 'D': death(prop.pos1, prop.mass1); break;
        case 'M': move(prop.pos1, prop.mass1, prop.pos2); break;
        case 'E': exchange(prop.pos1, prop.mass1, prop.pos2, prop.mass2); break;
    }
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::birth(uint64_t pos)
{
    unsigned row = impl()->getRow(pos);
    unsigned col = impl()->getCol(pos);
    float mass = impl()->canUseGibbs(row, col) ? gibbsMass(row, col)
        : gaps::random::exponential(mLambda);

    mDomain.insert(pos, mass);
    mMatrix(row, col) += mass;
    //impl()->updateAPMatrix(row, col, mass);
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::death(uint64_t pos, float mass)
{
    mQueue.rejectDeath();
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::move(uint64_t src, float mass, uint64_t dest)
{
    unsigned r1 = impl()->getRow(src);
    unsigned c1 = impl()->getCol(src);
    unsigned r2 = impl()->getRow(dest);
    unsigned c2 = impl()->getCol(dest);
    if (r1 == r2 && c1 == c2)
    {
        mDomain.erase(src);
        mDomain.insert(dest, mass);
    }
    else
    {
/*
        if (deltaLL * mAnnealingTemp >= std::log(gaps::random::uniform()))
        {
            mDomain.deleteAtom(p1);
            mDomain.addAtom(p2, mass);
            mMatrix(r1, c1) += -mass;
            mMatrix(r2, c2) += mass;
            impl()->updateAPMatrix(r1, c1, -mass);
            impl()->updateAPMatrix(r2, c2, mass);
        }
*/
    }
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::exchange(uint64_t p1, float mass1, uint64_t p2, float mass2)
{
    mQueue.rejectDeath();
}

template <class T, class MatA, class MatB>
float GibbsSampler<T, MatA, MatB>::gibbsMass(unsigned row, unsigned col)
{        
    //AlphaParameters alpha = impl()->alphaParameters(row, col);
    AlphaParameters alpha(10.f, 10.f);
    alpha.s *= mAnnealingTemp / 2.f;
    alpha.su *= mAnnealingTemp / 2.f;
    float mean  = (2.f * alpha.su - mLambda) / (2.f * alpha.s);
    float sd = 1.f / std::sqrt(2.f * alpha.s);

    float plower = gaps::random::p_norm(0.f, mean, sd);

    //float newMass = death ? mass : 0.f;
    float newMass = 0.f;
    if (plower < 1.f && alpha.s > 0.00001f)
    {
        newMass = gaps::random::inverseNormSample(plower, 1.f, mean, sd);
    }
    return std::max(0.f, std::min(newMass, mMaxGibbsMass));
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::setAnnealingTemp(float temp)
{
    mAnnealingTemp = temp;
}
  
template <class T, class MatA, class MatB>
float GibbsSampler<T, MatA, MatB>::chi2() const
{   
    return 2.f * gaps::algo::loglikelihood(mDMatrix, mSMatrix, mAPMatrix);
}
  
template <class T, class MatA, class MatB>
float GibbsSampler<T, MatA, MatB>::nAtoms() const
{   
    return mDomain.size();
}


//

//

//
////  template <class T, class MatA, class MatB>
////  void GibbsSampler<T, MatA, MatB>::processProposal(const AtomicProposal &prop)
////  {
////      GAPS_ASSERT(prop.type == 'B' || prop.type == 'D' || prop.type == 'M' || prop.type == 'E');
////      switch (prop.type)
////      {
////          case 'B': birth(prop.pos1); break;
////          //case 'D': death(prop.pos1, prop.mass1); break;
////          //case 'M': move(prop.pos1, prop.mass1, prop.pos2); break;
////          //case 'E': exchange(prop.pos1, prop.mass1, prop.pos2, prop.mass2); break;
////      }
////  }
////  
//
////  
////  template <class T, class MatA, class MatB>
////  void GibbsSampler<T, MatA, MatB>::death(uint64_t pos, float mass)
////  {
////      /*unsigned row = impl()->getRow(pos);
////      unsigned col = impl()->getCol(pos);
////      mMatrix(row, col) += -mass;
////      impl()->updateAPMatrix(row, col, -mass);
////      mDomain.erase(pos);
////      mQueue.acceptDeath();
////  
////      float newMass = impl()->canUseGibbs(row, col) ? gibbsMass(row, col)
////      : mass;
//      float deltaLL = impl()->computeDeltaLL(row, col, newMass);
//  
//      if (deltaLL * mAnnealingTemp >= std::log(gaps::random::uniform()))
//      {
//      float delta = mDomain.updateAtomMass(pos, newMass - mass);
//      mMatrix(row, col) += delta;
//      impl()->updateAPMatrix(row, col, delta);
//      if (mDomain.at(pos))
//      {
//          mQueue.rejectDeath();
//      }
//      }
//      else
//      {
//      mDomain.deleteAtom(pos);
//      mQueue.maxAtoms--;
//      }*/
//  }
//  
//  template <class T, class MatA, class MatB>
//  void GibbsSampler<T, MatA, MatB>::move(uint64_t src, float mass, uint64_t dest)
//  {
//      /*
//      if (r1 == r2 && c1 == c2)
//      {
//      mDomain.deleteAtom(p1);
//      mDomain.addAtom(p2, mass);
//      }
//      else
//      {
//      if (deltaLL * mAnnealingTemp >= std::log(gaps::random::uniform()))
//      {
//          mDomain.deleteAtom(p1);
//          mDomain.addAtom(p2, mass);
//          mMatrix(r1, c1) += -mass;
//          mMatrix(r2, c2) += mass;
//          impl()->updateAPMatrix(r1, c1, -mass);
//          impl()->updateAPMatrix(r2, c2, mass);
//      }
//      }
//      */
//  }
//  
//  template <class T, class MatA, class MatB>
//  void GibbsSampler<T, MatA, MatB>::exchange(uint64_t p1, float mass1, uint64_t p2, float mass2)
//  {
//      /*
//      float newMass = gaps::random::inverseGammaSample(0.f, mass1 + mass2, 2.f, 1.f / mLambda);
//      float delta1 = mass1 > mass2 ? newMass - mass1 : mass2 - newMass;
//      float delta2 = mass2 > mass1 ? newMass - mass2 : mass1 - newMass;
//      unsigned r1 = impl()->getRow(p1);
//      unsigned c1 = impl()->getCol(p1);
//      unsigned r2 = impl()->getRow(p2);
//      unsigned c2 = impl()->getCol(p2);
//  
//      if (r1 == r2 && c1 == c2)
//      {
//      mDomain.updateAtomMass(p1, delta1);
//      mDomain.updateAtomMass(p2, delta2);
//      }
//      else
//      {
//      // all the exchange code
//      }
//      */
//      //mQueue.rejectDeath();
//  }
//  
//  // don't return mass less than epsilon
//  template <class T, class MatA, class MatB>
//  float GibbsSampler<T, MatA, MatB>::gibbsMass(unsigned row, unsigned col)
//  {        
//  /*
//      AlphaParameters alpha = impl()->alphaParameters(row, col);
//      alpha.s *= mAnnealingTemp / 2.f;
//      alpha.su *= mAnnealingTemp / 2.f;
//      float mean  = (2.f * alpha.su - mLambda) / (2.f * alpha.s);
//      float sd = 1.f / std::sqrt(2.f * alpha.s);
//  
//      float plower = gaps::random::p_norm(0.f, mean, sd);
//  
//      float newMass = death ? mass : 0.f;
//      if (plower < 1.f && alpha.s > 0.00001f)
//      {
//      newMass = gaps::random::inverseNormSample(plower, 1.f, mean, sd);
//      }
//      return std::max(0.f, std::min(newMax, mMaxGibbsMass));
//  */
//  }

#endif
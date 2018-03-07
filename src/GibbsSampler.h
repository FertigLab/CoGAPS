#ifndef __GAPS_GIBBS_SAMPLER_H__
#define __GAPS_GIBBS_SAMPLER_H__

#include "GapsAssert.h"
#include "GapsStatistics.h"
#include "Archive.h"
#include "Matrix.h"
#include "Random.h"
#include "Algorithms.h"
#include "ProposalQueue.h"
#include "AtomicDomain.h"

#include <Rcpp.h>

// stats should be friend class
// CRTP not really needed here, free friend functions are clear enough,
// it would reduce some repeated code - but only a few lines and probably
// not worth the headache
// maybe convert to CRTP once everything is working cleanly

class AmplitudeGibbsSampler;
class PatternGibbsSampler;
class GapsStatistics;

class AmplitudeGibbsSampler
{
private:

    friend PatternGibbsSampler;
    friend GapsStatistics;

    ColMatrix mAMatrix;
    RowMatrix *mPMatrix;
    RowMatrix mDMatrix;
    RowMatrix mSMatrix;
    RowMatrix mAPMatrix;

    ProposalQueue mQueue;
    AtomicDomain mDomain;

    float mLambda;
    float mMaxGibbsMass;
    float mAnnealingTemp;
    
    unsigned mNumRows;
    unsigned mNumCols;
    unsigned mBinSize;

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
    float nAtoms() const;
    void setAnnealingTemp(float temp);
    float chi2() const;

    template <class Sampler>
    friend void update(Sampler&, unsigned);

    template <class Sampler>
    friend void processProposal(Sampler &sampler, const AtomicProposal &prop);

    template <class Sampler>
    friend void birth(Sampler &sampler, uint64_t pos);

    template <class Sampler>
    friend void death(Sampler &sampler, uint64_t pos, float mass);

    template <class Sampler>
    friend void exchange(Sampler &sampler, uint64_t p1, float mass1, uint64_t p2, float mass2);
};

class PatternGibbsSampler
{
private:

    friend AmplitudeGibbsSampler;
    friend GapsStatistics;

    RowMatrix mPMatrix;
    ColMatrix *mAMatrix;
    ColMatrix mDMatrix;
    ColMatrix mSMatrix;
    ColMatrix mAPMatrix;

    ProposalQueue mQueue;
    AtomicDomain mDomain;

    float mLambda;
    float mMaxGibbsMass;
    float mAnnealingTemp;
    
    unsigned mNumRows;
    unsigned mNumCols;
    unsigned mBinSize;

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
    float nAtoms() const;
    void setAnnealingTemp(float temp);
    float chi2() const;

    template <class Sampler>
    friend void update(Sampler&, unsigned);

    template <class Sampler>
    friend void processProposal(Sampler &sampler, const AtomicProposal &prop);

    template <class Sampler>
    friend void birth(Sampler &sampler, uint64_t pos);

    template <class Sampler>
    friend void death(Sampler &sampler, uint64_t pos, float mass);

    template <class Sampler>
    friend void exchange(Sampler &sampler, uint64_t p1, float mass1, uint64_t p2, float mass2);
};

template <class Sampler>
void update(Sampler &sampler, unsigned nSteps)
{
    unsigned n = 0;
    while (n < nSteps)
    {
        /*
        assert(nSteps - (queue.size() + n) >= 0);
        mQueue.populate(mDomain, nSteps - (mQueue.size() + n))

        // would making this a mulitple of nCores be better?
        unsigned nJobs = mQueue.size();
        for (unsigned i = 0; i < nJobs; ++i) // can be run in parallel
        {
            processProposal(mDomain, mQueue[i]);
        }
        mQueue.clear();
        n += nJobs;
        assert(n <= nSteps);
        */
        sampler.mQueue.populate(sampler.mDomain, 1);
        GAPS_ASSERT(sampler.mQueue.size() == 1);
        processProposal(sampler, sampler.mQueue[0]);
        sampler.mQueue.clear();
        n++;
    }
}

template <class Sampler>
void processProposal(Sampler &sampler, const AtomicProposal &prop)
{
    GAPS_ASSERT(prop.type == 'B' || prop.type == 'D' || prop.type == 'M' || prop.type == 'E');
    switch (prop.type)
    {
        case 'B': birth(sampler, prop.pos1); break;
        case 'D': death(sampler, prop.pos1, prop.mass1); break;
        //case 'M': move(prop.pos1, prop.mass1, prop.pos2); break;
        case 'E': exchange(sampler, prop.pos1, prop.mass1, prop.pos2, prop.mass2); break;
    }
}

template <class Sampler>
void birth(Sampler &sampler, uint64_t pos)
{
    unsigned row = sampler.getRow(pos);
    unsigned col = sampler.getCol(pos);
    float mass = gaps::random::exponential(sampler.mLambda);
    sampler.mDomain.insert(pos, mass);
}

template <class Sampler>
void death(Sampler &sampler, uint64_t pos, float mass)
{
    sampler.mQueue.rejectDeath();
}

template <class Sampler>
void exchange(Sampler &sampler, uint64_t p1, float mass1, uint64_t p2, float mass2)
{
    sampler.mQueue.rejectDeath();
}

//  template <class T, class MatA, class MatB>
//  void GibbsSampler<T, MatA, MatB>::processProposal(const AtomicProposal &prop)
//  {
//      GAPS_ASSERT(prop.type == 'B' || prop.type == 'D' || prop.type == 'M' || prop.type == 'E');
//      switch (prop.type)
//      {
//          case 'B': birth(prop.pos1); break;
//          //case 'D': death(prop.pos1, prop.mass1); break;
//          //case 'M': move(prop.pos1, prop.mass1, prop.pos2); break;
//          //case 'E': exchange(prop.pos1, prop.mass1, prop.pos2, prop.mass2); break;
//      }
//  }
//  
//  template <class T, class MatA, class MatB>
//  void GibbsSampler<T, MatA, MatB>::birth(uint64_t pos)
//  {
//      //GAPS_ASSERT(impl());
//      //unsigned row = impl()->getRow(pos);
//      //unsigned col = impl()->getCol(pos);
//      //float mass = impl()->canUseGibbs(row, col) ? gibbsMass(row, col)
//      //    : gaps::random::exponential(mLambda);
//      /*float mass = gaps::random::exponential(mLambda);
//  
//      mDomain.insert(pos, mass);
//      mMatrix(row, col) += mass;
//      impl()->updateAPMatrix(row, col, mass);*/
//  }
//  
//  template <class T, class MatA, class MatB>
//  void GibbsSampler<T, MatA, MatB>::death(uint64_t pos, float mass)
//  {
//      /*unsigned row = impl()->getRow(pos);
//      unsigned col = impl()->getCol(pos);
//      mMatrix(row, col) += -mass;
//      impl()->updateAPMatrix(row, col, -mass);
//      mDomain.erase(pos);
//      mQueue.acceptDeath();
//  
//      float newMass = impl()->canUseGibbs(row, col) ? gibbsMass(row, col)
//      : mass;
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
//  
//  template <class T, class MatA, class MatB>
//  void GibbsSampler<T, MatA, MatB>::syncAP(const MatA &otherAP)
//  {   
//      mAPMatrix = otherAP;
//  }
//  
//  template <class T, class MatA, class MatB>
//  const MatB& GibbsSampler<T, MatA, MatB>::APMatrix() const
//  {   
//      return mAPMatrix;
//  }
//  
//  template <class T, class MatA, class MatB>
//  void GibbsSampler<T, MatA, MatB>::setAnnealingTemp(float temp)
//  {
//      mAnnealingTemp = temp;
//  }
//  
//  template <class Sampler>
//  float chi2(const Sampler &sampler)
//  {   
//      //return 2.f * gaps::algo::loglikelihood(mDMatrix, mSMatrix, mAPMatrix);
//      return 0.f;
//  }
//  
//  template <class T, class MatA, class MatB>
//  float GibbsSampler<T, MatA, MatB>::nAtoms() const
//  {   
//      //return 2.f * gaps::algo::loglikelihood(mDMatrix, mSMatrix, mAPMatrix);
//      return 0.f;
//  }

#endif
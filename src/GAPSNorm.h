#ifndef __COGAPS_GAPSNORM_H__
#define __COGAPS_GAPSNORM_H__

#include "Matrix.h"

#include <iostream>
#include <vector>
#include <utility>

namespace GAPSNorm
{

    void computeMock(Matrix &M, const Matrix &A, const Matrix &P,
                     unsigned int nFactor);

    double calChi2(const Matrix &D, const Matrix &S, const Matrix &A,
                   const Matrix &P, unsigned int nFactor);

    double calcDeltaLL1E(char matrix_label, const Matrix &D, const Matrix &S,
                         const Matrix &A, const Matrix &P,
                         const std::vector<ElementChange> ElemChange,
                         unsigned int nFactor);

    double calcDeltaLL2E(char matrix_label, const Matrix &D, const Matrix &S,
                         const Matrix &A, const Matrix &P, 
                         const std::vector<ElementChange> ElemChange,
                         unsigned int nFactor);
                                

    double calcDeltaLLGen(char matrix_label, const Matrix &D, const Matrix &S,
                          const Matrix &A, const Matrix &P,
                          const std::vector<ElementChange> ElemChange,
                          unsigned int nFactor);

    // For use when there are fixed patterns and an entire row is to be changed.
    double calcDeltaLLMap(char matrix_label, const Matrix &D, const Matrix &S,
                          const Matrix &A, const Matrix &P,
                          std::vector<double> &newPat, unsigned int chPat,
                          unsigned int nFactor);

    // For move exchange / multiple pattern changes
    double calcDeltaLL2Map(char matrix_label, const Matrix &D, const Matrix &S,
                           const Matrix &A, const Matrix &P,
                           std::vector<double> &newPat1, unsigned int chPat1,
                           std::vector <double> &newPat2, unsigned int chPat2,
                          unsigned int nFactor);

    /** NEW METHOD
    *   @short Calculate the parameters involved in exchange action for further use
    *   @return a pair <s, su> of alpha distribution parameters for further calculations.
    */
    std::pair<double, double> calcAlphaParameters(char the_matrix_label, unsigned int nFactor,
        const Matrix &D, const Matrix &S, const Matrix &AOrig, const Matrix &POrig,
        unsigned int iGene1, unsigned int iPattern1,
        unsigned int iGene2, unsigned int iPattern2,
        unsigned int iSample1, unsigned int iSample2);
}

#endif

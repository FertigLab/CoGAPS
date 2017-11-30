#ifndef __COGAPS_GAPSNORM_H__
#define __COGAPS_GAPSNORM_H__

#include <iostream>
#include <vector>
#include <boost/tuple/tuple.hpp>

typedef boost::tuple<unsigned int, unsigned int, double> ElementChange;

namespace GAPSNorm
{

    void computeMock(Matrix &M, const Matrix &A, const Matrix &P,
                     unsigned int nFactor);

    double calChi2(const Matrix &D, const Matrix &S, const Matrix &A,
                   const Matrix &P, unsigned int nFactor);

    double calcDeltaLL1E(char matrix_label, const Matrix &D, const Matrix &S,
                         const Matrix &A, const Matrix &P, unsigned int nFactor,
                         const std::vector<ElementChange> ElemChange);

    double calcDeltaLL2E(char matrix_label, const Matrix &D, const Matrix &S,
                         const Matrix &A, const Matrix &P, unsigned int nFactor,
                         const std::vector<ElementChange> ElemChange);
                                

    double calcDeltaLLGen(char matrix_label, const Matrix &D, const Matrix &S,
                          const Matrix &A, const Matrix &P, unsigned int nFactor,
                          const std::vector<ElementChange> ElemChange);

    // For use when there are fixed patterns and an entire row is to be changed.
    double calcDeltaLLMap(char matrix_label, const Matrix &D, const Matrix &S,
                          const Matrix &A, const Matrix &P, unsigned int nFactor,
                          std::vector<double> &newPat, unsigned int chPat);

    // For move exchange / multiple pattern changes
    double calcDeltaLL2Map(char matrix_label, const Matrix &D, const Matrix &S,
                           const Matrix &A, const Matrix &P, unsigned int nFactor,
                           std::vector<double> &newPat1, unsigned int chPat1,
                           std::vector <double> &newPat2, unsigned int chPat2)

    /** NEW METHOD
    *   @short Calculate the parameters involved in exchange action for further use
    *   @return a pair <s, su> of alpha distribution parameters for further calculations.
    */
    std::pair<double, double> calcAlphaParameters(char the_matrix_label, unsigned int nFactor,
        const Matrix &D, const Matrix &S, double **AOrig, double **POrig, unsigned int iGene1,
        unsigned int iPattern1, unsigned int iGene2, unsigned int iPattern2, unsigned int iSample1,
        unsigned int iSample2);
}

#endif

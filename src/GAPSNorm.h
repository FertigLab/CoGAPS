#ifndef __COGAPS_GAPSNORM_H__
#define __COGAPS_GAPSNORM_H__

#include "Matrix.h"

#include <iostream>
#include <vector>
#include <utility>

namespace gaps
{

namespace norm
{
    double calChi2(const Matrix &D, const Matrix &S, const Matrix &A,
        const Matrix &P, unsigned int nFactor);

    double calcDeltaLL1E(const Matrix &D, const Matrix &S, const Matrix &AP,
        const MatrixChange &change);

    double calcDeltaLL2E(const Matrix &D, const Matrix &S, const Matrix &AP,
        const MatrixChange &change);
                                
    /** NEW METHOD
    *   @short Calculate the parameters involved in exchange action for further use
    *   @return a pair <s, su> of alpha distribution parameters for further calculations.
    */
    /*std::pair<double, double> calcAlphaParameters(char the_matrix_label, unsigned int nFactor,
        const Matrix &D, const Matrix &S, const Matrix &AOrig, const Matrix &POrig,
        unsigned int iGene1, unsigned int iPattern1,
        unsigned int iGene2, unsigned int iPattern2,
        unsigned int iSample1, unsigned int iSample2);*/

} // namespace norm

} // namespace gaps

#endif

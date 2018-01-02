#include "Algorithms.h"
#include "Matrix.h"

#include <algorithm>

#define GAPS_SQ(x) ((x) * (x))

matrix_data_t gaps::algo::mean(const TwoWayMatrix &mat)
{
    return gaps::algo::sum(mat) / (mat.nRow() * mat.nCol());
}

matrix_data_t gaps::algo::sum(const Vector &vec)
{
    matrix_data_t sum = 0.0;
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        sum += vec(i);
    }
    return sum;
}

matrix_data_t gaps::algo::sum(const TwoWayMatrix &mat)
{
    matrix_data_t sum = 0.0;
    for (unsigned r = 0; r < mat.nRow(); ++r)
    {
        sum += gaps::algo::sum(mat.getRow(r));
    }
    return sum;
}

matrix_data_t gaps::algo::sum(const RowMatrix &mat)
{
    matrix_data_t sum = 0.0;
    for (unsigned r = 0; r < mat.nRow(); ++r)
    {
        sum += gaps::algo::sum(mat.getRow(r));
    }
    return sum;
}

matrix_data_t gaps::algo::sum(const ColMatrix &mat)
{
    matrix_data_t sum = 0.0;
    for (unsigned c = 0; c < mat.nCol(); ++c)
    {
        sum += gaps::algo::sum(mat.getCol(c));
    }
    return sum;
}

matrix_data_t gaps::algo::nonZeroMean(const TwoWayMatrix &mat)
{
    matrix_data_t sum = 0.0;
    unsigned int nNonZeros = 0;
    for (unsigned int r = 0; r < mat.nRow(); ++r)
    {
        for (unsigned int c = 0; c < mat.nCol(); ++c)
        {
            if (mat(r,c) != 0.0)
            {
                nNonZeros++;
                sum += mat(r,c);
            }
        }
    }
    return sum / nNonZeros;
}

ColMatrix gaps::algo::computeStdDev(const ColMatrix &stdMat,
const ColMatrix &meanMat, unsigned nUpdates)
{
    ColMatrix retMat(stdMat.nRow(), stdMat.nCol());
    double meanTerm = 0.0;
    for (unsigned r = 0; r < retMat.nRow(); ++r)
    {
        for (unsigned c = 0; c < retMat.nCol(); ++c)
        {
            meanTerm = meanMat(r,c) * meanMat(r,c) / (double)nUpdates;
            retMat(r,c) = std::sqrt((stdMat(r,c) - meanTerm) / ((double)nUpdates - 1.0));
        }
    }
    return retMat;
}

RowMatrix gaps::algo::computeStdDev(const RowMatrix &stdMat,
const RowMatrix &meanMat, unsigned nUpdates)
{
    RowMatrix retMat(stdMat.nRow(), stdMat.nCol());
    double meanTerm = 0.0;
    for (unsigned r = 0; r < retMat.nRow(); ++r)
    {
        for (unsigned c = 0; c < retMat.nCol(); ++c)
        {
            meanTerm = meanMat(r,c) * meanMat(r,c) / (double)nUpdates;
            retMat(r,c) = std::sqrt((stdMat(r,c) - meanTerm) / ((double)nUpdates - 1.0));
        }
    }
    return retMat;
}

// horribly slow, don't call often
void gaps::algo::matrixMultiplication(TwoWayMatrix &C, const ColMatrix &A,
const RowMatrix &B)
{
    for (unsigned i = 0; i < C.nRow(); ++i)
    {
        for (unsigned j = 0; j < C.nCol(); ++j)
        {
            double sum = 0.0;
            for (unsigned k = 0; k < A.nCol(); ++k)
            {
                sum += A(i,k) * B(k,j);
            }
            C.set(i, j, sum);
        }
    }
}

Vector gaps::algo::scalarMultiple(const Vector &vec, double n)
{
    Vector retVec(vec.size());
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        retVec(i) = vec(i) * n;
    }
    return retVec;
}

Vector gaps::algo::squaredScalarMultiple(const Vector &vec, double n)
{
    Vector retVec(vec.size());
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        retVec(i) = GAPS_SQ(vec(i) * n);
    }
    return retVec;
}

Vector gaps::algo::squaredScalarDivision(const Vector &vec, double n)
{
    Vector retVec(vec.size());
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        retVec(i) = GAPS_SQ(vec(i) / n);
    }
    return retVec;
}

ColMatrix gaps::algo::scalarDivision(const ColMatrix &mat, double n)
{
    ColMatrix retMat(mat.nRow(), mat.nCol());
    for (unsigned r = 0; r < retMat.nRow(); ++r)
    {
        for (unsigned c = 0; c < retMat.nCol(); ++c)
        {
            retMat(r,c) = mat(r,c) / n;
        }
    }
    return retMat;
}

Vector gaps::algo::scalarDivision(const Vector &vec, matrix_data_t n)
{
    Vector temp(vec.size());
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        temp(i) = vec(i) / n;
    }
    return temp;
}

RowMatrix gaps::algo::scalarDivision(const RowMatrix &mat, double n)
{
    RowMatrix retMat(mat.nRow(), mat.nCol());
    for (unsigned r = 0; r < retMat.nRow(); ++r)
    {
        for (unsigned c = 0; c < retMat.nCol(); ++c)
        {
            retMat(r,c) = mat(r,c) / n;
        }
    }
    return retMat;
}

matrix_data_t gaps::algo::loglikelihood(const TwoWayMatrix &D,
const TwoWayMatrix &S, const TwoWayMatrix &AP)
{
    double chi2 = 0.0;
    for (unsigned r = 0; r < D.nRow(); ++r)
    {
        for (unsigned c = 0; c < D.nCol(); ++c)
        {
            chi2 += GAPS_SQ(D(r,c) - AP(r,c)) / GAPS_SQ(S(r,c));
        }
    }
    return chi2;
}

bool gaps::algo::isRowZero(const RowMatrix &mat, unsigned row)
{
    return gaps::algo::sum(mat.getRow(row)) == 0;
}
    
bool gaps::algo::isColZero(const ColMatrix &mat, unsigned col)
{
    return gaps::algo::sum(mat.getCol(col)) == 0;
}

static double deltaLL_A(const TwoWayMatrix &D, const TwoWayMatrix &S,
const RowMatrix &P, const TwoWayMatrix &AP, unsigned row,
unsigned col1, double delta1, unsigned col2=0, double delta2=0.0,
bool twoChanges=false)
{
    double numer = 0.0, denom = 0.0, delLL = 0.0, d1 = 0.0, d2 = 0.0;
    for (unsigned j = 0; j < D.nCol(); ++j)
    {
        d1 = delta1 * P(col1,j);
        d2 = twoChanges ? delta2 * P(col2,j) : 0.0;
        numer = 2.0 * (D(row,j) - AP(row,j)) * (d1 + d2) - GAPS_SQ(d1 + d2);
        denom = 2.0 * GAPS_SQ(S(row,j));
        delLL += numer / denom;
    }
    return delLL;
}

static double deltaLL_P(const TwoWayMatrix &D, const TwoWayMatrix &S,
const ColMatrix &A, const TwoWayMatrix &AP, unsigned col,
unsigned row1, double delta1, unsigned row2=0, double delta2=0.0,
bool twoChanges=false)
{
    double numer = 0.0, denom = 0.0, delLL = 0.0, d1 = 0.0, d2 = 0.0;
    for (unsigned i = 0; i < D.nRow(); ++i)
    {
        d1 = delta1 * A(i,row1);
        d2 = twoChanges ? delta2 * A(i,row2) : 0.0;
        numer = 2.0 * (D(i,col) - AP(i,col)) * (d1 + d2) - GAPS_SQ(d1 + d2);
        denom = 2.0 * GAPS_SQ(S(i,col));
        delLL += numer / denom;
    }
    return delLL;
}

double gaps::algo::deltaLL(const MatrixChange &ch, const TwoWayMatrix &D,
const TwoWayMatrix &S, const ColMatrix &A, const RowMatrix &P,
const TwoWayMatrix &AP)
{
    // change in A matrix
    if (ch.label == 'A' && ch.nChanges == 2 && ch.row1 != ch.row2)
    {
        return deltaLL_A(D, S, P, AP, ch.row1, ch.col1, ch.delta1)
            + deltaLL_A(D, S, P, AP, ch.row2, ch.col2, ch.delta2);
    }
    else if (ch.label == 'A')
    {
        return deltaLL_A(D, S, P, AP, ch.row1, ch.col1, ch.delta1, ch.col2,
            ch.delta2, ch.nChanges == 2);
    }

    // change in P matrix
    if (ch.label == 'P' && ch.nChanges == 2 && ch.col1 != ch.col2)
    {
        return deltaLL_P(D, S, A, AP, ch.col1, ch.row1, ch.delta1)
            + deltaLL_P(D, S, A, AP, ch.col2, ch.row2, ch.delta2);
    }
    else if (ch.label == 'P')
    {
        return deltaLL_P(D, S, A, AP, ch.col1, ch.row1, ch.delta1, ch.row2,
            ch.delta2, ch.nChanges == 2);
    }
}

static AlphaParameters alphaParameters_A(const TwoWayMatrix &D,
const TwoWayMatrix &S, const RowMatrix &P, const TwoWayMatrix &AP,
unsigned row, unsigned col1, unsigned col2=0, bool twoChanges=false)
{
    double s = 0.0, su = 0.0, ratio = 0.0, p2 = 0.0;
    for (unsigned j = 0; j < D.nCol(); ++j)
    {
        p2 = twoChanges ? P(col2,j) : 0.0;
        ratio = (P(col1,j) - p2) / S(row,j);
        s += GAPS_SQ(ratio);
        su += ratio * (D(row,j) - AP(row,j)) / S(row,j);
    }
    return AlphaParameters(s, su);
}

static AlphaParameters alphaParameters_P(const TwoWayMatrix &D,
const TwoWayMatrix &S, const ColMatrix &A, const TwoWayMatrix &AP,
unsigned col, unsigned row1, unsigned row2=0, bool twoChanges=false)
{
    double s = 0.0, su = 0.0, ratio = 0.0, a2 = 0.0;
    for (unsigned i = 0; i < D.nRow(); ++i)
    {
        a2 = twoChanges ? A(i,row2) : 0.0;
        ratio = (A(i,row1) - a2) / S(i,col);
        s += GAPS_SQ(ratio);
        su += ratio * (D(i,col) - AP(i,col)) / S(i,col);
    }
    return AlphaParameters(s, su);
}

AlphaParameters gaps::algo::alphaParameters(const MatrixChange &ch,
const TwoWayMatrix &D, const TwoWayMatrix &S, const ColMatrix &A,
const RowMatrix &P, const TwoWayMatrix &AP)
{
    // change in A matrix
    AlphaParameters p(0,0);
    if (ch.label == 'A' && ch.nChanges == 2 && ch.row1 != ch.row2)
    {
        p = alphaParameters_A(D, S, P, AP, ch.row1, ch.col1)
            + alphaParameters_A(D, S, P, AP, ch.row2, ch.col2);
    }
    else if (ch.label == 'A')
    {
        p = alphaParameters_A(D, S, P, AP, ch.row1, ch.col1, ch.col2,
            ch.nChanges == 2);
    }

    // change in P matrix
    if (ch.label == 'P' && ch.nChanges == 2 && ch.col1 != ch.col2)
    {
        p = alphaParameters_P(D, S, A, AP, ch.col1, ch.row1)
            + alphaParameters_P(D, S, A, AP, ch.col2, ch.row2);
    }
    else if (ch.label == 'P')
    {
        p = alphaParameters_P(D, S, A, AP, ch.col1, ch.row1, ch.row2,
            ch.nChanges == 2);
    }
    return p;
}
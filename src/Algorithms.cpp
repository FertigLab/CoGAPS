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
            meanTerm = meanMat(r,c) * meanMat(r,c) / nUpdates;
            retMat(r,c) = std::sqrt((stdMat(r,c) - meanTerm) / (nUpdates - 1));
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
            meanTerm = meanMat(r,c) * meanMat(r,c) / nUpdates;
            retMat(r,c) = std::sqrt((stdMat(r,c) - meanTerm) / (nUpdates - 1));
        }
    }
    return retMat;
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

bool gaps::algo::isRowZero(const RowMatrix &mat, unsigned row)
{
    return gaps::algo::sum(mat.getRow(row)) == 0;
}
    
bool gaps::algo::isColZero(const ColMatrix &mat, unsigned col)
{
    return gaps::algo::sum(mat.getCol(col)) == 0;
}

static double deltaLL_A(unsigned row, unsigned col, double delta,
const TwoWayMatrix &D, const TwoWayMatrix &S, const ColMatrix &A,
const RowMatrix &P, const TwoWayMatrix &AP)
{
    double numer = 0.0, denom = 0.0, delLL = 0.0;
    for (unsigned j = 0; j < D.nCol(); ++j)
    {
        numer = 2 * delta * (D(row,j) - AP(row,j)) * P(col,j) - GAPS_SQ(delta * P(col,j));
        denom = 2 * GAPS_SQ(S(row,j));
        delLL += numer / denom;
    }
    return delLL;
}

static double deltaLL_P(unsigned row, unsigned col, double delta,
const TwoWayMatrix &D, const TwoWayMatrix &S, const ColMatrix &A,
const RowMatrix &P, const TwoWayMatrix &AP)
{
    double numer = 0.0, denom = 0.0, delLL = 0.0;
    for (unsigned i = 0; i < D.nRow(); ++i)
    {
        numer = 2 * delta * (D(i,col) - AP(i,col)) * A(i,row) - GAPS_SQ(delta * A(i,row));
        denom = 2 * GAPS_SQ(S(i,col));
        delLL += numer / denom;
    }
    return delLL;
}

double gaps::algo::deltaLL(const MatrixChange &ch, const TwoWayMatrix &D,
const TwoWayMatrix &S, const ColMatrix &A, const RowMatrix &P,
const TwoWayMatrix &AP)
{
    if (ch.nChanges == 1 && ch.label == 'A')
    {
        return deltaLL_A(ch.row1, ch.col1, ch.delta1, D, S, A, P, AP);
    }
    else if (ch.nChanges == 1 && ch.label == 'P')
    {
        return deltaLL_P(ch.row1, ch.col1, ch.delta1, D, S, A, P, AP);
    }
    else if (ch.nChanges == 2 && ch.label == 'A')
    {
        return deltaLL_A(ch.row1, ch.col1, ch.delta1, D, S, A, P, AP)
            + deltaLL_A(ch.row2, ch.col2, ch.delta2, D, S, A, P, AP);
    }
    else
    {
        return deltaLL_P(ch.row1, ch.col1, ch.delta1, D, S, A, P, AP)
            + deltaLL_P(ch.row2, ch.col2, ch.delta2, D, S, A, P, AP);
    }
}

static AlphaParameters alphaParameters_A(unsigned row, unsigned col,
const TwoWayMatrix &D, const TwoWayMatrix &S, const RowMatrix &P,
const TwoWayMatrix &AP)
{
    double s = 0.0, su = 0.0, ratio = 0.0;
    for (unsigned j = 0; j < D.nCol(); ++j)
    {
        ratio = P(col,j) / S(row,j);
        s += GAPS_SQ(ratio);
        su += P(col,j) * (D(row,j) - AP(row,j)) / GAPS_SQ(S(row,j));
    }
    return AlphaParameters(s, su);
}

static AlphaParameters alphaParameters_P(unsigned row, unsigned col,
const TwoWayMatrix &D, const TwoWayMatrix &S, const ColMatrix &A,
const TwoWayMatrix &AP)
{
    double s = 0.0, su = 0.0, ratio = 0.0;
    for (unsigned i = 0; i < D.nRow(); ++i)
    {
        ratio = A(i,row) / S(i,col);
        s += GAPS_SQ(ratio);
        su += A(i,row) * (D(i,col) - AP(i,col)) / GAPS_SQ(S(i,col));
    }
    return AlphaParameters(s, su);
}

AlphaParameters gaps::algo::alphaParameters(const MatrixChange &ch,
const TwoWayMatrix &D, const TwoWayMatrix &S, const ColMatrix &A,
const RowMatrix &P, const TwoWayMatrix &AP)
{
    if (ch.nChanges == 1 && ch.label == 'A')
    {
        return alphaParameters_A(ch.row1, ch.col1, D, S, P, AP);
    }
    else if (ch.nChanges == 1 && ch.label == 'P')
    {
        return alphaParameters_P(ch.row1, ch.col1, D, S, A, AP);
    }
    else if (ch.nChanges == 2 && ch.label == 'A')
    {
        return alphaParameters_A(ch.row1, ch.col1, D, S, P, AP)
            + alphaParameters_A(ch.row2, ch.col2, D, S, P, AP);
    }
    else
    {
        return alphaParameters_P(ch.row1, ch.col1, D, S, A, AP)
            + alphaParameters_P(ch.row2, ch.col2, D, S, A, AP);
    }
}

bool gaps::algo::checkAPMatrix(const TwoWayMatrix &AP, const ColMatrix &A,
const RowMatrix &P)
{
    double sum = 0.0, diff = 0.0;
    for (unsigned r = 0; r < AP.nRow(); ++r)
    {
        for (unsigned c = 0; c < AP.nCol(); ++c)
        {
            sum = 0.0;
            for (unsigned k = 0; k < A.nCol(); ++k)
            {
                sum += A(r,k) * P(k,c);
            }

            diff = std::abs(AP(r,c) - sum);
            if (diff > 1E-10)
            {
                return false;
            }
        }
    }
    return true;
}
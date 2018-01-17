#include "Algorithms.h"
#include "Matrix.h"

#include <algorithm>

#define GAPS_SQ(x) ((x) * (x))

float gaps::algo::sum(const Vector &vec)
{
    float sum = 0.f;
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        sum += vec[i];
    }
    return sum;
}

Vector gaps::algo::scalarMultiple(const Vector &vec, float n)
{
    Vector retVec(vec.size());
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        retVec[i] = vec[i] * n;
    }
    return retVec;
}

Vector gaps::algo::squaredScalarMultiple(const Vector &vec, float n)
{
    Vector retVec(vec.size());
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        retVec[i] = GAPS_SQ(vec[i] * n);
    }
    return retVec;
}

Vector gaps::algo::scalarDivision(const Vector &vec, float n)
{
    Vector temp(vec.size());
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        temp[i] = vec[i] / n;
    }
    return temp;
}

Vector gaps::algo::squaredScalarDivision(const Vector &vec, float n)
{
    Vector retVec(vec.size());
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        retVec[i] = GAPS_SQ(vec[i] / n);
    }
    return retVec;
}

bool gaps::algo::isRowZero(const RowMatrix &mat, unsigned row)
{
    return gaps::algo::sum(mat.getRow(row)) == 0;
}
    
bool gaps::algo::isColZero(const ColMatrix &mat, unsigned col)
{
    return gaps::algo::sum(mat.getCol(col)) == 0;
}

// horribly slow, don't call often
void gaps::algo::matrixMultiplication(TwoWayMatrix &C, const ColMatrix &A,
const RowMatrix &B)
{
    for (unsigned i = 0; i < C.nRow(); ++i)
    {
        for (unsigned j = 0; j < C.nCol(); ++j)
        {
            float sum = 0.0;
            for (unsigned k = 0; k < A.nCol(); ++k)
            {
                sum += A(i,k) * B(k,j);
            }
            C.set(i, j, sum);
        }
    }
}

float gaps::algo::loglikelihood(const TwoWayMatrix &D,
const TwoWayMatrix &S, const TwoWayMatrix &AP)
{
    float chi2 = 0.f;
    for (unsigned i = 0; i < D.nRow(); ++i)
    {
        for (unsigned j = 0; j < D.nCol(); ++j)
        {
            chi2 += GAPS_SQ(D(i,j) - AP(i,j)) / GAPS_SQ(S(i,j));
        }
    }
    return chi2;
}

static float deltaLL_comp(const Vector &D, const Vector &S,
const Vector &AP, const Vector &other, float delta)
{
    float d = 0.f, delLL = 0.f;
    for (unsigned i = 0; i < D.size(); ++i)
    {
        d = delta * other[i];
        delLL += (2.f * d * (D[i] - AP[i]) - d * d) / (2.f * S[i] * S[i]);
    }
    return delLL;
}

static float deltaLL_comp(const Vector &D, const Vector &S,
const Vector &AP, const Vector &other1, float delta1,
const Vector &other2, float delta2)
{
    float d = 0.f, delLL = 0.f;
    for (unsigned i = 0; i < D.size(); ++i)
    {
        d = delta1 * other1[i] + delta2 * other2[i];
        delLL += (2.f * d * (D[i] - AP[i]) - d * d) / (2.f * S[i] * S[i]);
    }
    return delLL;
}

float gaps::algo::deltaLL(const MatrixChange &ch, const TwoWayMatrix &D,
const TwoWayMatrix &S, const ColMatrix &A, const RowMatrix &P,
const TwoWayMatrix &AP)
{
    if (ch.label == 'A' && ch.nChanges == 1)
    {
        return deltaLL_comp(D.getRow(ch.row1), S.getRow(ch.row1),
            AP.getRow(ch.row1), P.getRow(ch.col1), ch.delta1);
    }
    else if (ch.label == 'A' && ch.row1 != ch.row2)
    {
        return deltaLL_comp(D.getRow(ch.row1), S.getRow(ch.row1),
            AP.getRow(ch.row1), P.getRow(ch.col1), ch.delta1)
            + deltaLL_comp(D.getRow(ch.row2), S.getRow(ch.row2),
            AP.getRow(ch.row2), P.getRow(ch.col2), ch.delta2);
    }
    else if (ch.label == 'A')
    {
        return deltaLL_comp(D.getRow(ch.row1), S.getRow(ch.row1),
            AP.getRow(ch.row1), P.getRow(ch.col1), ch.delta1,
            P.getRow(ch.col2), ch.delta2);
    }
    else if (ch.label == 'P' && ch.nChanges == 1)
    {
        return deltaLL_comp(D.getCol(ch.col1), S.getCol(ch.col1),
            AP.getCol(ch.col1), A.getCol(ch.row1), ch.delta1);
    }
    else if (ch.label == 'P' && ch.col1 != ch.col2)
    {
        return deltaLL_comp(D.getCol(ch.col1), S.getCol(ch.col1),
            AP.getCol(ch.col1), A.getCol(ch.row1), ch.delta1)
            + deltaLL_comp(D.getCol(ch.col2), S.getCol(ch.col2),
            AP.getCol(ch.col2), A.getCol(ch.row2), ch.delta2);
    }
    else
    {
        return deltaLL_comp(D.getCol(ch.col1), S.getCol(ch.col1),
            AP.getCol(ch.col1), A.getCol(ch.row1), ch.delta1,
            A.getCol(ch.row2), ch.delta2);
    }
}

static AlphaParameters alphaParameters_comp(const Vector &D,
const Vector &S, const Vector &AP, const Vector &other)
{
    float s = 0.f, su = 0.f, ratio = 0.f;
    for (unsigned i = 0; i < D.size(); ++i)
    {
        ratio = other[i] / S[i];
        s += ratio * ratio;
        su += ratio * (D[i] - AP[i]) / S[i];
    }
    return AlphaParameters(s, su);
}

static AlphaParameters alphaParameters_comp(const Vector &D,
const Vector &S, const Vector &AP, const Vector &other1,
const Vector &other2)
{
    float s = 0.f, su = 0.f, ratio = 0.f;
    for (unsigned i = 0; i < D.size(); ++i)
    {
        ratio = (other1[i] - other2[i]) / S[i];
        s += ratio * ratio;
        su += ratio * (D[i] - AP[i]) / S[i];
    }
    return AlphaParameters(s, su);
}

AlphaParameters gaps::algo::alphaParameters(const MatrixChange &ch,
const TwoWayMatrix &D, const TwoWayMatrix &S, const ColMatrix &A,
const RowMatrix &P, const TwoWayMatrix &AP)
{
    if (ch.label == 'A' && ch.nChanges == 1)
    {
        return alphaParameters_comp(D.getRow(ch.row1), S.getRow(ch.row1),
            AP.getRow(ch.row1), P.getRow(ch.col1));
    }
    else if (ch.label == 'A' && ch.row1 != ch.row2)
    {
        return alphaParameters_comp(D.getRow(ch.row1), S.getRow(ch.row1),
            AP.getRow(ch.row1), P.getRow(ch.col1))
            + alphaParameters_comp(D.getRow(ch.row2), S.getRow(ch.row2),
            AP.getRow(ch.row2), P.getRow(ch.col2));
    }
    else if (ch.label == 'A')
    {
        return alphaParameters_comp(D.getRow(ch.row1), S.getRow(ch.row1),
            AP.getRow(ch.row1), P.getRow(ch.col1), P.getRow(ch.col2));
    }
    else if (ch.label == 'P' && ch.nChanges == 1)
    {
        return alphaParameters_comp(D.getCol(ch.col1), S.getCol(ch.col1),
            AP.getCol(ch.col1), A.getCol(ch.row1));
    }
    else if (ch.label == 'P' && ch.col1 != ch.col2)
    {
        return alphaParameters_comp(D.getCol(ch.col1), S.getCol(ch.col1),
            AP.getCol(ch.col1), A.getCol(ch.row1))
            + alphaParameters_comp(D.getCol(ch.col2), S.getCol(ch.col2),
            AP.getCol(ch.col2), A.getCol(ch.row2));
    }
    else
    {
        return alphaParameters_comp(D.getCol(ch.col1), S.getCol(ch.col1),
            AP.getCol(ch.col1), A.getCol(ch.row1), A.getCol(ch.row2));
    }
}
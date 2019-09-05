#include "HybridMatrix.h"
#include "HybridVector.h"
#include "Matrix.h"
#include "../utils/Archive.h"
#include "Vector.h"

HybridMatrix::HybridMatrix(unsigned nrow, unsigned ncol)
    :
mRows(nrow, Vector(ncol)),
mCols(ncol, HybridVector(nrow)),
mNumRows(nrow),
mNumCols(ncol)
{}
    
unsigned HybridMatrix::nRow() const
{
    return mNumRows;
}

unsigned HybridMatrix::nCol() const
{
    return mNumCols;
}

void HybridMatrix::add(unsigned i, unsigned j, float v)
{
    GAPS_ASSERT(i < mNumRows);
    GAPS_ASSERT(j < mNumCols);
    mRows[i][j] += v;
    mCols[j].add(i, v);
}

void HybridMatrix::set(unsigned i, unsigned j, float v)
{
    GAPS_ASSERT(i < mNumRows);
    GAPS_ASSERT(j < mNumCols);
    mRows[i][j] = v;
    mCols[j].set(i, v);
}

float HybridMatrix::operator()(unsigned i, unsigned j) const
{
    return mRows[i][j];
}

const Vector& HybridMatrix::getRow(unsigned n) const
{
    return mRows[n];
}

const HybridVector& HybridMatrix::getCol(unsigned n) const
{
    return mCols[n];
}

Matrix HybridMatrix::getMatrix() const
{
    Matrix mat(mNumRows, mNumCols);
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            mat(i,j) = this->operator()(i,j);
        }
    }
    return mat;
}

void HybridMatrix::operator=(const Matrix &mat)
{
    GAPS_ASSERT(mNumRows == mat.nRow());
    GAPS_ASSERT(mNumCols == mat.nCol());

    for (unsigned i = 0; i < mNumRows; ++i)
    {
        for (unsigned j = 0; j < mNumCols; ++j)
        {
            mRows[i][j] = mat(i,j);
            mCols[j].add(i, -1.f * mCols[j][i]);
            mCols[j].add(i, mat(i,j));
        }
    }    
}

Archive& operator<<(Archive &ar, const HybridMatrix &vec)
{
    ar << vec.mNumRows << vec.mNumCols;
    for (unsigned i = 0; i < vec.mRows.size(); ++i)
    {
        ar << vec.mRows[i];
    }
    for (unsigned i = 0; i < vec.mCols.size(); ++i)
    {
        ar << vec.mCols[i];
    }
    return ar;
}

Archive& operator>>(Archive &ar, HybridMatrix &vec)
{
    unsigned nr = 0, nc = 0;
    ar >> nr >> nc;

    GAPS_ASSERT(vec.mNumRows == nr);
    GAPS_ASSERT(vec.mNumCols == nc);

    for (unsigned i = 0; i < vec.mRows.size(); ++i)
    {
        ar >> vec.mRows[i];
    }
    for (unsigned i = 0; i < vec.mCols.size(); ++i)
    {
        ar >> vec.mCols[i];
    }
    return ar;
}


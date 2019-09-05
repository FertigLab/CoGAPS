#include "Matrix.h"
#include "SparseVector.h"
#include "../file_parser/FileParser.h"
#include "../file_parser/MatrixElement.h"
#include "../utils/Archive.h"
#include "../utils/GapsAssert.h"
#include "Vector.h"

#include <algorithm>
#include <iterator>

Matrix::Matrix() : mNumRows(0), mNumCols(0) {}

Matrix::Matrix(unsigned nrow, unsigned ncol)
    :
mCols(ncol, Vector(nrow)),
mNumRows(nrow),
mNumCols(ncol)
{}

void Matrix::pad(float val)
{
    for (unsigned j = 0; j < mNumCols; ++j)
    {
        mCols[j].pad(val);
    }
}

// constructor from data set read in as a matrix
Matrix::Matrix(const Matrix &mat, bool genesInCols, bool subsetGenes,
std::vector<unsigned> indices)
{
#ifdef GAPS_DEBUG
    for (unsigned i = 0; i < indices.size(); ++i)
    {
        GAPS_ASSERT_MSG(indices[i] > 0,
            "index 0 detected in subset: R indices should start at 1\n");
    }
#endif

    bool subsetData = !indices.empty();

    unsigned nGenes = (subsetData && subsetGenes)
        ? indices.size()
        : genesInCols ? mat.nCol() : mat.nRow();
    unsigned nSamples = (subsetData && !subsetGenes)
        ? indices.size()
        : genesInCols ? mat.nRow() : mat.nCol();
    
    for (unsigned j = 0; j < nSamples; ++j)
    {
        mCols.push_back(Vector(nGenes));
        for (unsigned i = 0; i < nGenes; ++i)
        {
            unsigned dataRow = (subsetData && (subsetGenes != genesInCols))
                ? indices[genesInCols ? j : i] - 1
                : genesInCols ? j : i;

            unsigned dataCol = (subsetData && (subsetGenes == genesInCols))
                ? indices[genesInCols ? i : j] - 1
                : genesInCols ? i : j;

            GAPS_ASSERT(i < nGenes);
            mCols[j][i] = mat(dataRow, dataCol);
        }
    }
    mNumRows = nGenes;
    mNumCols = nSamples;
}

// constructor from data set given as a file path
Matrix::Matrix(const std::string &path, bool genesInCols, bool subsetGenes,
std::vector<unsigned> indices)
{
#ifdef GAPS_DEBUG
    for (unsigned i = 0; i < indices.size(); ++i)
    {
        GAPS_ASSERT_MSG(indices[i] > 0,
            "index 0 detected in subset: R indices should start at 1\n");
    }
#endif

    FileParser fp(path);

    // calculate the number of rows and columns
    bool subsetData = !indices.empty();
    mNumRows = (subsetData && subsetGenes) // nGenes
        ? indices.size()
        : genesInCols ? fp.nCol() : fp.nRow();
    mNumCols = (subsetData && !subsetGenes) // nSamples
        ? indices.size()
        : genesInCols ? fp.nRow() : fp.nCol();

    // allocate space for the data
    for (unsigned j = 0; j < mNumCols; ++j)
    {
        mCols.push_back(Vector(mNumRows));
    }

    // read from file
    if (!subsetData)
    {
        while (fp.hasNext())
        {
            MatrixElement e(fp.getNext());
            unsigned row = genesInCols ? e.col : e.row;
            unsigned col = genesInCols ? e.row : e.col;
            this->operator()(row, col) = e.value;
        }
    }
    else
    {
        std::sort(indices.begin(), indices.end());
        while (fp.hasNext())
        {
            MatrixElement e(fp.getNext());
            unsigned searchIndex = 1 + ((subsetGenes != genesInCols) ? e.row : e.col);
            std::vector<unsigned>::iterator pos = 
                std::lower_bound(indices.begin(), indices.end(), searchIndex);
        
            // this index is included in the subset
            if (pos != indices.end() && *pos == searchIndex)
            {
                unsigned row = subsetGenes
                    ? std::distance(indices.begin(), pos)
                    : genesInCols ? e.col : e.row;
                unsigned col = !subsetGenes
                    ? std::distance(indices.begin(), pos)
                    : genesInCols ? e.row : e.col;
                this->operator()(row, col) = e.value;
            }
        }
    }
}

unsigned Matrix::nRow() const
{
    return mNumRows;
}

unsigned Matrix::nCol() const
{
    return mNumCols;
}

float Matrix::operator()(unsigned i, unsigned j) const
{
    GAPS_ASSERT_MSG(i < mNumRows, i << " : " << mNumRows);
    GAPS_ASSERT_MSG(j < mNumCols, j << " : " << mNumCols);
    return mCols[j][i];
}

float& Matrix::operator()(unsigned i, unsigned j)
{
    GAPS_ASSERT_MSG(i < mNumRows, i << " : " << mNumRows);
    GAPS_ASSERT_MSG(j < mNumCols, j << " : " << mNumCols);
    return mCols[j][i];
}

Vector& Matrix::getCol(unsigned col)
{
    GAPS_ASSERT_MSG(col < mNumCols, col << " : " << mNumCols);
    return mCols[col];
}

const Vector& Matrix::getCol(unsigned col) const
{
    GAPS_ASSERT_MSG(col < mNumCols, col << " : " << mNumCols);
    return mCols[col];
}

bool Matrix::empty() const
{
    return mNumRows == 0;
}

Matrix Matrix::getMatrix() const
{
    return *this;
}

Archive& operator<<(Archive &ar, const Matrix &mat)
{
    ar << mat.mNumRows << mat.mNumCols;
    for (unsigned j = 0; j < mat.mNumCols; ++j)
    {
        ar << mat.mCols[j];
    }
    return ar;
}

Archive& operator>>(Archive &ar, Matrix &mat)
{
    unsigned nr = 0, nc = 0;
    ar >> nr >> nc;
    GAPS_ASSERT(nr == mat.mNumRows);
    GAPS_ASSERT(nc == mat.mNumCols);

    for (unsigned j = 0; j < mat.mNumCols; ++j)
    {
        ar >> mat.mCols[j];
    }
    return ar;
}
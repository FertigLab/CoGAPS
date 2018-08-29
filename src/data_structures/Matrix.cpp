#include "Matrix.h"
#include "../file_parser/CsvParser.h"
#include "../file_parser/MatrixElement.h"
#include "../utils/GapsAssert.h"

ColMatrix::ColMatrix(unsigned nrow, unsigned ncol)
    :
mNumRows(nrow), mNumCols(ncol)
{
    allocate();
}

ColMatrix::ColMatrix(const ColMatrix &mat, bool transpose, bool partitionRows,
const std::vector<unsigned> &indices)
{
    if (indices.size() <= 1)
    {
        mNumRows = transpose ? mat.nCol() : mat.nRow();
        mNumCols = transpose ? mat.nRow() : mat.nCol();
        allocate();

        for (unsigned i = 0; i < mNumRows; ++i)
        {
            for (unsigned j = 0; j < mNumCols; ++j)
            {
                this->operator()(i,j) = transpose ? mat(j,i) : mat(i,j);
            }
        }
    }
    else
    {
        // create matrix
        mNumRows = partitionRows ? indices.size() : (transpose ? mat.nCol() : mat.nRow());
        mNumCols = partitionRows ? (transpose ? mat.nRow() : mat.nCol()) : indices.size();
        allocate();

        // fill matrix, TODO use binary search on indices
        for (unsigned i = 0; i < mat.nRow(); ++i)
        {
            for (unsigned j = 0; j < mat.nCol(); ++j)
            {
                unsigned dataIndex = transpose ? (partitionRows ? j : i) : (partitionRows ? i : j);
                std::vector<unsigned>::const_iterator pos = std::find(indices.begin(), indices.end(), dataIndex);
                if (pos != indices.end())
                {
                    unsigned index = std::distance(indices.begin(), pos);
                    unsigned row = partitionRows ? index : (transpose ? j : i);
                    unsigned col = partitionRows ? (transpose ? i : j) : index;
                    this->operator()(row,col) = mat(i,j);
                }
            }
        }
    }
}

// apply transpose first, then partition either rows or columns
ColMatrix::ColMatrix(const std::string &path, bool transpose, bool partitionRows,
const std::vector<unsigned> &indices)
{
    FileParser parser(path);
    if (indices.size() <= 1)
    {
        mNumRows = transpose ? parser.nCol() : parser.nRow();
        mNumCols = transpose ? parser.nRow() : parser.nCol();
        allocate();

        while (parser.hasNext())
        {
            MatrixElement e(parser.getNext());
            unsigned row = transpose ? e.col : e.row;
            unsigned col = transpose ? e.row : e.col;
            this->operator()(row,col) = e.value;
        }
    }
    else
    {
        // create matrix
        mNumRows = partitionRows ? indices.size() : (transpose ? parser.nCol() : parser.nRow());
        mNumCols = partitionRows ? (transpose ? parser.nRow() : parser.nCol()) : indices.size();
        allocate();

        // fill matrix
        while (parser.hasNext())
        {
            MatrixElement e(parser.getNext());
            unsigned dataIndex = transpose ? (partitionRows ? e.col : e.row) : (partitionRows ? e.row : e.col);
            std::vector<unsigned>::const_iterator pos = std::find(indices.begin(), indices.end(), dataIndex);
            if (pos != indices.end())
            {
                unsigned index = std::distance(indices.begin(), pos);
                unsigned row = partitionRows ? index : (transpose ? e.col : e.row);
                unsigned col = partitionRows ? (transpose ? e.row : e.col) : index;
                this->operator()(row, col) = e.value;
            }
        }
    }
}

void ColMatrix::allocate()
{
    for (unsigned i = 0; i < mNumCols; ++i)
    {
        mData.push_back(Vector(mNumRows));
    }
}

unsigned ColMatrix::nRow() const
{
    return mNumRows;
}

unsigned ColMatrix::nCol() const
{
    return mNumCols;
}

ColMatrix ColMatrix::operator*(float val) const
{
    ColMatrix mat(*this);
    for (unsigned i = 0; i < mNumCols; ++i)
    {
        mat.mData[i] *= val;
    }
    return mat;
}

ColMatrix ColMatrix::operator/(float val) const
{
    ColMatrix mat(*this);
    for (unsigned i = 0; i < mNumCols; ++i)
    {
        mat.mData[i] /= val;
    }
    return mat;
}

float& ColMatrix::operator()(unsigned r, unsigned c)
{
    return mData[c][r];
}

float ColMatrix::operator()(unsigned r, unsigned c) const
{
    return mData[c][r];
}

Vector& ColMatrix::getCol(unsigned col)
{
    return mData[col];
}

const Vector& ColMatrix::getCol(unsigned col) const
{
    return mData[col];
}

float* ColMatrix::colPtr(unsigned col)
{
    GAPS_ASSERT(col < mNumCols);
    return mData[col].ptr();
}

const float* ColMatrix::colPtr(unsigned col) const
{
    GAPS_ASSERT(col < mNumCols);
    return mData[col].ptr();
}

Archive& operator<<(Archive &ar, ColMatrix &mat)
{
    ar << mat.mNumRows << mat.mNumCols;
    for (unsigned i = 0; i < mat.mNumCols; ++i)
    {
        ar << mat.mData[i];
    }
    return ar;
}

Archive& operator>>(Archive &ar, ColMatrix &mat)
{
    // should already by allocated
    unsigned nr = 0, nc = 0;
    ar >> nr >> nc;
    GAPS_ASSERT(nr == mat.mNumRows);
    GAPS_ASSERT(nc == mat.mNumCols);

    // read in data
    for (unsigned i = 0; i < mat.mNumCols; ++i)
    {
        ar >> mat.mData[i];
    }
    return ar;   
}
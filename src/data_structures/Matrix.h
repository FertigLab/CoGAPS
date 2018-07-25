#ifndef __COGAPS_MATRIX_H__
#define __COGAPS_MATRIX_H__

#include "../Archive.h"
#include "../file_parser/FileParser.h"
#include "Vector.h"

#include <algorithm>
#include <vector>

/****************************** MATRIX INTERFACE ******************************/

class RowMatrix;
class ColMatrix;
typedef RowMatrix Matrix; // when we don't care about row/col major order

template <class T>
class GenericMatrix;

template <class T>
Archive& operator<<(Archive &ar, GenericMatrix<T> &samp);

template <class T>
Archive& operator>>(Archive &ar, GenericMatrix<T> &samp);


template <class T>
class GenericMatrix
{
protected:

    std::vector<Vector> mData;
    unsigned mNumRows, mNumCols;

    T* impl();
    const T* impl() const;

public:

    GenericMatrix(unsigned nrow, unsigned ncol);

    GenericMatrix(const Matrix &mat, bool transpose=false);
    GenericMatrix(const std::string &path, bool transpose=false);

    GenericMatrix(const Matrix &mat, bool transpose, bool partitionRows,
        const std::vector<unsigned> &indices);
    GenericMatrix(const std::string &path, bool transpose, bool partitionRows,
        const std::vector<unsigned> &indices);

    unsigned nRow() const;
    unsigned nCol() const;

    T operator*(float val) const;
    T operator/(float val) const;
    T pmax(float scale, float max) const;

    friend Archive& operator<< <T> (Archive &ar, GenericMatrix &mat);
    friend Archive& operator>> <T> (Archive &ar, GenericMatrix &mat);
};

class RowMatrix : public GenericMatrix<RowMatrix>
{
private:

    friend class GenericMatrix;
    void allocate();

public:

    RowMatrix(unsigned nrow, unsigned ncol);

    RowMatrix(const Matrix &mat, bool transpose=false);
    RowMatrix(const std::string &path, bool transpose=false);

    RowMatrix(const Matrix &mat, bool transpose, bool partitionRows,
        const std::vector<unsigned> &indices);
    RowMatrix(const std::string &path, bool transpose, bool partitionRows,
        const std::vector<unsigned> &indices);

    // simplifies assigning to/from standard Matrix (RowMatrix)
    RowMatrix& operator=(const ColMatrix &mat);

    // for single element access - do not loop over elements with this
    float& operator()(unsigned r, unsigned c);
    float operator()(unsigned r, unsigned c) const;

    // for convenience when doing non-performance critical math
    Vector& getRow(unsigned row);
    const Vector& getRow(unsigned row) const;

    // raw pointer allows for high performance looping over rows with SIMD
    float* rowPtr(unsigned row);
    const float* rowPtr(unsigned row) const;
};

class ColMatrix : public GenericMatrix<ColMatrix>
{
private:

    friend class GenericMatrix;
    void allocate();

public:

    ColMatrix(unsigned nrow, unsigned ncol);

    ColMatrix(const Matrix &mat, bool transpose=false);
    ColMatrix(const std::string &path, bool transpose=false);

    ColMatrix(const Matrix &mat, bool transpose, bool partitionRows,
        const std::vector<unsigned> &indices);
    ColMatrix(const std::string &path, bool transpose, bool partitionRows,
        const std::vector<unsigned> &indices);

    // simplifies assigning to/from standard Matrix (RowMatrix)
    ColMatrix& operator=(const RowMatrix &mat);

    // for single element access - do not loop over elements with this
    float& operator()(unsigned r, unsigned c);
    float operator()(unsigned r, unsigned c) const;

    // for convenience when doing non-performance critical math
    Vector& getCol(unsigned col);
    const Vector& getCol(unsigned col) const;

    // raw pointer allows for high performance looping over rows with SIMD
    float* colPtr(unsigned col);
    const float* colPtr(unsigned col) const;
};

/******************** IMPLEMENTATION OF TEMPLATED FUNCTIONS *******************/

template <class T>
T* GenericMatrix<T>::impl()
{
    return static_cast<T*>(this);
}

template <class T>
const T* GenericMatrix<T>::impl() const
{
    return static_cast<const T*>(this);
}

template <class T>
GenericMatrix<T>::GenericMatrix(unsigned nrow, unsigned ncol)
{
    mNumRows = nrow;
    mNumCols = ncol;
    impl()->allocate();
}

template <class T>
GenericMatrix<T>::GenericMatrix(const Matrix &mat, bool transpose)
{
    mNumRows = transpose ? mat.nCol() : mat.nRow();
    mNumCols = transpose ? mat.nRow() : mat.nCol();
    impl()->allocate();

    for (unsigned i = 0; i < mNumRows; ++i)
    {
        for (unsigned j = 0; j < mNumCols; ++j)
        {
            impl()->operator()(i,j) = transpose ? mat(j,i) : mat(i,j);
        }
    }
}

template <class T>
GenericMatrix<T>::GenericMatrix(const std::string &path, bool transpose)
{
    FileParser p(path);
    mNumRows = transpose ? p.nCol() : p.nRow();
    mNumCols = transpose ? p.nRow() : p.nCol();
    impl()->allocate();

    while (p.hasNext())
    {
        MatrixElement e(p.getNext());
        unsigned row = transpose ? e.col() : e.row();
        unsigned col = transpose ? e.row() : e.col();
        impl()->operator()(row,col) = e.value();
    }
}

// apply transpose first, then partition either rows or columns
template <class T>
GenericMatrix<T>::GenericMatrix(const Matrix &mat, bool transpose,
bool partitionRows, const std::vector<unsigned> &indices)
{
    // create matrix
    mNumRows = partitionRows ? indices.size() : (transpose ? mat.nCol() : mat.nRow());
    mNumCols = partitionRows ? (transpose ? mat.nRow() : mat.nCol()) : indices.size();
    impl()->allocate();

    // fill matrix
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
                impl()->operator()(row,col) = mat(i,j);
            }
        }
    }
}
    
// apply transpose first, then partition either rows or columns
template <class T>
GenericMatrix<T>::GenericMatrix(const std::string &path, bool transpose,
bool partitionRows, const std::vector<unsigned> &indices)
{
    // create matrix
    FileParser p(path);
    mNumRows = partitionRows ? indices.size() : (transpose ? p.nCol() : p.nRow());
    mNumCols = partitionRows ? (transpose ? p.nRow() : p.nCol()) : indices.size();
    impl()->allocate();

    // fill matrix
    while (p.hasNext())
    {
        MatrixElement e(p.getNext());
        unsigned dataIndex = transpose ? (partitionRows ? e.col : e.row) : (partitionRows ? e.row : e.col);
        std::vector<unsigned>::const_iterator pos = std::find(indices.begin(), indices.end(), dataIndex);
        if (pos != indices.end())
        {
            unsigned index = std::distance(indices.begin(), pos);
            unsigned row = partitionRows ? index : (transpose ? e.col : e.row);
            unsigned col = partitionRows ? (transpose ? e.row : e.col) : index;
            impl()->operator()(row, col) = e.value;
        }
    }
}

template <class T>
unsigned GenericMatrix<T>::nRow() const
{
    return mNumRows;
}

template <class T>
unsigned GenericMatrix<T>::nCol() const
{
    return mNumCols;
}

template <class T>
T GenericMatrix<T>::operator*(float val) const
{
    T mat(*this);
    for (unsigned i = 0; i < mData.size(); ++i)
    {
        mat.mData *= val;
    }
    return mat;
}

template <class T>
T GenericMatrix<T>::operator/(float val) const
{
    T mat(*this);
    for (unsigned i = 0; i < mData.size(); ++i)
    {
        mat.mData /= val;
    }
    return mat;
}

template <class T>
T GenericMatrix<T>::pmax(float scale, float max) const
{
    T mat(mNumRows, mNumCols);
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        for (unsigned j = 0; j < mNumCols; ++j)
        {
            mat(i,j) = std::max(impl()->operator()(i,j) * scale, max);
        }
    }
    return mat;
}

template <class T>
Archive& operator<<(Archive &ar, GenericMatrix<T> &mat)
{
    for (unsigned i = 0; i < mat.mData.size(); ++i)
    {
        ar << mat.mData[i];
    }
    return ar;
}

template <class T>
Archive& operator>>(Archive &ar, GenericMatrix<T> &mat)
{
    for (unsigned i = 0; i < mat.mData.size(); ++i)
    {
        ar >> mat.mData[i];
    }
    return ar;
}

#endif // __COGAPS_MATRIX_H__

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
typedef Matrix RowMatrix; // when we don't care about row/col major order

template <class T>
class GenericMatrix
{
protected:

    std::vector<Vector> mData;
    unsigned mNumRows, mNumCols;

    T* impl();

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

    friend Archive& operator<< <T>(Archive &ar, GenericMatrix &mat);
    friend Archive& operator>> <T>(Archive &ar, GenericMatrix &mat);
};

class RowMatrix : public Matrix<RowMatrix>
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
    float& operator()(unsigned r, unsigned c) {return mData[r][c];}
    float operator()(unsigned r, unsigned c) const {return mData[r][c];}

    // for convenience when doing non-performance critical math
    Vector& getRow(unsigned row) {return mData[row];}
    const Vector& getRow(unsigned row) const {return mData[row];}

    // raw pointer allows for high performance looping over rows with SIMD
    float* rowPtr(unsigned row) {return mData[row].ptr();}
    const float* rowPtr(unsigned row) const {return mData[row].ptr();}
};

class ColMatrix : public Matrix<ColMatrix>
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
    float& operator()(unsigned r, unsigned c) {return mData[c][r];}
    float operator()(unsigned r, unsigned c) const {return mData[c][r];}

    // for convenience when doing non-performance critical math
    Vector& getCol(unsigned col) {return mData[col];}
    const Vector& getCol(unsigned col) const {return mData[col];}

    // raw pointer allows for high performance looping over rows with SIMD
    float* colPtr(unsigned col) {return mData[col].ptr();}
    const float* colPtr(unsigned col) const {return mData[col].ptr();}
};

/******************** IMPLEMENTATION OF TEMPLATED FUNCTIONS *******************/

template <class T>
GenericMatrix::GenericMatrix(unsigned nrow, unsigned ncol)
    :
mNumRows(nrow), mNumCols(ncol)
{
    impl()->allocate();
}

GenericMatrix::GenericMatrix(const Matrix &mat, bool transpose)
    :
mNumRows(transpose ? mat.nCol() : mat.nRow()),
mNumCols(transpose ? mat.nRow() : mat.nCol())
{
    impl()->allocate();
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        for (unsigned j = 0; j < mNumCols; ++j)
        {
            impl()->operator()(i,j) = transpose ? mat(j,i) : mat(i,j);
        }
    }
}

GenericMatrix::GenericMatrix(const std::string &path, bool transpose)
{
    FileParser p(path);
    mNumRows = transpose ? p.nCol() : p.nRow();
    mNumCols = transpose ? p.nRow() : p.nCol();
    allocate();

    while (p.hasNext())
    {
        MatrixElement e(p.getNext());
        unsigned row = transpose ? e.col() : e.row();
        unsigned col = transpose ? e.row() : e.col();
        impl()->operator()(row,col) = e.value();
    }
}

    GenericMatrix(const Matrix &mat, bool transpose, bool partitionRows,
        const std::vector<unsigned> &indices);
    GenericMatrix(const std::string &path, bool transpose, bool partitionRows,
        const std::vector<unsigned> &indices);

    unsigned nRow() const;
    unsigned nCol() const;

    T operator*(float val) const;
    T operator/(float val) const;
    T& operator=(const ColMatrix &mat);
    T pmax(float scale, float max) const;

    friend Archive& operator<<(Archive &ar, T &mat);
    friend Archive& operator>>(Archive &ar, T &mat);


#if 0
// forward declarations
class RowMatrix;
class ColMatrix;

class RowMatrix
{
private:

    std::vector<Vector> mRows;


    void allocate();

public:

    RowMatrix(unsigned nrow, unsigned ncol);
    explicit RowMatrix(const ColMatrix &mat);

    RowMatrix(const RowMatrix &mat, bool transpose=false);
    RowMatrix(const std::string &path, bool transpose=false);

    RowMatrix(const RowMatrix &mat, bool transpose, bool partitionRows,
        const std::vector<unsigned> &indices);
    RowMatrix(const std::string &p, bool transpose, bool partitionRows,
        const std::vector<unsigned> &indices);

    RowMatrix operator*(float val) const;
    RowMatrix operator/(float val) const;
    RowMatrix& operator=(const ColMatrix &mat);
    RowMatrix pmax(float scale, float max) const;

    friend Archive& operator<<(Archive &ar, RowMatrix &mat);
    friend Archive& operator>>(Archive &ar, RowMatrix &mat);
};

class ColMatrix
{
private:

    std::vector<Vector> mCols;
    unsigned mNumRows, mNumCols;

public:

    ColMatrix(unsigned nrow, unsigned ncol);
    explicit ColMatrix(const RowMatrix &mat);

    ColMatrix(ColMatrix mat, bool transpose=false);
    ColMatrix(const std::string &path, bool transpose=false);

    ColMatrix(ColMatrix mat, bool transpose, bool partitionRows,
        const std::vector<unsigned> &indices);
    ColMatrix(const std::string &p, bool transpose, bool partitionRows,
        const std::vector<unsigned> &indices);

    ColMatrix operator*(float val) const;
    ColMatrix operator/(float val) const;
    ColMatrix& operator=(const RowMatrix &mat);
    ColMatrix pmax(float scale, float max) const;

    friend Archive& operator<<(Archive &ar, ColMatrix &mat);
    friend Archive& operator>>(Archive &ar, ColMatrix &mat);
};

#endif

#endif // __COGAPS_MATRIX_H__

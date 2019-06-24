#ifndef __COGAPS_SPARSE_MATRIX_H__
#define __COGAPS_SPARSE_MATRIX_H__

#include "SparseVector.h"

#include <string>
#include <vector>

class Archive;
class Matrix;

// no random access, all data is const, can only access with iterator
// over a given column
class SparseMatrix
{
public:
    SparseMatrix(const Matrix &mat, bool genesInCols, bool subsetGenes,
        std::vector<unsigned> indices);
    SparseMatrix(const std::string &path, bool genesInCols, bool subsetGenes,
        std::vector<unsigned> indices);
    unsigned nRow() const;
    unsigned nCol() const;
    const SparseVector& getCol(unsigned n) const;
    void operator=(const Matrix &mat);
    friend Archive& operator<<(Archive &ar, const SparseMatrix &vec);
    friend Archive& operator>>(Archive &ar, SparseMatrix &vec);
private:
    std::vector<SparseVector> mCols;
    unsigned mNumRows;
    unsigned mNumCols;
};

#endif
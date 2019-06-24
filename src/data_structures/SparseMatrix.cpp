#include "SparseMatrix.h"
#include "Matrix.h"
#include "../file_parser/FileParser.h"
#include "../utils/Archive.h"
#include "../utils/GapsAssert.h"

#include <algorithm>
#include <iterator>

// constructor from data set read in as a matrix
SparseMatrix::SparseMatrix(const Matrix &mat, bool genesInCols,
bool subsetGenes, std::vector<unsigned> indices)
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
        std::vector<float> values;
        for (unsigned i = 0; i < nGenes; ++i)
        {
            unsigned dataRow = (subsetData && (subsetGenes != genesInCols))
                ? indices[genesInCols ? j : i] - 1
                : genesInCols ? j : i;

            unsigned dataCol = (subsetData && (subsetGenes == genesInCols))
                ? indices[genesInCols ? i : j] - 1
                : genesInCols ? i : j;

            values.push_back(mat(dataRow, dataCol));
        }
        mCols.push_back(SparseVector(values));
    }
    mNumRows = nGenes;
    mNumCols = nSamples;
}

// constructor from data set given as a file path
SparseMatrix::SparseMatrix(const std::string &path, bool genesInCols,
bool subsetGenes, std::vector<unsigned> indices)
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
        mCols.push_back(SparseVector(mNumRows));
    }

    // read from file
    if (!subsetData)
    {
        while (fp.hasNext())
        {
            MatrixElement e(fp.getNext());
            unsigned row = genesInCols ? e.col : e.row;
            unsigned col = genesInCols ? e.row : e.col;
            if (e.value > 0.f)
            {
                mCols[col].insert(row, e.value);
            }
        }
    }
    else
    {
        std::sort(indices.begin(), indices.end());
        while (fp.hasNext())
        {
            MatrixElement e(fp.getNext());
            if (e.value > 0.f)
            {
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
                    mCols[col].insert(row, e.value);
                }
            }
        }
    }
}

unsigned SparseMatrix::nRow() const
{
    return mNumRows;
}

unsigned SparseMatrix::nCol() const
{
    return mNumCols;
}

const SparseVector& SparseMatrix::getCol(unsigned n) const
{
    return mCols[n];   
}

void SparseMatrix::operator=(const Matrix &mat)
{
    GAPS_ASSERT(mNumRows == mat.nRow());
    GAPS_ASSERT(mNumCols == mat.nCol());

    // no nice way to overwrite this
    mCols.clear();
    for (unsigned j = 0; j < mNumCols; ++j)
    {
        mCols.push_back(SparseVector(mat.getCol(j)));
    }
}

Archive& operator<<(Archive &ar, const SparseMatrix &vec)
{
    ar << vec.mNumRows << vec.mNumCols;
    for (unsigned j = 0; j < vec.mNumCols; ++j)
    {
        ar << vec.mCols[j];
    }
    return ar;
}

Archive& operator>>(Archive &ar, SparseMatrix &vec)
{
    unsigned nr = 0, nc = 0;
    ar >> nr >> nc;
    GAPS_ASSERT(nr == vec.mNumRows);
    GAPS_ASSERT(nc == vec.mNumCols);

    for (unsigned j = 0; j < vec.mNumCols; ++j)
    {
        ar >> vec.mCols[j];
    }
    return ar;
}

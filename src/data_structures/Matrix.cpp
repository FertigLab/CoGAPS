#include "Matrix.h"
#include "../file_parser/CsvParser.h"
#include "../file_parser/MatrixElement.h"

template<class Matrix>
static Rcpp::NumericMatrix convertToRMatrix(const Matrix &mat)
{
    Rcpp::NumericMatrix rmat(mat.nRow(), mat.nCol());
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            rmat(i,j) = mat(i,j);
        }
    }
    return rmat;
}

template<class MatA, class MatB>
static void copyMatrix(MatA &dest, const MatB &source)
{
    for (unsigned i = 0; i < source.nRow(); ++i)
    {
        for (unsigned j = 0; j < source.nCol(); ++j)
        {
            dest(i,j) = source(i,j);
        }
    }
}

/****************************** ROW MATRIX *****************************/

RowMatrix::RowMatrix(unsigned nrow, unsigned ncol)
: mNumRows(nrow), mNumCols(ncol)
{
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mRows.push_back(Vector(mNumCols));
    }
}

RowMatrix::RowMatrix(const Rcpp::NumericMatrix &rmat)
: mNumRows(rmat.nrow()), mNumCols(rmat.ncol())
{
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mRows.push_back(Vector(mNumCols));
        for (unsigned j = 0; j < mNumCols; ++j)
        {
            this->operator()(i,j) = rmat(i,j);
        }
    }
}

RowMatrix& RowMatrix::operator=(const ColMatrix &mat)
{
    copyMatrix(*this, mat);
    return *this;
}

Rcpp::NumericMatrix RowMatrix::rMatrix() const
{
    return convertToRMatrix(*this);
}

RowMatrix RowMatrix::operator/(float val) const
{
    RowMatrix mat(mNumRows, mNumCols);
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        mat.getRow(i) = mRows[i] / val;
    }
    return mat;
}

Archive& operator<<(Archive &ar, RowMatrix &mat)
{
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        ar << mat.mRows[i];
    }
    return ar;
}

Archive& operator>>(Archive &ar, RowMatrix &mat)
{
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        ar >> mat.mRows[i];
    }
    return ar;
}

void RowMatrix::writeToCsv(const std::string &path)
{
    std::ofstream outputFile;
    outputFile.open(path.c_str());
    outputFile << "\"\"";
    for (unsigned i = 0; i < mNumCols; ++i)
    {
        outputFile << ",\"Col" << i << "\"";
    }
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        outputFile << "\"Row" << i << "\"";
        for (unsigned j = 0; j < mNumCols; ++j)
        {
            outputFile << "," << mRows[i][j];
        }
        outputFile << "\n";
    }
    outputFile.close();
}

void RowMatrix::writeToTsv(const std::string &path)
{
    std::ofstream outputFile;
    outputFile.open(path.c_str());
    outputFile << "\"\"";
    for (unsigned i = 0; i < mNumCols; ++i)
    {
        outputFile << "\t\"Col" << i << "\"";
    }
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        outputFile << "\"Row" << i << "\"";
        for (unsigned j = 0; j < mNumCols; ++j)
        {
            outputFile << "\t" << mRows[i][j];
        }
        outputFile << "\n";
    }
    outputFile.close();
}

void RowMatrix::writeToMtx(const std::string &path)
{
    std::ofstream outputFile;
    outputFile.open(path.c_str());
    outputFile << "%%\n";
    outputFile << mNumRows << " " << mNumCols << " " << (mNumRows * mNumCols);
    outputFile << "\n";
    for (unsigned j = 1; j <= mNumCols; ++j)
    {
        for (unsigned i = 1; i <= mNumRows; ++i)
        {
            outputFile << i << " " << j << " " << mRows[i][j] << "\n";
        }
    }
    outputFile.close();
}

/**************************** COLUMN MATRIX ****************************/

ColMatrix::ColMatrix(unsigned nrow, unsigned ncol)
: mNumRows(nrow), mNumCols(ncol)
{
    for (unsigned i = 0; i < mNumCols; ++i)
    {
        mCols.push_back(Vector(mNumRows));
    }
}

ColMatrix::ColMatrix(const Rcpp::NumericMatrix &rmat)
: mNumRows(rmat.nrow()), mNumCols(rmat.ncol())
{
    for (unsigned j = 0; j < mNumCols; ++j)
    {
        mCols.push_back(Vector(mNumRows));
        for (unsigned i = 0; i < mNumRows; ++i)
        {
            this->operator()(i,j) = rmat(i,j);
        }
    }
}

ColMatrix ColMatrix::operator/(float val) const
{
    ColMatrix mat(mNumRows, mNumCols);
    for (unsigned j = 0; j < mNumCols; ++j)
    {
        mat.getCol(j) = mCols[j] / val;
    }
    return mat;
}

ColMatrix& ColMatrix::operator=(const RowMatrix &mat)
{
    copyMatrix(*this, mat);
    return *this;
}

Rcpp::NumericMatrix ColMatrix::rMatrix() const
{
    return convertToRMatrix(*this);
}

Archive& operator<<(Archive &ar, ColMatrix &mat)
{
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        ar << mat.mCols[j];
    }
    return ar;
}

Archive& operator>>(Archive &ar, ColMatrix &mat)
{
    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        ar >> mat.mCols[j];
    }
    return ar;
}

void ColMatrix::writeToCsv(const std::string &path)
{
    std::ofstream outputFile;
    outputFile.open(path.c_str());
    outputFile << "\"\"";
    for (unsigned i = 0; i < mNumCols; ++i)
    {
        outputFile << ",\"Col" << i << "\"";
    }
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        outputFile << "\"Row" << i << "\"";
        for (unsigned j = 0; j < mNumCols; ++j)
        {
            outputFile << "," << mCols[j][i];
        }
        outputFile << "\n";
    }
    outputFile.close();
}

void ColMatrix::writeToTsv(const std::string &path)
{
    std::ofstream outputFile;
    outputFile.open(path.c_str());
    outputFile << "\"\"";
    for (unsigned i = 0; i < mNumCols; ++i)
    {
        outputFile << "\t\"Col" << i << "\"";
    }
    for (unsigned i = 0; i < mNumRows; ++i)
    {
        outputFile << "\"Row" << i << "\"";
        for (unsigned j = 0; j < mNumCols; ++j)
        {
            outputFile << "\t" << mCols[j][i];
        }
        outputFile << "\n";
    }
    outputFile.close();
}

void ColMatrix::writeToMtx(const std::string &path)
{
    std::ofstream outputFile;
    outputFile.open(path.c_str());
    outputFile << "%%\n";
    outputFile << mNumRows << " " << mNumCols << " " << (mNumRows * mNumCols);
    outputFile << "\n";
    for (unsigned j = 1; j <= mNumCols; ++j)
    {
        for (unsigned i = 1; i <= mNumRows; ++i)
        {
            outputFile << i << " " << j << " " << mCols[j][i] << "\n";
        }
    }
    outputFile.close();
}

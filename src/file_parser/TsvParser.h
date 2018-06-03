#ifndef __COGAPS_TSV_PARSER_H__
#define __COGAPS_TSV_PARSER_H__

#include "MatrixElement.h"

#include <fstream>
#include <string>

class TsvParser
{
private:

    std::ifstream mFile;

    unsigned mNumRows;
    unsigned mNumCols;

    unsigned mCurrentRow;
    unsigned mCurrentCol;

public:

    explicit TsvParser(const std::string &path);

    unsigned nRow() const { return mNumRows; }
    unsigned nCol() const { return mNumCols; }

    bool hasNext();
    MatrixElement getNext();
};

#endif
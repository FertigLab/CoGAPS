#ifndef __COGAPS_MTX_PARSER_H__
#define __COGAPS_MTX_PARSER_H__

#include "MatrixElement.h"

#include <fstream>
#include <string>

class MtxParser
{
private:

    std::ifstream mFile;

    unsigned mNumRows;
    unsigned mNumCols;

public:

    explicit MtxParser(const std::string &path);

    unsigned nRow() const { return mNumRows; }
    unsigned nCol() const { return mNumCols; }

    bool hasNext();
    MatrixElement getNext();
};

#endif
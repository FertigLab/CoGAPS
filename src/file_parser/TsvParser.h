#ifndef __COGAPS_TSV_PARSER_H__
#define __COGAPS_TSV_PARSER_H__

#include "FileParser.h"
#include "MatrixElement.h"

#include <fstream>
#include <string>

class TsvParser : public AbstractFileParser
{
private:

    std::ifstream mFile;

    unsigned mNumRows;
    unsigned mNumCols;

    unsigned mCurrentRow;
    unsigned mCurrentCol;

    TsvParser(const TsvParser &p); // don't allow copies
    TsvParser& operator=(const TsvParser &p); // don't allow copies

public:

    explicit TsvParser(const std::string &path);
    ~TsvParser() {}

    unsigned nRow() const { return mNumRows; }
    unsigned nCol() const { return mNumCols; }

    bool hasNext();
    MatrixElement getNext();
};

#endif
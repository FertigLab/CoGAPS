#ifndef __COGAPS_GCT_PARSER_H__
#define __COGAPS_GCT_PARSER_H__

#include "FileParser.h"
#include "MatrixElement.h"

#include <fstream>
#include <string>

class GctParser : public AbstractFileParser
{
private:

    std::ifstream mFile;

    unsigned mNumRows;
    unsigned mNumCols;

    unsigned mCurrentRow;
    unsigned mCurrentCol;

    GctParser(const GctParser &p); // don't allow copies
    GctParser& operator=(const GctParser &p); // don't allow copies

public:

    explicit GctParser(const std::string &path);
    ~GctParser() {}

    unsigned nRow() const { return mNumRows; }
    unsigned nCol() const { return mNumCols; }

    bool hasNext();
    MatrixElement getNext();
};

#endif
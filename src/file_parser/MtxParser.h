#ifndef __COGAPS_MTX_PARSER_H__
#define __COGAPS_MTX_PARSER_H__

#include "FileParser.h"
#include "MatrixElement.h"

#include <cstdio>
#include <string>

class MtxParser : public AbstractFileParser
{
public:

    explicit MtxParser(const std::string &path);
    ~MtxParser() { mFile.close(); }

    unsigned nRow() const { return mNumRows; }
    unsigned nCol() const { return mNumCols; }

    bool hasNext();
    MatrixElement getNext();

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    FILE *mFile;

    unsigned mNumRows;
    unsigned mNumCols;

    MtxParser(const MtxParser &p); // don't allow copies
    MtxParser& operator=(const MtxParser &p); // don't allow copies
};

#endif
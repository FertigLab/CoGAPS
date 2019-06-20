#ifndef __COGAPS_MTX_PARSER_H__
#define __COGAPS_MTX_PARSER_H__

#include "FileParser.h"
#include "MatrixElement.h"

#include <fstream>
#include <string>

class MtxParser : public AbstractFileParser
{
public:

    explicit MtxParser(const std::string &path);
    ~MtxParser();

    unsigned nRow() const;
    unsigned nCol() const;

    bool hasNext();
    MatrixElement getNext();

private:

    std::ifstream mFile;

    unsigned mNumRows;
    unsigned mNumCols;

    MtxParser(const MtxParser &p); // don't allow copies
    MtxParser& operator=(const MtxParser &p); // don't allow copies

    void checkFileState() const;
};

#endif
#ifndef __COGAPS_CSV_PARSER_H__
#define __COGAPS_CSV_PARSER_H__

#include "FileParser.h"
#include "MatrixElement.h"

#include <fstream>
#include <string>

class CsvParser : public AbstractFileParser
{
private:

    std::ifstream mFile;

    unsigned mNumRows;
    unsigned mNumCols;

    unsigned mCurrentRow;
    unsigned mCurrentCol;

public:

    explicit CsvParser(const std::string &path);
    ~CsvParser() {}

    unsigned nRow() const { return mNumRows; }
    unsigned nCol() const { return mNumCols; }

    bool hasNext();
    MatrixElement getNext();
};

#endif
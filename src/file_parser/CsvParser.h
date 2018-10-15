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

    CsvParser(const CsvParser &p); // don't allow copies
    CsvParser& operator=(const CsvParser &p); // don't allow copies

public:

    explicit CsvParser(const std::string &path);
    ~CsvParser() { mFile.close(); }

    unsigned nRow() const { return mNumRows; }
    unsigned nCol() const { return mNumCols; }

    bool hasNext();
    MatrixElement getNext();
};

#endif
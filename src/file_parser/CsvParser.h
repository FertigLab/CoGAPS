#ifndef __COGAPS_CSV_PARSER_H__
#define __COGAPS_CSV_PARSER_H__

#include "MatrixElement.h"

#include <fstream>
#include <vector>
#include <string>

class CsvParser
{
private:

    std::ifstream mFile;

    unsigned mCurrentRow;
    unsigned mCurrentCol;

public:

    CsvParser(const std::string &path);

    bool hasNext();
    MatrixElement getNext();

    static MatrixDimensions getDimensions(const std::string &path);
};

#endif
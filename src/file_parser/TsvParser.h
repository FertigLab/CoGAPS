#ifndef __COGAPS_TSV_PARSER_H__
#define __COGAPS_TSV_PARSER_H__

#include "MatrixElement.h"

#include <fstream>
#include <vector>
#include <string>

class TsvParser
{
private:

    std::ifstream mFile;

    unsigned mCurrentRow;
    unsigned mCurrentCol;

public:

    TsvParser(const std::string &path);

    bool hasNext();
    MatrixElement getNext();

    static MatrixDimensions getDimensions(const std::string &path);
};

#endif
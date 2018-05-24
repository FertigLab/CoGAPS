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

    std::vector<std::string> mRowNames;
    std::vector<std::string> mColNames;

    unsigned mCurrentRow;
    unsigned mCurrentCol;

public:

    CsvParser(const std::string &path);

    bool hasNext() const;
    MatrixElement getNext() const; 
};

#endif
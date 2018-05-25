#ifndef __COGAPS_CSV_PARSER_H__
#define __COGAPS_CSV_PARSER_H__

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

    // read through whole file once, store row/col names - gives dimensions
    CsvParser(const std::string &path);

    unsigned nRow() const { return mRowNames.size(); }
    unsigned nCol() const { return mColNames.size(); }

    bool hasNextRow();
    std::vector<float> getNextRow(); 
    void skipNextRow();

    bool hasNextCol();
    std::vector<float> getNextCol();
    void skipNextCol();
};

#endif
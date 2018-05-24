#include "CsvParser.h"

#include <iostream>

// TODO need to parse by rows - otherwise it would be neccesary to
// know dimensions beforehand

// open file, read column names
CsvParser::CsvParser(const std::string &path) : mCurrentRow(0), mCurrentCol(0)
{
    mFile.open(path.c_str());

    std::string line;
    std::getline(mFile, line, ','); // read first entry (blank)

    std::size_t pos;
    std::getline(mFile, line, ',');
    while ((pos = line.find('\n')) == std::string::npos)
    {
        mColNames.push_back(line);
        std::getline(mFile, line, ',');
    }
    mColNames.push_back(line.substr(0,pos));
    mRowNames.push_back(line.substr(pos+1));
}

bool CsvParser::hasNext()
{
    return mFile.peek() != EOF;
}

MatrixElement CsvParser::getNext()
{
    std::string line;
    std::getline(mFile, line, ',');

    std::size_t pos;
    if ((pos = line.find('\n')) != std::string::npos)
    {
        if (pos + 1 < line.size())
        {
            mRowNames.push_back(line.substr(pos + 1));
        }
        unsigned col = mCurrentCol;
        mCurrentCol = 0;
        return MatrixElement(mCurrentRow++, col, line.substr(0, pos));
    }
    else
    {
        return MatrixElement(mCurrentRow, mCurrentCol++, line);
    }
}
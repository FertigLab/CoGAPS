#include "CsvParser.h"

// open file, read column names
CsvParser::CsvParser(const std::string &path)
{
    mFile.open(path.c_str());

    std::string line;
    std::getline(mFile, line); // read first entry (blank)
    
    while (!std::isdigit(mFile.peek()))
    {
        std::getline(mFile, line)
        mColNames.push_back(line);
    }
    mRowNames.push_back(mColNames.back());
    mColNames.pop_back(); // read one too far
}

bool CsvParser::hasNext() const
{
    return mFile.peek() != EOF;
}

MatrixElement CsvParser::getNext() const
{
    int c = mFile.peek();
    
    std::string line;
    std::getline(mFile, line);
    if (std::isdigit(c)) // matrix element
    {
        Rcout << line << '\n';
        return MatrixElement(0,0,0.f);
    }
    else // row name
    {
        mRowNames.push_back(line);
        return getNext();
    }
}
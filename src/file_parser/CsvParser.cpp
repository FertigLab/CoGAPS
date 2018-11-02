#include "CsvParser.h"
#include "../utils/GapsAssert.h"

#include <fstream>
#include <iostream>

// read through whole file once, store row/col names - gives dimensions
// open file, throw away column names
CsvParser::CsvParser(const std::string &path) :  mNumRows(0), mNumCols(0),
mCurrentRow(0), mCurrentCol(0)
{
    // open file stream
    std::ifstream file_str(path.c_str());

    // read first entry (blank)
    std::string line;
    std::getline(file_str, line, ',');
    if (file_str.eof() || file_str.fail())
    {
        GAPS_ERROR("Invalid CSV file");
    }

    // get col size
    std::size_t pos;
    do
    {
        std::getline(file_str, line, ',');
        if (file_str.eof() || file_str.fail())
        {
            GAPS_ERROR("Invalid CSV file");
        }
        ++mNumCols;
    }
    while ((pos = line.find('\n')) == std::string::npos);

    // get row size
    ++mNumRows; // acount for current row
    while (file_str.peek() != EOF)
    {
        // throw away data
        do
        {
            std::getline(file_str, line, ',');
        }
        while ((pos = line.find('\n')) == std::string::npos);

        // increment row number, ignore last newline in file
        if (pos + 1 < line.size())
        {
            ++mNumRows;
        }
    }

    // open member file stream and set at first element
    mFile.open(path.c_str());
    std::getline(mFile, line); // get rid of first line (column names)
}

bool CsvParser::hasNext()
{
    mFile >> std::ws;
    return mFile.peek() != EOF;
}

MatrixElement CsvParser::getNext()
{
    std::string line;
    std::size_t pos;
    std::getline(mFile, line, ',');
    if ((pos = line.find('\n')) != std::string::npos) // end of line
    {
        unsigned col = mCurrentCol;
        mCurrentCol = 0;
        return MatrixElement(mCurrentRow++, col, line.substr(0, pos));
    }
    if (std::isdigit(line[0]) != 0) // data
    {
        return MatrixElement(mCurrentRow, mCurrentCol++, line);
    }
    return getNext();
}

#include "CsvParser.h"
#include "../math/Algorithms.h"

#include <iostream>

// get the number of rows and cols in a csv file
MatrixDimension CsvParser::getDimensions(const std::string &path)
{
    // initialize struct that holds dimensions
    MatrixDimension dim(0,0);

    // open file stream
    std::ifstream str(path);

    // read first entry (blank)
    std::string line;
    std::getline(mFile, line, ',');

    // get col size
    std::size_t pos;
    do
    {
        std::getline(mFile, line, ',');
        dim.nCol++;
    }
    while ((pos = line.find('\n')) == std::string::npos);

    // get row size
    dim.nRow++; // acount for current row
    while (mFile.peek() != EOF)
    {
        // throw away data
        do
        {
            std::getline(mFile, line, ',');
        }
        while ((pos = line.find('\n')) == std::string::npos);

        // increment row number, ignore last newline in file
        if (pos + 1 < line.size())
        {
            dim.nRow++;
        }
    }
    return dim;
}

// read through whole file once, store row/col names - gives dimensions
// open file, read column names
CsvParser::CsvParser(const std::string &path) : mCurrentRow(0), mCurrentCol(0)
{
    mFile.open(path);
}

bool CsvParser::hasNext()
{
    return mFile.peek() != EOF;
}

MatrixElement CsvParser::getNext()
{
    std::string line;
    std::size_t pos;
    std::getline(mFile, line, ',');
    if ((pos = line.find('\n')) != std::string::npos) // end of line
    {
        return MatrixElement(mCurrentRow, mCurrentCol, line.substr(0, pos));
    }
    else if (std::isdigit(line[0])) // data
    {
        return MatrixElement(mCurrentRow, mCurrentCol, line);
    }
    else // row/col name
    {
        return getNext();
    }
}

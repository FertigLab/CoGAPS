#include "CsvParser.h"
#include "../math/Algorithms.h"

#include <fstream>
#include <iostream>

// get the number of rows and cols in a csv file
MatrixDimensions CsvParser::getDimensions(const std::string &path)
{
    // initialize struct that holds dimensions
    MatrixDimensions dim(0,0);

    // open file stream
    std::ifstream file_str(path.c_str());

    // read first entry (blank)
    std::string line;
    std::getline(file_str, line, ',');

    // get col size
    std::size_t pos;
    do
    {
        std::getline(file_str, line, ',');
        dim.nCol++;
    }
    while ((pos = line.find('\n')) == std::string::npos);

    // get row size
    dim.nRow++; // acount for current row
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
            dim.nRow++;
        }
    }
    return dim;
}

// read through whole file once, store row/col names - gives dimensions
// open file, read column names
CsvParser::CsvParser(const std::string &path) : mCurrentRow(0), mCurrentCol(0)
{
    mFile.open(path.c_str());
    std::string line;
    std::getline(mFile, line); // get rid of first line (column names)
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
        unsigned col = mCurrentCol;
        mCurrentCol = 0;
        return MatrixElement(mCurrentRow++, col, line.substr(0, pos));
    }
    else if (std::isdigit(line[0])) // data
    {
        return MatrixElement(mCurrentRow, mCurrentCol++, line);
    }
    else // row/col name
    {
        return getNext();
    }
}

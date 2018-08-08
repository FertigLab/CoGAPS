#include "GctParser.h"
#include "../GapsAssert.h"
#include "../math/Algorithms.h"

#include <fstream>
#include <iostream>

// read through whole file once, store row/col names - gives dimensions
// open file, throw away column names
GctParser::GctParser(const std::string &path) :  mNumRows(0), mNumCols(0),
mCurrentRow(0), mCurrentCol(0)
{
    mFile.open(path.c_str());

    // read first line
    std::string line;
    std::getline(mFile, line);
    if (mFile.eof() || mFile.fail())
    {
        GAPS_ERROR("Invalid GCT file");
    }
        
    // read the second line containing the dimensions
    std::getline(mFile, line);
    std::stringstream ss(line);

    // store dimensions
    ss >> mNumRows >> mNumCols;

    // throw away columns names
    std::getline(mFile, line);
}

bool GctParser::hasNext()
{
    mFile >> std::ws;
    return mFile.peek() != EOF;
}

MatrixElement GctParser::getNext()
{
    std::string line;
    std::size_t pos;
    std::getline(mFile, line, '\t');
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

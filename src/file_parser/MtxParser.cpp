#include "MatrixElement.h"
#include "MtxParser.h"

#include "../utils/GapsAssert.h"

#include <sstream>

MtxParser::MtxParser(const std::string &path) : mNumRows(0), mNumCols(0)
{
    mFile.open(path.c_str());

    // read first line
    std::string line;
    std::getline(mFile, line);
    if (mFile.eof() || mFile.fail())
    {
        GAPS_ERROR("Invalid MTX file");
    }

    // skip over comments
    while (line.find('%') != std::string::npos)
    {
        std::getline(mFile, line);
        if (mFile.eof() || mFile.fail())
        {
            GAPS_ERROR("Invalid MTX file");
        }
    }
    std::stringstream ss(line); // this line contains dimensions

    // store dimensions
    ss >> mNumRows >> mNumCols;
}

bool MtxParser::hasNext()
{
    mFile >> std::ws; // get rid of whitespace
    return mFile.peek() != EOF;
}

MatrixElement MtxParser::getNext()
{
    unsigned row = 0, col = 0;
    float val = 0.f;

    mFile >> row;
    mFile >> col;
    mFile >> val;

    return MatrixElement(row - 1, col - 1, val);
}

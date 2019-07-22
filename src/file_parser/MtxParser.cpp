#include "MatrixElement.h"
#include "MtxParser.h"

#include "../utils/GapsAssert.h"

#include <sstream>

MtxParser::MtxParser(const std::string &path) : mNumRows(0), mNumCols(0)
{
    mFile.open(path.c_str());

    // skip over comments
    std::string line = "%";
    while (line.find('%') != std::string::npos)
    {
        std::getline(mFile, line);
        checkFileState();
    }

    std::stringstream ss(line); // this line contains dimensions
    ss >> mNumRows >> mNumCols;
}

MtxParser::~MtxParser()
{
    mFile.close();
}

unsigned MtxParser::nRow() const
{
    return mNumRows;
}

unsigned MtxParser::nCol() const
{
    return mNumCols;
}

void MtxParser::checkFileState() const
{
    if (mFile.eof() || mFile.fail())
    {
        GAPS_ERROR("Invalid MTX file");
    }
}

bool MtxParser::hasNext()
{
    mFile >> std::ws; // get rid of whitespace
    return mFile.peek() != EOF;
}

MatrixElement MtxParser::getNext()
{
    unsigned row = 0;
    unsigned col = 0;
    std::string val;
    mFile >> row;
    mFile >> col;
    mFile >> val;
    return MatrixElement(row - 1, col - 1, val);
}

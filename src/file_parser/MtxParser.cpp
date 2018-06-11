#include "MtxParser.h"
#include "MatrixElement.h"

#include <sstream>

MtxParser::MtxParser(const std::string &path) : mNumRows(0), mNumCols(0)
{
    mFile.open(path.c_str());

    // skip over comments
    std::string line = "%";
    while (line.find('%') != std::string::npos)
    {
        std::getline(mFile, line);
    }
    std::stringstream ss(line); // this line contains dimensions

    // store dimensions
    ss >> mNumRows >> mNumCols;
}

bool MtxParser::hasNext()
{
    return mFile.peek() != EOF;
}

MatrixElement MtxParser::getNext()
{
    MatrixElement e(0, 0, 0.f);
    unsigned buff;
    mFile >> buff;
    e.dim[0] = buff - 1;
    mFile >> buff;
    e.dim[1] = buff - 1;
    mFile >> e.value;
    return e;
}

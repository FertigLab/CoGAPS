#include "CharacterDelimitedParser.h"
#include "../utils/GapsAssert.h"

#include <fstream>
#include <iostream>

static const std::string whitespace = " \r\n";

std::string ltrim(const std::string &s)
{
    std::size_t start = s.find_first_not_of(whitespace);
    return (start == std::string::npos) ? "" : s.substr(start);
}

std::string rtrim(const std::string &s)
{
    std::size_t end = s.find_last_not_of(whitespace);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

std::string trim(const std::string &s)
{
    return rtrim(ltrim(s));
}

std::string trimNewline(const std::string &s)
{
    std::size_t pos;
    if ((pos = s.find('\n')) == std::string::npos)
        return s;
    return s.substr(0, pos);
}

// read through whole file once, store row/col names - gives dimensions
// open file, throw away column names
CharacterDelimitedParser::CharacterDelimitedParser(const std::string &path, char delimiter)
    :
mNumRows(0), mNumCols(0), mCurrentRow(0), mCurrentCol(0), mRowNamesPresent(false),
mDelimiter(delimiter)
{
    // open file stream
    std::ifstream file_str(path.c_str());

    // read first entry
    std::string line;
    std::getline(file_str, line, mDelimiter);
    if (file_str.eof() || file_str.fail())
    {
        GAPS_ERROR("Invalid character delimited file");
    }

    // check if row names are given, if not count this first read as a column
    mRowNamesPresent = trim(line).empty();
    if (!mRowNamesPresent)
    {
        ++mNumCols;
        mColNames.push_back(trim(trimNewline(line)));
    }

    // get col size
    std::size_t pos;
    do
    {
        std::getline(file_str, line, mDelimiter);
        if (file_str.eof() || file_str.fail())
        {
            GAPS_ERROR("Invalid character delimited file");
        }
        mColNames.push_back(trim(trimNewline(line)));
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
            std::getline(file_str, line, mDelimiter);
        }
        while ((pos = line.find('\n')) == std::string::npos);

        if (pos + 1 < line.size())
        {
            ++mNumRows; // increment row number, ignore last newline in file
        }
    }

    // open member file stream and set at first element
    mFile.open(path.c_str());
    std::getline(mFile, line); // get rid of first line (column names)
}

CharacterDelimitedParser::~CharacterDelimitedParser()
{
    mFile.close();
}

bool CharacterDelimitedParser::hasNext()
{
    mFile >> std::ws;
    return mFile.peek() != EOF;
}

MatrixElement CharacterDelimitedParser::getNext()
{
    if (!mBuffer.empty())
    {
        std::string temp = mBuffer;
        mBuffer = "";
        return MatrixElement(mCurrentRow, mCurrentCol++, temp);
    }

    std::string line;
    std::size_t pos;
    std::getline(mFile, line, mDelimiter);
    if ((pos = line.find('\n')) != std::string::npos) // end of line
    {
        unsigned col = mCurrentCol;
        mCurrentCol = 0;
        if (mRowNamesPresent)
        {
            return MatrixElement(mCurrentRow++, col, line.substr(0, pos));
        }
        else
        {
            mBuffer = line.substr(pos + 1);
            return MatrixElement(mCurrentRow++, col, line.substr(0, pos));
        }
    }
    else
    {
        return MatrixElement(mCurrentRow, mCurrentCol++, line);
    }
    return getNext();
}

unsigned CharacterDelimitedParser::nRow() const
{
    return mNumRows;
}

unsigned CharacterDelimitedParser::nCol() const
{
    return mNumCols;
}

std::vector<std::string> CharacterDelimitedParser::rowNames() const
{
    return mRowNames;
}

std::vector<std::string> CharacterDelimitedParser::colNames() const
{
    return mColNames;
}

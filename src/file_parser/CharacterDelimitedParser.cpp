#include "CharacterDelimitedParser.h"
#include "../utils/GapsAssert.h"

#include <fstream>
#include <iostream>
#include <sstream>

static const std::string trimChars = " \r\n\"";

static std::string ltrim(const std::string &s)
{
    std::size_t start = s.find_first_not_of(trimChars);
    return (start == std::string::npos) ? "" : s.substr(start);
}

static std::string rtrim(const std::string &s)
{
    std::size_t end = s.find_last_not_of(trimChars);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

static std::string trim(const std::string &s)
{
    return rtrim(ltrim(s));
}

static std::string trimNewline(const std::string &s)
{
    std::size_t pos;
    if ((pos = s.find('\n')) == std::string::npos)
    {
        return s;
    }
    return s.substr(0, pos);
}

static std::vector<std::string> split(const std::string &s, char delimiter)
{
    std::vector<std::string> tokens;
    std::string temp;
    std::stringstream ss(s);
    while (std::getline(ss, temp, delimiter))
    {
        tokens.push_back(trim(temp));
    }
    return tokens;    
}

void CharacterDelimitedParser::checkFileState() const
{
    if (mFile.eof() || mFile.fail())
    {
        GAPS_ERROR("Invalid character delimited file");
    }
}

// read through whole file once, store row/col names and dimensions
CharacterDelimitedParser::CharacterDelimitedParser(const std::string &path, char delimiter, bool gctFormat)
    :
mNumRows(0), mNumCols(0), mCurrentRow(0), mCurrentCol(0), mRowNamesPresent(false),
mDelimiter(delimiter), mGctFormat(gctFormat)
{
    // read first entry
    mFile.open(path.c_str());
    std::string line;

    // if stored in gct format, read the second line with the dimensions
    if (mGctFormat)
    {
        std::getline(mFile, line);
        std::getline(mFile, line);
        checkFileState();
        std::stringstream ss(line);
        ss >> mNumRows >> mNumCols;
    }
    else // other formats need to walk the whole file to get the info
    {
        std::getline(mFile, line, mDelimiter);
        checkFileState();    

        // check if row names are given, if not count this first read as a column
        mRowNamesPresent = trim(line).empty();
        if (!mRowNamesPresent)
        {
            ++mNumCols;
            mColNames.push_back(trim(trimNewline(line)));
        }

        // find number of columns and record column names
        std::size_t pos;
        do
        {
            std::getline(mFile, line, mDelimiter);
            checkFileState();
            ++mNumCols;
            mColNames.push_back(trim(trimNewline(line)));
        }
        while ((pos = line.find('\n')) == std::string::npos);

        // find number of rows
        while (std::getline(mFile, line))
        {
            ++mNumRows;
        }
    }

    // reset file stream to beginning
    mFile.clear();
    mFile.seekg(0, std::ios::beg);
    std::getline(mFile, line); // get rid of first line (column names for csv/tsv)
    if (mGctFormat)
    {
        std::getline(mFile, line); // get rid of the file dimensions
        std::getline(mFile, line); // get rid of the column names
    }
    parseNextLine();
}

CharacterDelimitedParser::~CharacterDelimitedParser()
{
    mFile.close();
}

bool CharacterDelimitedParser::hasNext()
{
    if (mCurrentCol < mCurrentLine.size())
    {
        return true;
    }
    mFile >> std::ws;
    return mFile.peek() != EOF;
}

void CharacterDelimitedParser::parseNextLine()
{
    std::string fullLine;
    std::getline(mFile, fullLine);
    mCurrentLine = split(fullLine, mDelimiter);
    if (mRowNamesPresent)
    {
        mCurrentLine.erase(mCurrentLine.begin());
    }
    if (mGctFormat)
    {
        mCurrentLine.erase(mCurrentLine.begin(), mCurrentLine.begin() + 2);
    }
}

MatrixElement CharacterDelimitedParser::getNext()
{
    if (mCurrentCol < mCurrentLine.size())
    {
        unsigned colIndex = mCurrentCol;
        ++mCurrentCol;
        return MatrixElement(mCurrentRow, colIndex, mCurrentLine[colIndex]);
    }
    parseNextLine();
    ++mCurrentRow;
    mCurrentCol = 0;
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

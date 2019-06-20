#ifndef __COGAPS_CHARACTER_DELIMITED_PARSER_H__
#define __COGAPS_CHARACTER_DELIMITED_PARSER_H__

#include "FileParser.h"
#include "MatrixElement.h"

#include <fstream>
#include <vector>
#include <string>

class CharacterDelimitedParser : public AbstractFileParser
{
public:

    explicit CharacterDelimitedParser(const std::string &path, char delimiter, bool gctFormat=false);
    ~CharacterDelimitedParser();

    std::vector<std::string> rowNames() const;
    std::vector<std::string> colNames() const;
    unsigned nRow() const;
    unsigned nCol() const;
    bool hasNext();
    MatrixElement getNext();

private:

    std::ifstream mFile;

    std::vector<std::string> mRowNames;
    std::vector<std::string> mColNames;

    std::vector<std::string> mCurrentLine;

    unsigned mNumRows;
    unsigned mNumCols;

    unsigned mCurrentRow;
    unsigned mCurrentCol;

    bool mRowNamesPresent;
    char mDelimiter;
    bool mGctFormat;

    void parseNextLine();
    void checkFileState() const;

    CharacterDelimitedParser(const CharacterDelimitedParser &p); // don't allow copies
    CharacterDelimitedParser& operator=(const CharacterDelimitedParser &p); // don't allow copies
};

#endif // __COGAPS_CHARACTER_DELIMITED_PARSER_H__
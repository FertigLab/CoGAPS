#ifndef __COGAPS_CHARACTER_DELIMITED_PARSER_H__
#define __COGAPS_CHARACTER_DELIMITED_PARSER_H__

#include "FileParser.h"

#include <fstream>
#include <string>
#include <vector>

struct MatrixElement;

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
    CharacterDelimitedParser(const CharacterDelimitedParser &p); // don't allow copies
    CharacterDelimitedParser& operator=(const CharacterDelimitedParser &p); // don't allow copies
    void parseNextLine();
    void checkFileState() const;

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
};

#endif // __COGAPS_CHARACTER_DELIMITED_PARSER_H__
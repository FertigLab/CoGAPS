#ifndef __COGAPS_FILE_PARSER_H__
#define __COGAPS_FILE_PARSER_H__

#include "MatrixElement.h"

// file parser interface
class AbstractFileParser
{
public:

    static AbstractFileParser* create(const std::string &path);

    virtual ~AbstractFileParser() = 0;

    virtual unsigned nRow() const = 0;
    virtual unsigned nCol() const = 0;

    virtual bool hasNext() = 0;
    virtual MatrixElement getNext() = 0;
};

// wrap the pointer to the parser implementation
class FileParser
{
private:

    AbstractFileParser *mParser;

    FileParser(const FileParser &p); // don't allow copies
    FileParser& operator=(const FileParser &p); // don't allow copies

public:

    explicit FileParser(const std::string &path)
    {
        mParser = AbstractFileParser::create(path);
    }

    ~FileParser() { delete mParser; }

    unsigned nRow() const { return mParser->nRow(); }
    unsigned nCol() const { return mParser->nCol(); }

    bool hasNext() { return mParser->hasNext(); }
    MatrixElement getNext() { return mParser->getNext(); }
    
};

#endif
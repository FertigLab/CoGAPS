#include "../utils/GapsAssert.h"
#include "CharacterDelimitedParser.h"
#include "FileParser.h"
#include "MtxParser.h"

#include <string>

AbstractFileParser* AbstractFileParser::create(const std::string &path)
{
    switch (FileParser::fileType(path))
    {
        case GAPS_MTX: return new MtxParser(path);
        case GAPS_CSV: return new CharacterDelimitedParser(path, ',');
        case GAPS_TSV: return new CharacterDelimitedParser(path, '\t');
        case GAPS_GCT: return new CharacterDelimitedParser(path, '\t', true);
        default: GAPS_ERROR("Invalid file type\n");
    }
}

AbstractFileParser::AbstractFileParser() {}
AbstractFileParser::~AbstractFileParser() {}

std::vector<std::string> AbstractFileParser::rowNames() const
{
    return std::vector<std::string>();
}

std::vector<std::string> AbstractFileParser::colNames() const
{
    return std::vector<std::string>();
}

FileParser::FileParser(const std::string &path)
{
    mParser = AbstractFileParser::create(path);
}

FileParser::~FileParser()
{
    delete mParser;
}

std::vector<std::string> FileParser::rowNames() const
{
    return mParser->rowNames();
}

std::vector<std::string> FileParser::colNames() const
{
    return mParser->colNames();
}

unsigned FileParser::nRow() const
{
    return mParser->nRow();
}

unsigned FileParser::nCol() const
{
    return mParser->nCol();
}

bool FileParser::hasNext()
{
    return mParser->hasNext();
}

MatrixElement FileParser::getNext()
{
    return mParser->getNext();
}

GapsFileType FileParser::fileType(const std::string &path)
{
    std::size_t pos = path.find_last_of('.');
    std::string ext = path.substr(pos);

    if (ext.find('/') != std::string::npos) { return GAPS_INVALID_FILE_TYPE; }
    if (ext == ".mtx")  { return GAPS_MTX; }
    if (ext == ".csv")  { return GAPS_CSV; }
    if (ext == ".tsv")  { return GAPS_TSV; }
    if (ext == ".gct")  { return GAPS_GCT; }

    return GAPS_INVALID_FILE_TYPE;
}


#include "../utils/GapsAssert.h"
#include "CsvParser.h"
#include "FileParser.h"
#include "MtxParser.h"
#include "TsvParser.h"
#include "GctParser.h"
#include "Hdf5Parser.h"

#include <string>

AbstractFileParser* AbstractFileParser::create(const std::string &path)
{
    switch (FileParser::fileType(path))
    {
        case GAPS_MTX: return new MtxParser(path);
        case GAPS_CSV: return new CsvParser(path);
        case GAPS_TSV: return new TsvParser(path);
        case GAPS_GCT: return new GctParser(path);
        case GAPS_HDF5: return new Hdf5Parser(path);
        default: GAPS_ERROR("Invalid file type\n");
    }
}

AbstractFileParser::~AbstractFileParser() {}

FileParser::FileParser(const std::string &path)
{
    mParser = AbstractFileParser::create(path);
}

FileParser::~FileParser()
{
    delete mParser;
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
    if (ext == ".h5")   { return GAPS_HDF5; }
    if (ext == ".hdf5") { return GAPS_HDF5; }

    return GAPS_INVALID_FILE_TYPE;
}


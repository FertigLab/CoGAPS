#include "../GapsAssert.h"
#include "FileParser.h"
#include "MtxParser.h"
#include "CsvParser.h"
#include "TsvParser.h"

#include <string>

enum GapsFileType
{
    GAPS_MTX,
    GAPS_CSV,
    GAPS_TSV,
    GAPS_INVALID_FILE_TYPE
};

static GapsFileType fileType(const std::string &path)
{
    std::size_t pos = path.find_last_of(".");
    std::string ext = path.substr(pos);

    if (ext.find("/") != std::string::npos) { return GAPS_INVALID_FILE_TYPE; }
    else if (ext == ".mtx") { return GAPS_MTX; }
    else if (ext == ".csv") { return GAPS_CSV; }
    else if (ext == ".tsv") { return GAPS_TSV; }
    else { return GAPS_INVALID_FILE_TYPE; }
}

AbstractFileParser* AbstractFileParser::create(const std::string &path)
{
    switch (fileType(path))
    {
        case GAPS_MTX: return new MtxParser(path);
        case GAPS_CSV: return new CsvParser(path);
        case GAPS_TSV: return new TsvParser(path);
        default: GAPS_ERROR("Invalid file type\n");
    }
}

AbstractFileParser::~AbstractFileParser() {}
#include "GapsResult.h"

#include "GapsStatistics.h"
#include "file_parser/FileParser.h"

GapsResult::GapsResult(const GapsStatistics &stat)
    :
Amean(stat.Amean()), Asd(stat.Asd()), Pmean(stat.Pmean()),
Psd(stat.Psd()), meanChiSq(0.f), seed(0)
{}

void GapsResult::writeToFile(const std::string &fullPath)
{
    std::size_t pos = fullPath.find_last_of('.');
    std::string base = fullPath.substr(0, pos);

    switch (FileParser::fileType(fullPath))
    {
        case GAPS_CSV: return writeCsv(base);
        case GAPS_TSV: return writeTsv(base);
        case GAPS_GCT: return writeGct(base);
        default : GAPS_ERROR("Invalid File Type");
    }
}

void GapsResult::writeCsv(const std::string &path)
{
    unsigned nPatterns = Amean.nCol();
    std::string label("_" + gaps::to_string(nPatterns) + "_");
    FileParser::writeToCsv(path + label + "Amean.csv", Amean);
    FileParser::writeToCsv(path + label + "Pmean.csv", Pmean);
    FileParser::writeToCsv(path + label + "Asd.csv", Asd);
    FileParser::writeToCsv(path + label + "Psd.csv", Psd);
}

void GapsResult::writeTsv(const std::string &path)
{
    unsigned nPatterns = Amean.nCol();
    std::string label("_" + gaps::to_string(nPatterns) + "_");
    FileParser::writeToCsv(path + label + "Amean.tsv", Amean);
    FileParser::writeToCsv(path + label + "Pmean.tsv", Pmean);
    FileParser::writeToCsv(path + label + "Asd.tsv", Asd);
    FileParser::writeToCsv(path + label + "Psd.tsv", Psd);
}

void GapsResult::writeGct(const std::string &path)
{
    unsigned nPatterns = Amean.nCol();
    std::string label("_" + gaps::to_string(nPatterns) + "_");
    FileParser::writeToCsv(path + label + "Amean.gct", Amean);
    FileParser::writeToCsv(path + label + "Pmean.gct", Pmean);
    FileParser::writeToCsv(path + label + "Asd.gct", Asd);
    FileParser::writeToCsv(path + label + "Psd.gct", Psd);
}
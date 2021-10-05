#include "GapsResult.h"
#include "GapsStatistics.h"
#include "file_parser/FileParser.h"

template <class T>
static std::string to_string(T a)
{
    std::stringstream ss;
    ss << a;
    return ss.str();
}

GapsResult::GapsResult(){}

GapsResult::GapsResult(const GapsStatistics &stat)
    :
Amean(stat.Amean()), Asd(stat.Asd()), Pmean(stat.Pmean()),
Psd(stat.Psd()),
equilibrationSnapshotsA(stat.getEquilibrationSnapshots('A')),
equilibrationSnapshotsP(stat.getEquilibrationSnapshots('P')),
samplingSnapshotsA(stat.getSamplingSnapshots('A')),
samplingSnapshotsP(stat.getSamplingSnapshots('P')),
chisqHistory(stat.chisqHistory()), atomHistoryA(stat.atomHistory('A')),
atomHistoryP(stat.atomHistory('P')), seed(0), meanChiSq(0.f)
{}

void GapsResult::writeToFile(const std::string &path)
{
    unsigned nPatterns = Amean.nCol();
    std::string label("_" + to_string(nPatterns) + "_");
    FileParser::writeToCsv(path + label + "Amean.csv", Amean);
    FileParser::writeToCsv(path + label + "Pmean.csv", Pmean);
    FileParser::writeToCsv(path + label + "Asd.csv", Asd);
    FileParser::writeToCsv(path + label + "Psd.csv", Psd);
}

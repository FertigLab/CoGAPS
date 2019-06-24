#ifndef __COGAPS_FILE_PARSER_H__
#define __COGAPS_FILE_PARSER_H__

#include "MatrixElement.h"
#include "../utils/GapsAssert.h"

#include <fstream>
#include <vector>

enum GapsFileType
{
    GAPS_MTX,
    GAPS_CSV,
    GAPS_TSV,
    GAPS_GCT,
    GAPS_INVALID_FILE_TYPE
};

// file parser interface
class AbstractFileParser
{
public:
    static AbstractFileParser* create(const std::string &path);
    AbstractFileParser();
    virtual ~AbstractFileParser() = 0;
    virtual std::vector<std::string> rowNames() const;
    virtual std::vector<std::string> colNames() const;
    virtual unsigned nRow() const = 0;
    virtual unsigned nCol() const = 0;
    virtual bool hasNext() = 0;
    virtual MatrixElement getNext() = 0;
private:
    AbstractFileParser(const AbstractFileParser &p); // don't allow copies
    AbstractFileParser& operator=(const AbstractFileParser &p); // don't allow copies
};

// wrap the pointer to the parser implementation
class FileParser
{
public:
    explicit FileParser(const std::string &path);
    ~FileParser();
    std::vector<std::string> rowNames() const;
    std::vector<std::string> colNames() const;
    unsigned nRow() const;
    unsigned nCol() const;
    bool hasNext();
    MatrixElement getNext();
    static GapsFileType fileType(const std::string &path);
    template <class MatrixType>
    static void writeToCsv(const std::string &path, const MatrixType &mat);
private:
    FileParser(const FileParser &p); // don't allow copies
    FileParser& operator=(const FileParser &p); // don't allow copies

    AbstractFileParser *mParser;
};

template <class MatrixType>
void FileParser::writeToCsv(const std::string &path, const MatrixType &mat)
{
    if (FileParser::fileType(path) != GAPS_CSV)
    {
        GAPS_ERROR("output file must be a csv");
    }
    std::ofstream outputFile;
    outputFile.open(path.c_str());
    outputFile << "\"\"";

    // write column names
    for (unsigned i = 0; i < mat.nCol(); ++i)
    {
        outputFile << ",\"Col" << i << "\"";
    }
    outputFile << "\n";
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        // write row names
        outputFile << "\"Row" << i << "\"";
        
        // write data
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            outputFile << "," << mat(i,j);
        }
        outputFile << "\n";
    }
    outputFile.close();
}

#endif
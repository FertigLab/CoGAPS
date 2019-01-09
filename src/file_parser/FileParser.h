#ifndef __COGAPS_FILE_PARSER_H__
#define __COGAPS_FILE_PARSER_H__

#include "../utils/GapsAssert.h"
#include "MatrixElement.h"

#include <fstream>

enum GapsFileType
{
    GAPS_MTX,
    GAPS_CSV,
    GAPS_TSV,
    GAPS_GCT,
    GAPS_HDF5,
    GAPS_INVALID_FILE_TYPE
};

// file parser interface
class AbstractFileParser
{
public:

    AbstractFileParser() {}

    static AbstractFileParser* create(const std::string &path);

    virtual ~AbstractFileParser() = 0;

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

    unsigned nRow() const;
    unsigned nCol() const;

    bool hasNext();
    MatrixElement getNext();

    static GapsFileType fileType(const std::string &path);

    template <class MatrixType>
    static void writeToTsv(const std::string &path, const MatrixType &mat);

    template <class MatrixType>
    static void writeToCsv(const std::string &path, const MatrixType &mat);

    template <class MatrixType>
    static void writeToMtx(const std::string &path, const MatrixType &mat);

    template <class MatrixType>
    static void writeToGct(const std::string &path, const MatrixType &mat);

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    AbstractFileParser *mParser;

    FileParser(const FileParser &p); // don't allow copies
    FileParser& operator=(const FileParser &p); // don't allow copies
};

// temporary solution - should be moved into specific file parsers, ok for now
// since writing to file not exposed to user, only used for internal testing

template <class MatrixType>
void FileParser::writeToTsv(const std::string &path, const MatrixType &mat)
{
    GAPS_ASSERT(FileParser::fileType(path) == GAPS_TSV);

    std::ofstream outputFile;
    outputFile.open(path.c_str());
    outputFile << "\"\"";

    // write column names
    for (unsigned i = 0; i < mat.nCol(); ++i)
    {
        outputFile << "\t\"Col" << i << "\"";
    }
    outputFile << "\n";

    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        // write row names
        outputFile << "\"Row" << i << "\"";
        
        // write data
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            outputFile << "\t" << mat(i,j);
        }
        outputFile << "\n";
    }
    outputFile.close();
}

template <class MatrixType>
void FileParser::writeToCsv(const std::string &path, const MatrixType &mat)
{
    GAPS_ASSERT(FileParser::fileType(path) == GAPS_CSV);

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

template <class MatrixType>
void FileParser::writeToMtx(const std::string &path, const MatrixType &mat)
{
    GAPS_ASSERT(FileParser::fileType(path) == GAPS_MTX);

    std::ofstream outputFile;
    outputFile.open(path.c_str());
    outputFile << "%%\n";
    outputFile << mat.nRow() << " " << mat.nCol() << " " << mat.nRow() * mat.nCol();
    outputFile << "\n";
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            if (mat(i,j) > 0.f)
            {
                outputFile << i + 1 << " " << j + 1 << " " << mat(i,j) << "\n";
            }
        }
    }
    outputFile.close();
}

template <class MatrixType>
void FileParser::writeToGct(const std::string &path, const MatrixType &mat)
{
    GAPS_ASSERT(FileParser::fileType(path) == GAPS_GCT);

    std::ofstream outputFile;
    outputFile.open(path.c_str());
    outputFile << "#1.2\n";
    outputFile << mat.nRow() << "\t" << mat.nCol() << "\n";

    outputFile << "\"NAME\"\t\"Description\"";

    // write column names
    for (unsigned i = 0; i < mat.nCol(); ++i)
    {
        outputFile << "\t\"Col" << i << "\"";
    }
    outputFile << "\n";

    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        // write row names
        outputFile << "\"Row" << i << "\"\t\"BLANK\"\t";
        
        // write data
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            outputFile << "\t" << mat(i,j);
        }
        outputFile << "\n";
    }
    outputFile.close();
}

#endif
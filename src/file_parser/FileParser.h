#ifndef __COGAPS_FILE_PARSER_H__
#define __COGAPS_FILE_PARSER_H__

#include "MatrixElement.h"

enum GapsFileType
{
    GAPS_MTX,
    GAPS_CSV,
    GAPS_TSV,
    GAPS_INVALID_FILE_TYPE
};

// file parser interface
class AbstractFileParser
{
private:

    AbstractFileParser(const AbstractFileParser &p); // don't allow copies
    AbstractFileParser& operator=(const AbstractFileParser &p); // don't allow copies

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

    explicit FileParser(const std::string &path);
    ~FileParser();

    unsigned nRow();
    unsigned nCol();

    bool hasNext();
    MatrixElement getNext();

    static GapsFileType fileType(const std::string &path);

    template <class MatrixType>
    static void writeToTsv(const std::string &path, const MatrixType &mat);

    template <class MatrixType>
    static void writeToCsv(const std::string &path, const MatrixType &mat);

    template <class MatrixType>
    static void writeToMtx(const std::string &path, const MatrixType &mat);
};

// temporary solution - should be moved into specific file parsers, ok for now
// since writing to file not exposed to user, only used for internal testing

template <class MatrixType>
void FileParser::writeToTsv(const std::string &path, const MatrixType &mat)
{
    std::ofstream outputFile;
    outputFile.open(path.c_str());
    outputFile << "\"\"";

    // write column names
    for (unsigned i = 0; i < mat.nCol(); ++i)
    {
        outputFile << "\t\"Col" << i << "\"";
    }

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
    std::ofstream outputFile;
    outputFile.open(path.c_str());
    outputFile << "\"\"";

    // write column names
    for (unsigned i = 0; i < mat.nCol(); ++i)
    {
        outputFile << ",\"Col" << i << "\"";
    }

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
    std::ofstream outputFile;
    outputFile.open(path.c_str());
    outputFile << "%%\n";
    outputFile << mat.nRow() << " " << mat.nCol() << " " << mat.nRow() * mat.nCol();
    outputFile << "\n";
    for (unsigned j = 0; j < mat.nRow(); ++j)
    {
        for (unsigned i = 0; i < mat.nCol(); ++i)
        {
            if (mat(i,j) > 0.f)
            {
                outputFile << i + 1 << " " << j + 1 << " " << mat(i,j) << "\n";
            }
        }
    }
    outputFile.close();
}

#endif
#ifndef __COGAPS_HDF5_PARSER_H__
#define __COGAPS_HDF5_PARSER_H__

#include "FileParser.h"
#include "MatrixElement.h"

#include "H5Cpp.h"

#include <fstream>
#include <string>

class Hdf5Parser : public AbstractFileParser
{
public:

    explicit Hdf5Parser(const std::string &path);
    ~Hdf5Parser() { mFile.close(); }

    unsigned nRow() const { return mNumRows; }
    unsigned nCol() const { return mNumCols; }

    bool hasNext();
    MatrixElement getNext();

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    H5::H5File mFile;

    unsigned mNumRows;
    unsigned mNumCols;

    Hdf5Parser(const Hdf5Parser &p); // don't allow copies
    Hdf5Parser& operator=(const Hdf5Parser &p); // don't allow copies
};

#endif
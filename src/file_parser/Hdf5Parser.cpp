#include "Hdf5Parser.h"

#include "H5Cpp.h"

Hdf5Parser::Hdf5Parser(const std::string &path)
: mFile(path, H5F_ACC_RDONLY)
{
    
}

bool Hdf5Parser::hasNext()
{
    return false;
}

MatrixElement Hdf5Parser::getNext()
{
    return MatrixElement(0,0,0);
}
#ifndef __COGAPS_MATRIX_ELEMENT_H__
#define __COGAPS_MARRIX_ELEMENT_H__

#include <sstream>
#include <string>

struct MatrixElement
{
    unsigned row;
    unsigned col;
    float value;

    MatrixElement(unsigned r, unsigned c, float v)
        : row(r), col(c), value(v)
    {}

    MatrixElement(unsigned r, unsigned c, const std::string &v)
        : row(r), col(c), value(0.f)
    {
        std::stringstream ss(v);
        ss >> value;
    }
};

#endif
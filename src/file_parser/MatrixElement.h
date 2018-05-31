#ifndef __COGAPS_MATRIX_ELEMENT_H__
#define __COGAPS_MATRIX_ELEMENT_H__

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

    MatrixElement(unsigned r, unsigned c, const std::string &s)
        : row(r), col(c), value(0.f)
    {
        std::stringstream ss(s);
        ss >> value;
    }
};

#endif
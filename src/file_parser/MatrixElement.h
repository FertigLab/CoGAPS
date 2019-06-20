#ifndef __COGAPS_MATRIX_ELEMENT_H__
#define __COGAPS_MATRIX_ELEMENT_H__

#include <string>

struct MatrixElement
{
    unsigned row;
    unsigned col;
    float value;

    MatrixElement(unsigned r, unsigned c, float v);
    MatrixElement(unsigned r, unsigned c, const std::string &s);
};

#endif // __COGAPS_MATRIX_ELEMENT_H__

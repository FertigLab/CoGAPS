#ifndef __COGAPS_MATRIX_ELEMENT_H__
#define __COGAPS_MATRIX_ELEMENT_H__

#include <string>

struct MatrixElement
{
    MatrixElement(unsigned r, unsigned c, float v);
    MatrixElement(unsigned r, unsigned c, const std::string &s);

    unsigned row;
    unsigned col;
    float value;
};

#endif // __COGAPS_MATRIX_ELEMENT_H__

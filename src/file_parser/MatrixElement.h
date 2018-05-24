#ifndef __COGAPS_MATRIX_ELEMENT_H__
#define __COGAPS_MARRIX_ELEMENT_H__

struct MatrixElement
{
    unsigned row;
    unsigned col;
    float value;

    MatrixElement(unsigned r, unsigned c, float v)
        : row(r), col(c), value(v)
    {}
};

#endif
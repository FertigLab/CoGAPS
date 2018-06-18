#ifndef __COGAPS_MATRIX_ELEMENT_H__
#define __COGAPS_MATRIX_ELEMENT_H__

#include <sstream>
#include <string>

struct MatrixElement
{
    unsigned dim[2];
    float val;

    MatrixElement(unsigned r, unsigned c, float v)
        : val(v)
    {
        dim[0] = r;
        dim[1] = c;
    }

    MatrixElement(unsigned r, unsigned c, const std::string &s)
        :  val(0.f)
    {
        dim[0] = r;
        dim[1] = c;
        std::stringstream ss(s);
        ss >> val;
    }

    unsigned operator[](unsigned i)
    {
        return dim[i];
    }

    unsigned row() const
    {
        return dim[0];
    }

    unsigned col() const
    {
        return dim[1];
    }

    float value() const
    {
        return val;
    }
};

#endif

#ifndef __COGAPS_MATRIX_ELEMENT_H__
#define __COGAPS_MATRIX_ELEMENT_H__

#include "../utils/GapsAssert.h"
#include "../utils/GapsPrint.h"

#include <sstream>
#include <string>

inline bool containsNumber(const std::string &s)
{
    return !s.empty() && s.find_first_not_of("0123456789.-") == std::string::npos;
}

struct MatrixElement
{
    unsigned row;
    unsigned col;
    float value;

    MatrixElement(unsigned r, unsigned c, float v) // NOLINT
        : row(r), col(c), value(v)
    {}

    MatrixElement(unsigned r, unsigned c, const std::string &s) // NOLINT
        :  row(r), col(c), value(0.f)
    {
        std::stringstream ss(s);
        ss >> value;
        if (!containsNumber(s))
        {
            gaps_printf("\nError: Invalid entry found in input data: %s\n", s.c_str());
            gaps_stop();
        }
    }
};

#endif

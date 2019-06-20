#include "MatrixElement.h"

#include "../utils/GapsAssert.h"
#include "../utils/GapsPrint.h"

#include <sstream>
#include <string>

static bool containsNumber(const std::string &s)
{
    return !s.empty() && s.find_first_not_of("0123456789.-") == std::string::npos;
}

MatrixElement::MatrixElement(unsigned r, unsigned c, float v) // NOLINT
    :
row(r), col(c), value(v)
{}

MatrixElement::MatrixElement(unsigned r, unsigned c, const std::string &s) // NOLINT
    :
row(r), col(c), value(0.f)
{
    std::stringstream ss(s);
    ss >> value;
    if (!containsNumber(s))
    {
        gaps_printf("\nError: Invalid entry found in input data: %s\n", s.c_str());
        gaps_stop();
    }
}

#include "MatrixElement.h"

#include "../utils/GapsAssert.h"
#include "../utils/GapsPrint.h"

#include <sstream>
#include <string>
#include <cmath>

static bool isNumber(const std::string &s)
{
    return !s.empty() && s.find_first_not_of("0123456789.-") == std::string::npos;
}

static float processValue(const std::string &s)
{
    if (isNumber(s)) // most common case should be the fastest
    {
        std::stringstream ss(s);
        float tempValue;
        ss >> tempValue;
        return tempValue;
    }
    std::size_t pos = s.find("e"); // scientific notation only reason isNumber should fail
    if (pos == std::string::npos)
    {
        gaps_printf("\nError: Invalid entry found in input data: %s\n", s.c_str());
        gaps_stop();
    }
    else // handle scientific notation
    {
        std::string sBase = s.substr(0, pos);
        std::string sExp = s.substr(pos + 1);
        if (!isNumber(sBase) || !isNumber(sExp))
        {
            gaps_printf("\nError: Invalid entry found in input data: %s\n", s.c_str());
            gaps_stop();
        }
        std::stringstream ssExp(sExp);
        std::stringstream ssBase(sBase);
        float fBase;
        float fExp;
        ssBase >> fBase;
        ssExp >> fExp;
        return fBase * std::pow(10.f, fExp);
    }
}

MatrixElement::MatrixElement(unsigned r, unsigned c, const std::string &s) // NOLINT
    : row(r), col(c), value(processValue(s))
{}
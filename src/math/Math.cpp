#include "Math.h"

float gaps::min(float a, float b)
{
    return a < b ? a : b;
}

unsigned gaps::min(unsigned a, unsigned b)
{
    return a < b ? a : b;
}

uint64_t gaps::min(uint64_t a, uint64_t b)
{
    return a < b ? a : b;
}

float gaps::max(float a, float b)
{
    return a < b ? b : a;
}

unsigned gaps::max(unsigned a, unsigned b)
{
    return a < b ? b : a;
}

uint64_t gaps::max(uint64_t a, uint64_t b)
{
    return a < b ? b : a;
}
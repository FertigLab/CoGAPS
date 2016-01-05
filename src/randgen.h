#ifndef _RANDGEN_H_
#define _RANDGEN_H_

#include <iostream>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random.hpp>
#include <ctime>
#include <boost/random/uniform_01.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>

extern boost::mt19937 rng;

double randgen(char rand_type, double para1 = 0, double para2 = 1);

#endif

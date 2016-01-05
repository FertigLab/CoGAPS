#include <iostream>
#include <vector>
#include <algorithm>

#include "randgen.h"

using namespace std;

double randgen(char rand_type, double para1, double para2) {
    switch (rand_type) {
        case 'U': {
            boost::random::uniform_01<boost::mt19937 &> zeroone(rng);
            return zeroone();
        }

        case 'N': {
            double mean = para1;
            double std_dev = para2;
            boost::normal_distribution<> nd(mean, std_dev);
            boost::variate_generator<boost::mt19937 &, boost::normal_distribution<> > norm_rnd(rng, nd);
            return norm_rnd();
        }

        case 'P': {
            double lambda = (para1 == 0 ? para2 : para1);
            boost::poisson_distribution<> pd(lambda);
            boost::variate_generator<boost::mt19937 &, boost::poisson_distribution<> > poisson_rnd(rng, pd);
            return poisson_rnd();
        }

        case 'E': {
            double lambda = (para1 == 0 ? para2 : para1);
            boost::exponential_distribution<> expd(lambda);
            boost::variate_generator<boost::mt19937 &, boost::exponential_distribution<> > exp_rnd(rng, expd);
            return exp_rnd();
        }

        default:
            return -9999.0;
    }

    // EJF -- return dummy value to avoid warning
    return -9999.0;
}

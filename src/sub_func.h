#ifndef _SUB_FUNC_H
#define _SUB_FUNC_H

#include <iostream>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/gamma.hpp> // for gamma distribution
#include <boost/math/distributions/normal.hpp> // for normal distribution
#include <boost/math/distributions/exponential.hpp> // for exponential distribution
#include "randgen.h" // for generating uniformly random variates in [0,1].

namespace gaps{

class sub_func

{
 private:

 public:
  sub_func(){};
  ~sub_func(){};

  static double pexp(double p, double rate, bool boolpara1, bool boolpara2);

  static double qexp(double q, double rate, bool boolpara1, bool boolpara2);

  static double dgamma(double newMass, double shape, double scale, bool boolpara);

  static double pgamma(double p, double shape, double scale, bool boolpara1, bool boolpara2);

  static double qgamma(double q, double shape, double scale, bool boolpara1, bool boolpara2);

  static double dnorm(double u, double mean, double sd, bool unknown);

  static double qnorm(double u, double mean, double sd, double INF_Ref,double unknown);

  static double pnorm(double u, double mean, double sd, double INF_Ref,double unknown);

  //static double runif(double a, double b, double rng); // old form, changed into the next line.
  static double runif(double a, double b);

};

}

#endif

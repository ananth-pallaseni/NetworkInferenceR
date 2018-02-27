#ifndef __GAMMA_FIT_H
#define __GAMMA_FIT_H

#include <RcppArmadillo.h>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/distributions/gamma.hpp>
using namespace Rcpp;
using namespace std;

boost::math::gamma_distribution<> fit_gamma(int n, double sumx, double logsumx);
boost::math::gamma_distribution<> fit_gamma(const arma::mat &vals);

#endif

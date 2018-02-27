#ifndef __MEASURES_H
#define __MEASURES_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;

arma::mat joint_probabilities(const arma::mat &xvals, const arma::mat &yvals, int num_bins_x, int num_bins_y);

double mi(const arma::mat &joint, const arma::vec &xprobs, const arma::rowvec &yprobs);

double mi(const arma::mat &joint);

arma::mat specific_information_xy(const arma::mat &joint, const arma::vec &x, const arma::rowvec &y);

arma::mat specific_information_yx(const arma::mat &joint, const arma::rowvec &y, const arma::vec &x);

arma::mat specific_information(const arma::mat &joint, int dimsum);

double redundancy(const arma::mat target_probs, const arma::mat spec1, const arma::mat spec2);

#endif

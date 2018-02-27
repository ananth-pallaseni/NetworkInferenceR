#include "measures.h"

// Returns element-wise divide of a and b.
// If b is a column vector then divides columns of a by b.
// If b is a row vector then divides rows of a by b.
arma::mat mdiv(const arma::mat a, const arma::mat b) {
  if (size(a) == size(b)) {
    return a / b;
  }
  else if (a.n_rows == b.n_rows) {
    return a.each_col() / b;
  }
  else if (a.n_cols == b.n_cols) {
    return a.each_row() / b;
  }
  else {
    stop("Dimension Mismatch");
    return a / b;
  }
}

// Returns matrix multiply of a and b. Invariant to order of a and b
arma::mat mmul(const arma::mat a, const arma::mat b) {
  if (a.n_cols == b.n_rows) {
    return a * b;
  }
  else if (a.n_rows == b.n_cols) {
    return b * a;
  }
  else {
    stop("Dimension Mismatch");
  }
}

// Returns an num_bins_x by num_bins_y matrix of joint probabilities. xvals and yvals are binned values.
arma::mat joint_probabilities(const arma::mat &xvals, const arma::mat &yvals, int num_bins_x, int num_bins_y) {
  if (xvals.size() != yvals.size()) stop("Dimension mismatch : binned value arrays are of different size.");
  int n = xvals.size();

  arma::mat jp = arma::zeros(num_bins_x, num_bins_y);

  for (int i = 0; i < n; i++) {
    int x = xvals[i];
    int y = yvals[i];
    jp(x, y)++;
  }

  return jp / n;
}

// Returns the sum of finite values in the matrix
double finite_accu(const arma::mat &x) {
  double ret = 0;
  for (size_t i = 0; i < x.n_rows; i++) {
    for (size_t j = 0; j < x.n_cols; j++) {
      double val = x(i, j);
      ret += ( std::isnan(val) || std::isinf(val) ? 0 : val );
    }
  }
  return ret;
}

void remove_non_finite(arma::mat &x) {
  int nrow = x.n_rows;
  int ncol = x.n_cols;
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      double val = x(i, j);
      x(i, j) = std::isnan(val) || std::isinf(val) ? 0.0 : val;
    }
  }
}

// Assumes xprobs is a col vector and yprobs is a row vector
double mi(const arma::mat &joint, const arma::vec &xprobs, const arma::rowvec &yprobs) {
  arma::mat interim = joint % log2(joint / (xprobs * yprobs));
  return finite_accu(interim);
}

// Performs sums to get individual vectors
double mi(const arma::mat &joint) {
  return mi(joint, arma::sum(joint, 1), arma::sum(joint, 0));
}

arma::mat specific_information_xy(const arma::mat &joint, const arma::vec &x, const arma::rowvec &y) {
  arma::mat interim = mdiv(joint, y) % arma::log2( joint / (x * y) );
  remove_non_finite(interim);
  return arma::sum(interim, 0);
}

arma::mat specific_information_yx(const arma::mat &joint, const arma::rowvec &y, const arma::vec &x) {
  arma::mat interim = mdiv(joint, x) % arma::log2( joint / (x * y) );
  remove_non_finite(interim);
  return arma::sum(interim, 1);
}

arma::mat specific_information(const arma::mat &joint, int dimsum) {
  if (dimsum == 0) {
    return specific_information_xy(joint, arma::sum(joint, 1), arma::sum(joint, 0));
  }
  else {
    return specific_information_yx(joint, arma::sum(joint, 0), arma::sum(joint, 1));
  }
}

double redundancy(const arma::mat target_probs, const arma::mat spec1, const arma::mat spec2) {
  if (spec1.size() != spec2.size()) stop("Dimension mismatch : specific information vectors of different size");
  int n = spec1.size();
  arma::mat minvec = arma::zeros(n);
  for (int i = 0; i < n; i++) {
    minvec[i] = min(spec1[i], spec2[i]);
  }

  return arma::accu(target_probs % minvec);
}


//////////////////////// Test ////////////////////////
// void measures_test() {
//
//   arma::mat jp = arma::zeros(2, 3);
//   jp(0, 0) = 1; jp(0, 1) = 2; jp(0, 2) = 3;
//   jp(1, 0) = 4; jp(1, 1) = 5; jp(1, 2) = 6;
//
//   arma::vec x = arma::sum(jp, 1);
//
//   arma::rowvec y = arma::sum(jp, 0);
//
//   double mival = mi(jp); // should be about -92.029
//   Rcout << mival << " should equal about -92.029" << endl;
//
//   arma::mat targetspecval = arma::zeros(3);
//   targetspecval(0) = -4.36443; targetspecval(1) = -4.39232; targetspecval(2) = -4.38454;
//   arma::mat specval = specific_information(jp, 0);
//   for (int i = 0; i < 3; i++) {
//     Rcout << specval(i) << " should equal about " << targetspecval(i) << endl;
//   }
//
//   Rcout << specific_information(jp, 1) << endl;
// }

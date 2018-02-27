#ifndef __COMMON_H
#define __COMMON_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;

class Node {
public:
  string label;
  int number_of_bins;
  arma::mat binned_values;
  arma::mat probabilities;
  Node();
  Node(const String lbl, const NumericVector &data);
};


vector<Node> get_nodes(const DataFrame &df);

class NodePair {
  public:
  double mi;
  arma::mat spec;
};

#endif

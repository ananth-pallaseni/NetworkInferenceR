#include <RcppArmadillo.h>
using namespace Rcpp;

#include "common.h"
#include "measures.h"
#include "gamma_fit.h"

arma::mat joint_probabilities(const Node &n1, const Node &n2) {
  return joint_probabilities(n1.binned_values, n2.binned_values, n1.number_of_bins, n2.number_of_bins);
}

double mi(const Node &n1, const Node &n2) {
  arma::mat joint_12 = joint_probabilities(n1.binned_values, n2.binned_values, n1.number_of_bins, n2.number_of_bins);
  return mi(joint_12);
}

NumericMatrix get_mi(const vector<Node> &nodes) {
  int num_nodes = nodes.size();
  NumericMatrix mat(num_nodes, num_nodes);
  mat.fill(0);

  for (int i = 0; i < num_nodes; i++) {
    for (int j = i+1; j < num_nodes; j++) {
      double pair_mi = mi(nodes[i], nodes[j]);
      mat(i, j) = pair_mi;
      mat(j, i) = pair_mi;
    }
  }

  return mat;
}

NumericMatrix get_mi(const DataFrame &df) {
  vector<Node> nodes = get_nodes(df);
  return get_mi(nodes);
}

// Add puc contribution from xyz triplet to puc_scores
void increment_puc_scores(const Node &znode, const NodePair &xzpair, const NodePair &yzpair, int x, int y, int z, NumericMatrix &puc_scores) {
  double redun = redundancy(znode.probabilities, xzpair.spec, yzpair.spec);

  double xzval = (xzpair.mi - redun) / xzpair.mi;
  xzval = isfinite(xzval) ? xzval : 0;
  puc_scores(x, z) += xzval;
  puc_scores(z, x) += xzval;

  double yzval = (yzpair.mi - redun) / yzpair.mi;
  yzval = isfinite(yzval) ? yzval : 0;
  puc_scores(y, z) += yzval;
  puc_scores(z, y) += yzval;
}

// Poorly implemented matrix class because I can't deal with c++ object matrix conventions.
class NodePairMat {
  int numrows;
  int numcols;
  vector<NodePair> arr;

  public:
  NodePair* get(int i, int j) {
    return &arr[i * numcols + j];
  }
  void setmi(int i, int j, double mi) {
    arr[i * numcols + j].mi = mi;
  }
  void setspec(int i, int j, arma::mat spec) {
    arr[i * numcols + j].spec = spec;
  }
  NodePairMat(int rows, int cols) : numrows(rows), numcols(cols) {
    arr.resize(rows * cols);
  }
};

// Returns a matrix where element (i,j) is the puc between nodes i and j
NumericMatrix get_puc(const vector<Node> &nodes) {
  int n = nodes.size();

  NumericMatrix puc_scores(n, n);
  puc_scores.fill(0);


  // Populate mi and si for node pairs
  //arma::Mat<NodePair> node_pairs(n, n);
  NodePairMat node_pairs(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = i+1; j < n; j++) {
      arma::mat jp = joint_probabilities(nodes[i], nodes[j]);
      arma::vec p1 = arma::sum(jp, 1);
      arma::rowvec p2 = arma::sum(jp, 0);
      double pairmi = mi(jp, p1, p2);
      arma::mat spec12 = specific_information(jp, 0);
      arma::mat spec21 = specific_information(jp, 1);
      node_pairs.setmi(i, j, pairmi);
      node_pairs.setspec(i, j, spec12);
      node_pairs.setmi(j, i, pairmi);
      node_pairs.setspec(j, i, spec21);
    }
  }


  // For every triplet, populate puc_scores
  for (int i = 0; i < n; i++) {
    for (int j = i+1; j < n; j++) {
      for (int k = j+1; k < n; k++) {
        increment_puc_scores(nodes[k], *node_pairs.get(i, k), *node_pairs.get(j, k), i, j, k, puc_scores);
        increment_puc_scores(nodes[j], *node_pairs.get(i, j), *node_pairs.get(k, j), i, k, j, puc_scores);
        increment_puc_scores(nodes[i], *node_pairs.get(j, i), *node_pairs.get(k, i), j, k, i, puc_scores);
      }
    }
  }

  return puc_scores;
}

// Apply context to raw puc scores - transforms puc scores to pidc scores
NumericMatrix get_weights(const NumericMatrix puc_scores) {
  int ncol = puc_scores.ncol();
  int nrow = puc_scores.nrow();
  if (ncol != nrow) stop("Dimension mismatch : puc_score matrix has unequal numbers of rows and columns");

  // Precompute gamma distributions of variables
  vector< boost::math::gamma_distribution<> > gammas(ncol, boost::math::gamma_distribution<>(1));
  for (int i = 0; i < ncol; i++) {
    double sumx = 0;
    double logsumx = 0;
    for (int j = 0; j < ncol; j++) {
      if (j != i) {
        sumx += puc_scores(j, i);
        logsumx += log(puc_scores(j, i));
      }
    }

    gammas[i] = fit_gamma(nrow-1, sumx, logsumx);
  }


  NumericMatrix weights(ncol, ncol);
  weights.fill(0);
  for (int i = 0; i < ncol; i++) {
    for (int j = i+1; j < ncol; j++) {
      double score = puc_scores(i, j);
      double w = boost::math::cdf(gammas[i], score) + boost::math::cdf(gammas[j], score);
      weights(i, j) = w;
      weights(j, i) = w;
    }
  }

  return weights;
}

// Converts a vector of nodes and a matrix of infromation values into a dataframe representing a network.
DataFrame to_df(const vector<Node> nodes, const NumericMatrix infovals) {
  int n = nodes.size();
  int num_edges = n * (n-1) / 2;

  CharacterVector node1(num_edges);
  CharacterVector node2(num_edges);
  NumericVector weights(num_edges);

  int k = 0;
  for (int i = 0; i < n; i++) {
    Node n1 = nodes[i];
    for (int j = i+1; j < n; j++) {
      Node n2 = nodes[j];
      node1[k] = n1.label;
      node2[k] = n2.label;
      weights[k] = infovals(i, j);
      k++;
    }
  }

  DataFrame ret = DataFrame::create(Named("N1") = node1, Named("N2") = node2, Named("Weights") = weights);
  return ret;
}

// Returns an unordered dataframe with columns: node1, node2 and weight.
// [[Rcpp::export]]
DataFrame infer_mi_network(const DataFrame &df) {
  vector<Node> nodes = get_nodes(df);
  NumericMatrix mivals = get_mi(nodes);
  return to_df(nodes, mivals);
}

// Returns an unordered dataframe with columns: node1, node2 and weight.
// [[Rcpp::export]]
DataFrame infer_puc_network(const DataFrame &df) {
  vector<Node> nodes = get_nodes(df);
  NumericMatrix pucvals = get_puc(nodes);
  return to_df(nodes, pucvals);
}

// Returns an unordered dataframe with columns: node1, node2 and weight.
// [[Rcpp::export]]
DataFrame infer_pidc_network(const DataFrame &df) {
  vector<Node> nodes = get_nodes(df);
  NumericMatrix pidvals = get_weights(get_puc(nodes));
  return to_df(nodes, pidvals);
}

#include "common.h"

// Populates unique_data with the unique values that appear data and populates unique_counts with how often they appear.
void get_uniques(const NumericVector &data, NumericVector* unique_data, NumericVector* unique_counts) {
  int datasize = data.size();
  std::map<float, int> unique_map;
  for (int i = 0; i < datasize; i++) {
    float d = data[i];
    if (unique_map.find(d) == unique_map.end()) unique_map[d] = 1;
    else unique_map[d]++;
  }

  int n = unique_map.size();

  NumericVector ret_data(n);
  std::map<float, int>::iterator it = unique_map.begin();
  for (int i = 0; i < n; i++) {
    ret_data[i] = it->first;
    it++;
  }
  std::sort(ret_data.begin(), ret_data.end());

  NumericVector ret_counts(n);
  for (int i = 0; i < n; i++) ret_counts[i] = unique_map[ret_data[i]];

  *unique_data = ret_data;
  *unique_counts = ret_counts;
}

// Returns index of maximum value
int indmax(const NumericVector &invec) {
  int imax = 0;
  float ival = invec[0];
  for (int i = 1; i < invec.size(); i++) {
    if (invec[i] > ival) {
      imax = i;
      ival = invec[i];
    }
  }
  return imax;
}

// Returns edges of discretized bins
NumericVector bayesian_blocks_discretizer(const NumericVector &data) {
  NumericVector unique_data, unique_counts;
  get_uniques(data, &unique_data, &unique_counts);
  int n = unique_data.size();

  NumericVector edges(n+1);
  edges[0] = unique_data[0];
  for (int i = 0; i < n-1; i++) {
    edges[i+1] = 0.5 * (unique_data[i] + unique_data[i+1]);
  }
  edges[n] = unique_data[n-1];

  NumericVector block_length = unique_data[n-1] - edges; // Array ops are vectorized using rcpp

  NumericVector count_vec(n);
  NumericVector best(n);
  NumericVector last(n);
  for (int i = 0; i < n; i++) {
    count_vec[i] = 0;
    best[i] = 0;
    last[i] = 0;
  }

  for (int k = 0; k < n; k++) {
    NumericVector width(k+1);
    for (int i = 0; i < k+1; i++) width[i] = block_length[i] - block_length[k+1];
    for (int i = 0; i < k+1; i++) count_vec[i] += unique_counts[k];

    // Fitness function
    NumericVector fit_vec(k+1);
    for (int i = 0; i < k+1; i++) {
      fit_vec[i] = count_vec[i] * log(count_vec[i] / width[i]);
      fit_vec[i] -= 4 - log(73.53 * 0.05 * pow(k+1, -0.478));
    }
    for (int i = 1; i < k+1; i++) fit_vec[i] += best[i-1];

    int imax = indmax(fit_vec);
    last[k] = imax;
    best[k] = fit_vec[imax];
  }

  NumericVector change_points(n);
  for (int i = 0; i < n; i++) change_points[i] = 0;

  int i_cp = n;
  int ind = n;
  while (i_cp >= 0) {
    i_cp--;
    change_points[i_cp] = ind;
    if (ind == 0) break;
    ind = last[ind-1];
  }

  NumericVector ret(n-i_cp);
  int j = 0;
  for (int i = i_cp; i < n; i++) {
    int change_ind = change_points[i];
    ret[j] = edges[change_ind];
    j++;
  }

  return ret;
}

// Maps a value to the bin it falls into. Uses binary search.
int value_to_bin(double val, const NumericVector &bin_edges) {
  int start = 1;
  int fin = bin_edges.size() - 1;
  int mid = (fin + start) / 2;

  while (start < fin) {
    mid = (start + fin) / 2;
    if (val <= bin_edges[mid]) {
      fin = mid;
    }
    else {
      start = mid + 1;
    }
  }
  return fin - 1;
}

void bin_data(const NumericVector &data, const NumericVector &bin_edges, arma::mat &binned_ids) {
  int ndata = data.size();

  // Binary search
  for (int i = 0; i < ndata; i++) {
    binned_ids[i] = value_to_bin(data[i], bin_edges);
  }
}

// Fills out bin_ids vector and returns number of bins
int get_bin_ids_bayesian_blocks(const NumericVector &data, arma::mat &bin_ids) {
  NumericVector bin_edges = bayesian_blocks_discretizer(data);
  bin_data(data, bin_edges, bin_ids);
  return bin_edges.size() - 1;
}


// Takes in a list of bin ids and returns the probability of each bin
void get_probabilities_max_likelihood(const arma::mat &bin_ids, int num_bins, arma::mat &probabilities) {
  for (int i = 0; i < num_bins; i++) probabilities[i] = 0;

  int countsum = 0;
  int n = bin_ids.size();
  for (int i = 0; i < n; i++) {
    probabilities[ bin_ids[i] ] += 1;
    countsum++;
  }

  for (int i = 0; i < num_bins; i++) probabilities[i] /= countsum;
}

// Default
Node::Node() {};

Node::Node(const String lbl, const NumericVector &data) {
  label = lbl;

  arma::mat bvals = arma::zeros(data.size());
  number_of_bins = get_bin_ids_bayesian_blocks(data, bvals);
  binned_values = bvals;

  arma::mat probs = arma::zeros(number_of_bins);
  get_probabilities_max_likelihood(binned_values, number_of_bins, probs);
  probabilities = probs;
}

vector<Node> get_nodes(const DataFrame &df) {
  int n = df.length();
  vector<Node> nodes(n);

  CharacterVector node_names = df.names();
  for (int i = 0; i < n; i++) {
    String cname = node_names[i];
    NumericVector col = df[cname];
    nodes[i] = Node(cname, col);
  }


  return nodes;
}


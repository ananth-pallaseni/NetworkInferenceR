\name{infer_network}
\alias{infer_network}
\title{Infer a network}
\usage{
  infer_network(filepath, method = "PIDC")
}
\arguments{
  \item{filepath}{Path to the data file}
  \item{method}{Algorithm for inferring the network}
}
\description{
  Infer a network of statistically dependent variables, given a dataset containing measurements for a set of variables over a set of conditions.
}
\examples{
  ## Use the default PIDC method
  infer_network("data.tsv")

  ## Use proportional unique contribution
  infer_network("data.tsv", method = "PUC")

  ## Use mutual information
  infer_network("data.tsv", method = "MI")
}
